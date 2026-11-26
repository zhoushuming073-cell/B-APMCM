"""Extract only the Q2 numeric columns from the supplied single-sheet XLSB archives."""
import hashlib
import io
import re
import struct
import sys
import zipfile
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]


def records(data):
    i = 0
    while i < len(data):
        tag = data[i]
        i += 1
        if tag & 128:
            tag = (tag & 127) + (data[i] << 7)
            i += 1
        size = 0
        for shift in (0, 7, 14, 21):
            byte = data[i]
            i += 1
            size += (byte & 127) << shift
            if byte < 128:
                break
        yield tag, i
        i += size


def rk(data, pos):
    bits, = struct.unpack_from('<I', data, pos)
    if bits & 2:
        value = struct.unpack_from('<i', data, pos)[0] >> 2
    else:
        value, = struct.unpack('<d', struct.pack('<II', 0, bits & ~3))
    return value / 100 if bits & 1 else value


def extract(path):
    dest = ROOT / 'data'
    dest.mkdir(exist_ok=True)
    blocks, manifest = [], []
    with zipfile.ZipFile(path) as batch:
        for name in sorted(batch.namelist()):
            match = re.search(r'Attachment (\d+)\. Heat No\. (\w+)\.xlsb$', name)
            if not match:
                continue
            number, heat = match.groups()
            payload = batch.read(name)
            with zipfile.ZipFile(io.BytesIO(payload)) as book:
                data = book.read('xl/worksheets/sheet1.bin')
                strings = book.read('xl/sharedStrings.bin')
            labels = [strings[p+5:p+5+2*struct.unpack_from('<I',strings,p+1)[0]].decode('utf-16-le')
                      for tag,p in records(strings) if tag == 19]
            cells, headers = {}, {}
            row = -1
            for tag, pos in records(data):
                if tag == 0:
                    row, = struct.unpack_from('<I', data, pos)
                elif tag in (2, 5, 7):
                    col, = struct.unpack_from('<I', data, pos)
                    if col not in (0, 2048, 2049, 2050):
                        continue
                    if row == 0 and tag == 7:
                        headers[col] = labels[struct.unpack_from('<I',data,pos+8)[0]]
                    if row > 0 and tag in (2, 5):
                        cells[row, col] = rk(data,pos+8) if tag == 2 else struct.unpack_from('<d',data,pos+8)[0]
            columns = [0, 2048, 2049, 2050]
            arr = np.array([[cells[i,c] for c in columns] for i in range(1,row+1)])
            assert np.array_equal(arr[:,0],np.arange(1,len(arr)+1)) and np.isfinite(arr).all()
            assert 'CO2' in headers[2049].upper().replace('₂','2') and 'CO' in headers[2050].upper()
            part = pd.DataFrame(arr,columns=['row','O','CO2','CO']).assign(attachment=int(number),heat=heat)
            blocks.append(part)
            manifest.append(dict(attachment=int(number),heat=heat,rows=len(arr),source=name,
                                 headers=' | '.join(headers[c] for c in columns),sha256=hashlib.sha256(payload).hexdigest()))
    assert len(blocks) == 10 and sum(map(len,blocks)) == 13589
    pd.concat(blocks).sort_values(['attachment','row']).to_csv(dest/'process.csv',index=False,float_format='%.15g')
    pd.DataFrame(manifest).sort_values('attachment').to_csv(dest/'sources.csv',index=False)
    (dest/'archive_sha256.txt').write_text(hashlib.sha256(Path(path).read_bytes()).hexdigest()+'  Data.zip\n')


if __name__ == '__main__':
    extract(sys.argv[1])
