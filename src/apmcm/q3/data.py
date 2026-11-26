"""Extract Q3 data directly from the verified official ZIP without modifying data/."""
from pathlib import Path
from tempfile import TemporaryDirectory
from zipfile import ZipFile
import hashlib
import re
import numpy as np
from pyxlsb import open_workbook

from apmcm.common.paths import SOURCE_SHA256


def load_archive(archive):
    archive=Path(archive)
    with archive.open('rb') as f:
        sha=hashlib.file_digest(f,'sha256').hexdigest()
    if sha!=SOURCE_SHA256:
        raise ValueError('Q3 requires the corrected official Data.zip; run git lfs pull.')
    raw={}; ids=[]
    with ZipFile(archive) as z, TemporaryDirectory() as tmp:
        for name in sorted(z.namelist()):
            if not name.lower().endswith('.xlsb'):
                continue
            match=re.search(r'Attachment (\d+)\. Heat No\. (\w+)\.xlsb$',name)
            if not match:
                raise ValueError(f'Unexpected workbook name: {name}')
            number,heat=int(match[1]),match[2]; ids.append(number)
            dest=Path(tmp)/(heat+'.xlsb');dest.write_bytes(z.read(name))
            with open_workbook(dest) as book:
                with book.get_sheet(1) as sheet:
                    rows=sheet.rows();header=[c.v for c in next(rows)]
                    a=np.array([[c.v for c in row] for row in rows],dtype=np.float64)
            wl=np.array(header[1:-3],float)
            if a.shape[1]!=2051 or len(wl)!=2047 or not np.isfinite(a).all() or not np.all(np.diff(wl)>0) or not np.all(np.diff(a[:,0])==1):
                raise ValueError(f'Invalid official data layout: {heat}')
            if np.any(a[:,1:]<0) or np.any(np.ptp(a[:,-2:],axis=0)==0) or heat in raw:
                raise ValueError(f'Invalid values or duplicate heat: {heat}')
            raw[heat]=dict(wl=wl,x=a[:,1:-3].astype('f4'),oxygen=a[:,-3],gas=a[:,-2:],sequence=a[:,0])
            dest.unlink()
    if sorted(ids)!=list(range(1,11)) or len(raw)!=10 or sum(len(v['gas']) for v in raw.values())!=13589:
        raise ValueError('Official ten-heat coverage check failed')
    np.testing.assert_allclose(raw['6B03205']['gas'].max(0),[21.74,59.04])
    return raw
