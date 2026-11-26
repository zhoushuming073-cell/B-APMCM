from pathlib import Path
import re
import logging
from collections import Counter
import tempfile
from python_calamine import CalamineWorkbook
from zipfile import ZipFile
import numpy as np
import pandas as pd

from apmcm.common.paths import ARCHIVE, DATA, DT, PROC, RAW, SOURCE_SHA256
from .verify import digest

CURRENT_DATA = DATA / "metadata"

def run(source=None):
    archive = Path(source) if source else ARCHIVE
    if not archive.is_file() or digest(archive) != SOURCE_SHA256:
        raise FileNotFoundError("The verified corrected Data.zip is missing from the source directory.")
    PROC.mkdir(parents=True, exist_ok=True)
    CURRENT_DATA.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory() as tmp:
        with ZipFile(archive) as z:
            names=z.namelist()
            if any(Path(n).is_absolute() or ".." in Path(n).parts for n in names):
                raise ValueError("Unsafe archive path")
            z.extractall(tmp)
        files=sorted(Path(tmp).rglob("*.xlsb"))
        ids=[int(re.search(r"Attachment (\d+)",p.name)[1]) for p in files]
        if sorted(ids)!=list(range(1,11)): raise ValueError(f"Incomplete attachment set: {ids}")
        audits=[]; stats=[]; manifests=[]
        for i,p in sorted(zip(ids,files)):
            dest=RAW/p.name
            sha=digest(p)
            if not dest.is_file(): raise FileNotFoundError(f"Official workbook is missing: {dest}")
            if digest(dest)!=sha: raise ValueError(f"Original differs: {dest}")
            book=CalamineWorkbook.from_path(str(p))
            sheet=book.sheet_names[0]
            records=book.get_sheet_by_index(0).to_python()
            header=records.pop(0)
            types=Counter(type(v).__name__ for row in records for v in row)
            frame=pd.DataFrame(records)
            del records
            num=frame.apply(pd.to_numeric,errors="coerce")
            a=num.to_numpy(dtype=float)
            seq=a[:,0]; wl=np.asarray(header[1:-3],float); x=a[:,1:-3]; other=a[:,-3:]
            if len(book.sheet_names)!=1 or len(wl)!=2047 or not np.all(np.diff(wl)>0):
                raise ValueError(f"Unexpected layout {p.name}")
            if not np.all(np.diff(seq)==1): raise ValueError(f"Sequence discontinuity {p.name}")
            t=(seq-seq[0])*DT
            manifests.append(dict(attachment=i,heat=sheet,file=p.name,sha256=sha,fields=" | ".join(map(str,[header[0],*header[-3:]]))))
            audits.append(dict(attachment=i,heat=sheet,rows=len(a),cols=a.shape[1],sheet=sheet,
                cell_types=str(dict(types)),sequence_start=seq[0],sequence_end=seq[-1],duration_s=t[-1],
                wavelength_min=wl.min(),wavelength_max=wl.max(),wavelength_count=len(wl),
                non_numeric=int((num.isna() & frame.notna() & frame.ne("")).sum().sum()),
                missing_spectrum=int((~np.isfinite(x)).sum()),missing_process=int((~np.isfinite(other)).sum()),
                duplicate_sequence=len(seq)-len(np.unique(seq)),duplicate_measurement_rows=int(num.iloc[:,1:].duplicated().sum()),
                negative_spectrum=int((x<0).sum()),oxygen_drops=int((np.diff(other[:,0])<0).sum()),
                oxygen_repeat=int((np.diff(other[:,0])==0).sum()),negative_gas=int((other[:,1:]<0).sum())))
            for field,values in zip(["intensity","oxygen","CO2","CO"],[x,*other.T]):
                vals=np.nanpercentile(values,[0,5,50,95,100])
                stats.append(dict(heat=sheet,field=field,**dict(zip(["min","p05","median","p95","max"],vals)),mean=np.nanmean(values),sd=np.nanstd(values)))
            np.savez(PROC/f"raw_{sheet}.npz",t=t,sequence=seq,wl=wl,x=x.astype("f4"),oxygen=other[:,0],gas=other[:,1:])
            logging.info("Audited %s: %s",sheet,a.shape)
        for name, rows in [("audit", audits), ("distributions", stats), ("manifest", manifests)]:
            pd.DataFrame(rows).to_csv(CURRENT_DATA / (name + ".csv"), index=False)

