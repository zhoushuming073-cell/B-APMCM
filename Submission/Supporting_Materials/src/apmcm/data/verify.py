"""Verify canonical raw files and their official ZIP provenance."""
import hashlib
from pathlib import Path
from zipfile import ZipFile

from apmcm.common.paths import ARCHIVE, DATA, RAW, ROOT, SOURCE_SHA256, config


def digest(path):
    with path.open("rb") as handle:
        return hashlib.file_digest(handle, "sha256").hexdigest()


def verify():
    known = dict(
        (name, sha) for sha, name in (
            line.split(maxsplit=1)
            for line in (DATA / "metadata/checksums.sha256").read_text(encoding="utf-8").splitlines()
        )
    )
    settings = config("data.yaml")
    names = {path.name for path in RAW.glob("*.xlsb")}
    if len(names) != settings["expected_heats"] or len(known) != len(names) + 1:
        raise ValueError("Official workbook set or checksum manifest is incomplete")
    raw_hashes = {name: digest(RAW / name) for name in names}
    if digest(ARCHIVE) != SOURCE_SHA256 or known.get(ARCHIVE.relative_to(ROOT).as_posix()) != SOURCE_SHA256:
        raise ValueError("Official Data.zip hash differs. If you cloned the repo, run git lfs pull.")
    if any(known.get((RAW / name).relative_to(ROOT).as_posix()) != sha for name, sha in raw_hashes.items()):
        raise ValueError("An extracted official workbook does not match its recorded hash")
    with ZipFile(ARCHIVE) as archive:
        entries = [name for name in archive.namelist() if name.lower().endswith(".xlsb")]
        if len(entries) != len(names) or {Path(name).name for name in entries} != names:
            raise ValueError("Archive and standalone workbook names differ")
        for entry in entries:
            with archive.open(entry) as handle:
                if hashlib.file_digest(handle, "sha256").hexdigest() != raw_hashes[Path(entry).name]:
                    raise ValueError(f"Archive workbook differs: {entry}")
    return len(names)
