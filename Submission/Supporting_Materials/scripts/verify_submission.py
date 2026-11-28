"""Verify the clean supporting-materials staging directory before ZIP creation.

This is a packaging check only; it does not execute or refit research models.
Run before preparing the separately supplied official raw data.
"""
from pathlib import Path
import csv
import hashlib
import io
import json
import zipfile

ROOT = Path(__file__).resolve().parents[1]
Q23 = ROOT / "experiments/q2_final_validation"


def git_blob_sha(path):
    raw = path.read_bytes()
    return hashlib.sha1(f"blob {len(raw)}\0".encode("ascii") + raw).hexdigest()


def require(condition, message):
    if not condition:
        raise ValueError(message)


def main():
    with (ROOT / "SOURCE_MANIFEST.tsv").open(encoding="utf-8", newline="") as f:
        rows = list(csv.DictReader(f, delimiter="\t"))
    require(bool(rows), "Empty source manifest")
    require(len({r["source_path"] for r in rows}) == len(rows), "Duplicate manifest paths")
    for row in rows:
        path = ROOT / row["source_path"]
        require(path.is_file(), f"Missing source copy: {path.relative_to(ROOT)}")
        require(path.stat().st_size == int(row["bytes"]), f"Size mismatch: {path.relative_to(ROOT)}")
        require(git_blob_sha(path) == row["git_blob_sha"], f"Git blob mismatch: {path.relative_to(ROOT)}")

    manifest = json.loads((Q23 / "source_manifest.json").read_text(encoding="utf-8"))
    parts = [Q23 / "archive_parts" / f"reproduction.zip.part{i:02d}" for i in range(1, 11)]
    raw = b"".join(path.read_bytes() for path in parts)
    require(len(raw) == manifest["source_zip_size"], "Q2.3 archive size mismatch")
    require(hashlib.sha256(raw).hexdigest() == manifest["source_zip_sha256"], "Q2.3 archive SHA-256 mismatch")
    with zipfile.ZipFile(io.BytesIO(raw)) as z:
        expected = manifest["archive_members"]
        require(len(z.namelist()) == len(expected) and set(z.namelist()) == set(expected), "Q2.3 archive member mismatch")
        require("inputs/process.csv" not in expected, "Raw-equivalent process copy included")
        for name, sha in expected.items():
            require(not Path(name).is_absolute() and ".." not in Path(name).parts, f"Unsafe archive member: {name}")
            require(not name.lower().endswith((".xlsb", ".npz", "data.zip")), f"Prohibited input in archive: {name}")
            require(hashlib.sha256(z.read(name)).hexdigest() == sha, f"Archive member mismatch: {name}")

    for path in ROOT.rglob("*"):
        if not path.is_file():
            continue
        rel = path.relative_to(ROOT).as_posix().lower()
        require(not rel.startswith(("data/raw/", "data/interim/extracted/")), f"Official or extracted data in staging: {rel}")
        require(not rel.endswith((".xlsb", ".npz", "/data.zip")) and rel != "data.zip", f"Excluded input file: {rel}")
        with path.open("rb") as f:
            prefix = f.read(64)
        require(not prefix.startswith(b"version https://git-lfs.github.com/spec/v1"), f"Git LFS pointer: {rel}")

    print(f"Verified {len(rows)} byte-identical source copies.")
    print(f"Verified Q2.3 archive and {len(expected)} retained members.")
    print("No direct official archive, workbooks, extracted raw caches or Git LFS pointers found.")
    print("Packaging checks only: full model execution and final paper/figure review remain separate.")


if __name__ == "__main__":
    main()
