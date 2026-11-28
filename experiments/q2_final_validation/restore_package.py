#!/usr/bin/env python3
"""Restore the exact archived Q2.3 source and frozen reference inputs."""
from __future__ import annotations

import argparse
import hashlib
import io
import json
from pathlib import Path
import zipfile

ROOT = Path(__file__).resolve().parent


def restore(out: Path) -> None:
    manifest = json.loads((ROOT / 'source_manifest.json').read_text(encoding='utf-8'))
    paths = [ROOT / 'archive_parts' / f'reproduction.zip.part{i:02d}' for i in range(1, 11)]
    raw = b''.join(p.read_bytes() for p in paths)
    if len(raw) != manifest['source_zip_size'] or hashlib.sha256(raw).hexdigest() != manifest['source_zip_sha256']:
        raise ValueError('Reproduction archive size or SHA-256 mismatch')
    with zipfile.ZipFile(io.BytesIO(raw)) as archive:
        expected = manifest['archive_members']
        if set(archive.namelist()) != set(expected) or len(archive.namelist()) != len(expected):
            raise ValueError('Archive member list mismatch')
        for name, digest in expected.items():
            path = Path(name)
            if path.is_absolute() or '..' in path.parts:
                raise ValueError(f'Unsafe archive path: {name}')
            if hashlib.sha256(archive.read(name)).hexdigest() != digest:
                raise ValueError(f'Archive member SHA-256 mismatch: {name}')
        if out.exists() and (not out.is_dir() or any(out.iterdir())):
            raise ValueError(f'Output directory must be empty: {out}')
        out.mkdir(parents=True, exist_ok=True)
        archive.extractall(out)
    print(f'Verified and restored {len(expected)} original files to {out}')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out', type=Path, required=True)
    restore(parser.parse_args().out.resolve())
