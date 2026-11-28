"""CLI entry point for the Q3 physical-chemistry measurement audit."""
from pathlib import Path
import gzip

_source = Path(__file__).with_name("measurement_audit_source.py.gz")
_code = gzip.decompress(_source.read_bytes())
exec(compile(_code, str(_source), "exec"), globals())
