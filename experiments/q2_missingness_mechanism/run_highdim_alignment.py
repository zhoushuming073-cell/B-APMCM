import gzip
from pathlib import Path

source = Path(__file__).with_suffix(".py.gz")
namespace = {"__name__": "__main__", "__file__": str(Path(__file__).resolve())}
exec(compile(gzip.decompress(source.read_bytes()), str(source), "exec"), namespace)
