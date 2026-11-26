"""Check official files. Optionally rebuild the extracted spectra and audit tables."""
import argparse

from apmcm.data.verify import verify


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--rebuild-cache", action="store_true")
    args = parser.parse_args()
    print(f"Verified {verify()} workbooks and Data.zip")
    if args.rebuild_cache:
        from apmcm.data.loader import run
        run()


if __name__ == "__main__":
    main()
