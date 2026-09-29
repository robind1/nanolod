#!/usr/bin/env python3

from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--accession2taxid", required=True, type=Path)
    ap.add_argument("--prefix", required=True)
    args = ap.parse_args()

    with open(args.accession2taxid, newline="", encoding="utf-8") as fh:
        sample = fh.read(4096)
        fh.seek(0)
        dialect = csv.Sniffer().sniff(sample, delimiters="\t,") if sample.strip() else csv.excel_tab
        reader = csv.DictReader(fh, dialect=dialect)
        if reader.fieldnames is None:
            sys.exit(f"error: {args.accession2taxid} has no header row")
        lower_to_actual = {f.lower(): f for f in reader.fieldnames}
        acc_col = lower_to_actual.get("accession.version") or lower_to_actual.get("accession")
        taxid_col = lower_to_actual.get("taxid")
        if acc_col is None or taxid_col is None:
            sys.exit(f"error: {args.accession2taxid} must have an accession and a taxid column; found {reader.fieldnames}")

        n = 0
        with open(f"{args.prefix}.seqid2taxid.map", "w", encoding="utf-8") as out:
            for row in reader:
                acc, taxid = row[acc_col].strip(), row[taxid_col].strip()
                if acc and taxid:
                    out.write(f"{acc}\t{taxid}\n")
                    n += 1

    print(f"make_seqid2taxid_map: wrote {n} mappings to {args.prefix}.seqid2taxid.map", file=sys.stderr)
    if n == 0:
        sys.exit("error: zero mappings written -- refusing to hand Kraken2 an empty seqid2taxid.map")
    return 0


if __name__ == "__main__":
    sys.exit(main())
