#!/usr/bin/env python3

from __future__ import annotations

import argparse
import gzip
import sys
from pathlib import Path


def open_maybe_gz(path: Path):
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8", errors="replace")
    return open(path, encoding="utf-8", errors="replace")


def fasta_headers(path: Path) -> list[str]:
    headers = []
    with open_maybe_gz(path) as fh:
        for line in fh:
            if line.startswith(">"):
                headers.append(line[1:].split()[0].strip() if line[1:].strip() else "")
    return headers


def base_accession(accession: str) -> str:
    return accession.rsplit(".", 1)[0] if "." in accession else accession


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--fasta", required=True, type=Path)
    ap.add_argument("--exclude-list", required=True, type=Path)
    ap.add_argument("--label", default="", help="cosmetic label for messages, e.g. the ablation level")
    ap.add_argument("--pass-marker", required=True, help="written on success, so this gates downstream steps")
    args = ap.parse_args()

    exclude_raw = [line.strip() for line in args.exclude_list.read_text().splitlines() if line.strip()]
    if not exclude_raw:
        Path(args.pass_marker).write_text(f"pass: exclude list is empty ({args.label})\n")
        print(f"verify_exclusion: exclude list empty, nothing to check ({args.label})", file=sys.stderr)
        return 0

    exclude_exact = set(exclude_raw)
    exclude_base = {base_accession(a) for a in exclude_raw}

    headers = fasta_headers(args.fasta)
    if not headers:
        sys.exit(f"error: {args.fasta} contains no FASTA headers -- refusing to call an empty file verified")

    violations = sorted(
        {h for h in headers if h in exclude_exact or base_accession(h) in exclude_base}
    )

    if violations:
        sys.exit(
            f"error: {args.fasta} ({args.label}) contains {len(violations)} sequence(s) that should have "
            f"been excluded by ablation: {', '.join(violations[:10])}"
            + (f" (+{len(violations) - 10} more)" if len(violations) > 10 else "")
            + ". The ablation did not take effect; this build must not proceed."
        )

    Path(args.pass_marker).write_text(
        f"pass: {len(headers)} sequences checked against {len(exclude_raw)} excluded accessions, "
        f"0 violations ({args.label})\n"
    )
    print(f"verify_exclusion: {len(headers)} sequences OK, 0 of {len(exclude_raw)} excluded accessions present ({args.label})", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
