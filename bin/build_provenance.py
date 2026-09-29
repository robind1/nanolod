#!/usr/bin/env python3

from __future__ import annotations

import argparse
import hashlib
import json
import sys
from pathlib import Path


def sha256_of_paths(paths: list[Path]) -> str:

    h = hashlib.sha256()
    for p in sorted(paths, key=lambda p: p.name):
        with open(p, "rb") as fh:
            for chunk in iter(lambda: fh.read(1024 * 1024), b""):
                h.update(chunk)
    return h.hexdigest()


def load_json(path: Path | None) -> dict:
    if path is None:
        return {}
    return json.loads(path.read_text())


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument(
        "--ablation-json",
        required=True,
        type=Path,
        help="ablation resolution over the PROTEIN accession space (paired with the identity measurement)",
    )
    ap.add_argument(
        "--ablation-json-nt",
        type=Path,
        default=None,
        help=(
            "ablation resolution over the NUCLEOTIDE accession space."
        ),
    )
    ap.add_argument("--kraken2-identity-json", type=Path, default=None)
    ap.add_argument("--diamond-identity-json", type=Path, default=None)
    ap.add_argument("--kraken2-db-files", nargs="+", type=Path, default=[])
    ap.add_argument("--diamond-db-files", nargs="+", type=Path, default=[])
    ap.add_argument("--kraken2-version", default="unknown")
    ap.add_argument("--diamond-version", default="unknown")
    ap.add_argument("--source-db-name", required=True)
    ap.add_argument("--source-db-release", required=True)
    ap.add_argument("--source-db-download-date", required=True)
    ap.add_argument("--prefix", required=True)
    args = ap.parse_args()

    ablation_protein = load_json(args.ablation_json)
    ablation_nucleotide = load_json(args.ablation_json_nt) if args.ablation_json_nt else None
    kraken2_identity = load_json(args.kraken2_identity_json)
    diamond_identity = load_json(args.diamond_identity_json)

    provenance = {
        "source_database": {
            "name": args.source_db_name,
            "release": args.source_db_release,
            "download_date": args.source_db_download_date,
        },
        "ablation_protein": ablation_protein,
        "ablation_nucleotide": ablation_nucleotide,
        "measured_identity": {
            "kraken2_reference": kraken2_identity or None,
            "diamond_reference": diamond_identity or None,
            "note": (
                "Both are reported from the same DIAMOND blastp-against-ablated-protein-DB "
                "measurement (docs/DESIGN.md section 3); duplicated here only if the two "
                "classifiers were built from different ablated protein sets in a given run."
            ),
        },
        "builds": [],
    }

    if args.kraken2_db_files:
        provenance["builds"].append(
            {
                "classifier": "kraken2",
                "tool_version": args.kraken2_version,
                "build_command": (
                    "kraken2-build --add-to-library <ablated + host + decoy FASTA> --db <db> ; "
                    "kraken2-build --build --db <db>  "
                    "(template reconstructed from modules/nf-core/kraken2/{add,build}/main.nf "
                    "with this run's task.ext.args; see note above)"
                ),
                "index_checksum_sha256": sha256_of_paths(args.kraken2_db_files),
                "index_files": sorted(str(p.name) for p in args.kraken2_db_files),
            }
        )

    if args.diamond_db_files:
        provenance["builds"].append(
            {
                "classifier": "diamond",
                "tool_version": args.diamond_version,
                "build_command": (
                    "diamond makedb --in <ablated + host + decoy protein FASTA> -d <prefix> "
                    "(template reconstructed from modules/nf-core/diamond/makedb/main.nf "
                    "with this run's task.ext.args; see note above)"
                ),
                "index_checksum_sha256": sha256_of_paths(args.diamond_db_files),
                "index_files": sorted(str(p.name) for p in args.diamond_db_files),
            }
        )

    if not provenance["builds"]:
        sys.exit("error: no database files were passed -- refusing to write a provenance record for nothing built")

    out_path = Path(f"{args.prefix}.provenance.json")
    out_path.write_text(json.dumps(provenance, indent=2, sort_keys=True) + "\n")
    print(f"build_provenance: wrote {out_path} covering {[b['classifier'] for b in provenance['builds']]}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
