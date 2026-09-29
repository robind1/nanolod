#!/usr/bin/env python3


from __future__ import annotations

import argparse
import json
import statistics
import sys
from pathlib import Path

DEFAULT_COLUMNS = [
    "qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
    "qstart", "qend", "sstart", "send", "evalue", "bitscore",
]


def load_all_query_ids(fasta_path: Path) -> list[str]:
    ids = []
    with open(fasta_path, encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if line.startswith(">"):
                ids.append(line[1:].split()[0].strip())
    return ids


def load_best_hits(tsv_path: Path, columns: list[str]) -> dict[str, dict]:
    qi, si, pi, bi = columns.index("qseqid"), columns.index("sseqid"), columns.index("pident"), columns.index("bitscore")
    best: dict[str, dict] = {}
    if not tsv_path.exists() or tsv_path.stat().st_size == 0:
        return best
    with open(tsv_path, encoding="utf-8", errors="replace") as fh:
        for line in fh:
            fields = line.rstrip("\n").split("\t")
            if len(fields) <= max(qi, si, pi, bi):
                continue
            qseqid, sseqid, pident, bitscore = fields[qi], fields[si], float(fields[pi]), float(fields[bi])
            current = best.get(qseqid)
            if current is None or bitscore > current["bitscore"]:
                best[qseqid] = {"sseqid": sseqid, "pident": pident, "bitscore": bitscore}
    return best


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--blastp-tsv", required=True, type=Path, help="DIAMOND blastp outfmt 6 output")
    ap.add_argument("--query-fasta", required=True, type=Path, help="the target protein FASTA that was searched")
    ap.add_argument("--columns", default=",".join(DEFAULT_COLUMNS), help="comma-separated blastp output column names")
    ap.add_argument("--fail-above", type=float, default=None, help="exit non-zero if max identity exceeds this (docs/DESIGN.md section 3: ~80 for L4)")
    ap.add_argument("--label", default="")
    ap.add_argument("--prefix", required=True)
    args = ap.parse_args()

    columns = [c.strip() for c in args.columns.split(",")]
    for required in ("qseqid", "sseqid", "pident", "bitscore"):
        if required not in columns:
            sys.exit(f"error: --columns must include '{required}'; got {columns}")

    query_ids = load_all_query_ids(args.query_fasta)
    if not query_ids:
        sys.exit(f"error: {args.query_fasta} contains no proteins to summarise")

    best_hits = load_best_hits(args.blastp_tsv, columns)

    per_protein = []
    for qid in query_ids:
        hit = best_hits.get(qid)
        if hit is None:
            per_protein.append({"query": qid, "identity": 0.0, "best_hit": None, "no_hit": True})
        else:
            per_protein.append({"query": qid, "identity": hit["pident"], "best_hit": hit["sseqid"], "no_hit": False})

    identities = [p["identity"] for p in per_protein]
    n_no_hit = sum(1 for p in per_protein if p["no_hit"])

    summary = {
        "label": args.label,
        "n_proteins": len(per_protein),
        "n_no_hit": n_no_hit,
        "max_identity": max(identities),
        "median_identity": statistics.median(identities),
    }

    Path(f"{args.prefix}.identity_summary.json").write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")

    with open(f"{args.prefix}.identity_per_protein.tsv", "w", encoding="utf-8") as fh:
        fh.write("query\tidentity\tbest_hit\tno_hit\n")
        for p in per_protein:
            fh.write(f"{p['query']}\t{p['identity']}\t{p['best_hit'] or ''}\t{p['no_hit']}\n")

    print(
        f"summarize_identity ({args.label}): max={summary['max_identity']:.2f} "
        f"median={summary['median_identity']:.2f} over {summary['n_proteins']} proteins "
        f"({n_no_hit} with no hit)",
        file=sys.stderr,
    )

    if args.fail_above is not None and summary["max_identity"] > args.fail_above:
        sys.exit(
   
        )

    return 0


if __name__ == "__main__":
    sys.exit(main())
