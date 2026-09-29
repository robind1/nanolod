#!/usr/bin/env python3

from __future__ import annotations

import argparse
import csv
import json
import sys
from collections import defaultdict
from pathlib import Path

RANK_LADDER = {
    "none": None,
    "species": "species",
    "genus": "genus",
    "subfamily": "subfamily",
    "family": "family",
}


def parse_dmp(path: Path) -> list[list[str]]:
    rows = []
    with open(path, encoding="utf-8", errors="replace") as fh:
        for line in fh:
            line = line.rstrip("\n").rstrip("\t|")
            fields = [f.strip() for f in line.split("\t|\t")]
            if fields and fields != [""]:
                rows.append(fields)
    return rows


def load_nodes(path: Path) -> tuple[dict[str, str], dict[str, str]]:
    """Return (taxid -> parent_taxid, taxid -> rank)."""
    parent: dict[str, str] = {}
    rank: dict[str, str] = {}
    for fields in parse_dmp(path):
        if len(fields) < 3:
            continue
        taxid, parent_taxid, node_rank = fields[0], fields[1], fields[2]
        parent[taxid] = parent_taxid
        rank[taxid] = node_rank
    return parent, rank


def load_names(path: Path | None) -> dict[str, str]:
    names: dict[str, str] = {}
    if path is None:
        return names
    for fields in parse_dmp(path):
        if len(fields) < 4:
            continue
        taxid, name, _unique, name_class = fields[0], fields[1], fields[2], fields[3]
        if name_class == "scientific name":
            names[taxid] = name
    return names


def load_accession2taxid(path: Path) -> dict[str, str]:

    with open(path, newline="", encoding="utf-8") as fh:
        sample = fh.read(4096)
        fh.seek(0)
        dialect = csv.Sniffer().sniff(sample, delimiters="\t,") if sample.strip() else csv.excel_tab
        reader = csv.DictReader(fh, dialect=dialect)
        if reader.fieldnames is None:
            sys.exit(f"error: {path} has no header row; expected an accession and a taxid column")
        lower_to_actual = {f.lower(): f for f in reader.fieldnames}
        acc_col = lower_to_actual.get("accession.version") or lower_to_actual.get("accession")
        taxid_col = lower_to_actual.get("taxid")
        if acc_col is None or taxid_col is None:
            sys.exit(
                f"error: {path} must have an 'accession' (or 'accession.version') column and a "
                f"'taxid' column; found {reader.fieldnames}"
            )
        mapping: dict[str, str] = {}
        for row in reader:
            acc = row[acc_col].strip()
            taxid = row[taxid_col].strip()
            if not acc or not taxid:
                continue
            if acc in mapping and mapping[acc] != taxid:
                sys.exit(
                    f"error: accession {acc!r} maps to two different taxids "
                    f"({mapping[acc]!r} and {taxid!r}) in {path}. This is a data problem "
                    "in the reference manifest, not something to silently resolve here."
                )
            mapping[acc] = taxid
    return mapping


def find_ancestor_at_rank(taxid: str, target_rank: str, parent: dict[str, str], rank: dict[str, str]) -> str | None:
    seen: set[str] = set()
    node = taxid
    while node and node not in seen:
        seen.add(node)
        if rank.get(node) == target_rank:
            return node
        if node == "1" or node not in parent:
            break
        next_node = parent[node]
        if next_node == node: 
            break
        node = next_node
    return None


def collect_subtree(root_taxid: str, parent: dict[str, str]) -> set[str]:
    children: dict[str, list[str]] = defaultdict(list)
    for child, par in parent.items():
        if child != par:
            children[par].append(child)

    subtree = set()
    stack = [root_taxid]
    while stack:
        node = stack.pop()
        if node in subtree:
            continue
        subtree.add(node)
        stack.extend(children.get(node, []))
    return subtree


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--nodes", required=True, type=Path)
    ap.add_argument("--names", type=Path, default=None)
    ap.add_argument("--accession2taxid", required=True, type=Path)
    ap.add_argument("--target-taxid", default=None)
    ap.add_argument("--holdout-rank", required=True, choices=sorted(RANK_LADDER))
    ap.add_argument("--holdout-level", default=None, help="cosmetic label (e.g. L4) carried into the JSON")
    ap.add_argument("--prefix", required=True)
    args = ap.parse_args()

    target_rank = RANK_LADDER[args.holdout_rank]

    if target_rank is not None and not args.target_taxid:
        sys.exit(f"error: --holdout-rank {args.holdout_rank} requires --target-taxid")

    accession2taxid = load_accession2taxid(args.accession2taxid)
    if not accession2taxid:
        sys.exit(f"error: {args.accession2taxid} produced zero accession/taxid pairs")

    all_accessions = set(accession2taxid)
    summary = {
        "target_taxid": args.target_taxid,
        "holdout_rank": args.holdout_rank,
        "holdout_level": args.holdout_level,
        "n_accessions_total": len(all_accessions),
    }

    if target_rank is None:
        include_accessions = all_accessions
        exclude_accessions: set[str] = set()
        summary.update(
            {
                "ancestor_taxid": None,
                "ancestor_rank": None,
                "ancestor_name": None,
                "n_excluded": 0,
                "n_included": len(include_accessions),
            }
        )
    else:
        parent, rank = load_nodes(args.nodes)
        names = load_names(args.names)

        if args.target_taxid not in rank:
            sys.exit(
                f"error: target taxid {args.target_taxid} is not present in {args.nodes}. "
                "Check the taxid against the taxonomy dump version actually in use."
            )

        ancestor = find_ancestor_at_rank(args.target_taxid, target_rank, parent, rank)
        if ancestor is None:
            sys.exit(

            )

        exclude_taxids = collect_subtree(ancestor, parent)
        exclude_accessions = {acc for acc, taxid in accession2taxid.items() if taxid in exclude_taxids}
        include_accessions = all_accessions - exclude_accessions

        if not exclude_accessions:
            sys.exit(

            )

        summary.update(
            {
                "ancestor_taxid": ancestor,
                "ancestor_rank": rank.get(ancestor),
                "ancestor_name": names.get(ancestor),
                "n_excluded": len(exclude_accessions),
                "n_included": len(include_accessions),
            }
        )

    prefix = args.prefix
    Path(f"{prefix}.include.accessions.txt").write_text(
        "\n".join(sorted(include_accessions)) + ("\n" if include_accessions else "")
    )
    Path(f"{prefix}.exclude.accessions.txt").write_text(
        "\n".join(sorted(exclude_accessions)) + ("\n" if exclude_accessions else "")
    )
    Path(f"{prefix}.ablation.json").write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")

    print(
        f"{args.holdout_rank}: {summary['n_included']} included, {summary['n_excluded']} excluded "
        f"(of {summary['n_accessions_total']} total)",
        file=sys.stderr,
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
