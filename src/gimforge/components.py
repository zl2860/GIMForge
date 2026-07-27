"""Turn an ordered conditional association matrix into GIM components."""

from __future__ import annotations

from collections import defaultdict
from statistics import median
from typing import Iterable, Mapping

from .frequency import classify_maf, maf_from_a1_frequency, normalise_allele_frequency, normalise_maf
from .io import as_float, as_int


class _UnionFind:
    def __init__(self) -> None:
        self.parent: dict[str, str] = {}

    def add(self, item: str) -> None:
        self.parent.setdefault(item, item)

    def find(self, item: str) -> str:
        parent = self.parent[item]
        if parent != item:
            self.parent[item] = self.find(parent)
        return self.parent[item]

    def union(self, left: str, right: str) -> None:
        self.add(left)
        self.add(right)
        left_root, right_root = self.find(left), self.find(right)
        if left_root != right_root:
            self.parent[right_root] = left_root


def _column(rows: list[Mapping[str, object]], choices: tuple[str, ...], label: str) -> str:
    if not rows:
        raise ValueError("matrix_out is empty.")
    available = set(rows[0])
    for candidate in choices:
        if candidate in available:
            return candidate
    raise ValueError(f"matrix_out requires {label}; accepted columns are {', '.join(choices)}.")


def components_from_matrix(
    matrix_out: Iterable[Mapping[str, object]], *, conditional_p: float
) -> dict[str, list[dict[str, object]]]:
    """Define GIMs as significant connected components of the bipartite matrix.

    ``matrix_out`` must be the ordered `V_R x M_R` matrix: every row is a
    test of one selected marker against one region member, conditioned on
    markers with lower ``marker_order``. This function intentionally does not
    re-order or re-fit that matrix; it is the final definition step.
    """

    rows = [dict(row) for row in matrix_out]
    p_column = _column(rows, ("p", "P", "p_value"), "a P-value column")
    snp_column = _column(rows, ("snp_id", "ID", "markername"), "a SNP column")
    metabolite_column = _column(rows, ("metabolite", "phenotype", "trait"), "a metabolite column")
    region_column = _column(rows, ("region_id", "region"), "a region column")
    order_column = next((name for name in ("marker_order", "order") if name in rows[0]), None)
    beta_column = next((name for name in ("beta", "BETA", "effect") if name in rows[0]), None)
    se_column = next((name for name in ("se", "SE", "standard_error") if name in rows[0]), None)
    n_column = next((name for name in ("n", "N", "OBS_CT") if name in rows[0]), None)
    a1_freq_column = next(
        (name for name in ("a1_freq", "A1_FREQ", "A1FREQ") if name in rows[0]),
        None,
    )
    maf_column = next((name for name in ("maf", "MAF") if name in rows[0]), None)

    by_region: dict[str, list[dict[str, object]]] = defaultdict(list)
    variant_mafs: dict[tuple[str, str], list[float]] = defaultdict(list)
    normalised_rows: list[dict[str, object]] = []
    for source_row in rows:
        row = dict(source_row)
        p_value = as_float(row.get(p_column))
        row["p"] = p_value if p_value is not None else ""
        row["snp_id"] = str(row[snp_column])
        row["metabolite"] = str(row[metabolite_column])
        row["region_id"] = str(row[region_column])
        row["marker_order"] = as_int(row.get(order_column)) if order_column else ""
        if beta_column is not None:
            row["beta"] = row.get(beta_column, "")
        if se_column is not None:
            row["se"] = row.get(se_column, "")
        if n_column is not None:
            row["n"] = row.get(n_column, "")
        a1_freq = (
            normalise_allele_frequency(row.get(a1_freq_column))
            if a1_freq_column is not None
            else None
        )
        maf = normalise_maf(row.get(maf_column)) if maf_column is not None else None
        if maf is None:
            maf = maf_from_a1_frequency(a1_freq)
        row["a1_freq"] = a1_freq if a1_freq is not None else ""
        row["maf"] = maf if maf is not None else ""
        row["maf_class"] = classify_maf(maf)
        if maf is not None:
            variant_mafs[(row["region_id"], row["snp_id"])].append(maf)
        normalised_rows.append(row)
        if p_value is None or p_value < 0 or p_value > conditional_p:
            continue
        by_region[row["region_id"]].append(row)

    edges: list[dict[str, object]] = []
    members: list[dict[str, object]] = []
    summary: list[dict[str, object]] = []
    for region_id in sorted(by_region):
        region_edges = by_region[region_id]
        graph = _UnionFind()
        marker_order: dict[str, int] = {}
        for row in region_edges:
            snp_node, metabolite_node = f"S:{row['snp_id']}", f"M:{row['metabolite']}"
            graph.union(snp_node, metabolite_node)
            if row["marker_order"] is not None:
                current = marker_order.get(snp_node)
                marker_order[snp_node] = min(current, row["marker_order"]) if current is not None else row["marker_order"]

        component_nodes: dict[str, list[str]] = defaultdict(list)
        for node in graph.parent:
            component_nodes[graph.find(node)].append(node)

        def component_key(nodes: list[str]) -> tuple[int, str]:
            orders = [marker_order[node] for node in nodes if node in marker_order]
            return (min(orders) if orders else 2**31 - 1, min(nodes))

        ordered_components = sorted(component_nodes.values(), key=component_key)
        node_to_gim: dict[str, str] = {}
        for index, nodes in enumerate(ordered_components, start=1):
            gim_id = f"{region_id}_GIM_{index:03d}"
            snps = sorted(node[2:] for node in nodes if node.startswith("S:"))
            metabolites = sorted(node[2:] for node in nodes if node.startswith("M:"))
            snp_maf = {
                snp: (
                    median(variant_mafs[(region_id, snp)])
                    if variant_mafs.get((region_id, snp))
                    else None
                )
                for snp in snps
            }
            snp_maf_class = {
                snp: classify_maf(snp_maf[snp])
                for snp in snps
            }
            for node in nodes:
                node_to_gim[node] = gim_id
                snp_id = node[2:] if node.startswith("S:") else None
                members.append(
                    {
                        "gim_id": gim_id,
                        "region_id": region_id,
                        "node_type": "SNP" if node.startswith("S:") else "metabolite",
                        "node_id": node[2:],
                        "marker_order": marker_order.get(node, ""),
                        "maf": (
                            snp_maf[snp_id]
                            if snp_id is not None and snp_maf[snp_id] is not None
                            else ""
                        ),
                        "maf_class": (
                            snp_maf_class[snp_id] if snp_id is not None else ""
                        ),
                    }
                )
            class_counts = {
                name: sum(value == name for value in snp_maf_class.values())
                for name in ("rare", "low_frequency", "common", "unknown")
            }
            summary.append(
                {
                    "gim_id": gim_id,
                    "region_id": region_id,
                    "n_snps": len(snps),
                    "n_rare_snps": class_counts["rare"],
                    "n_low_frequency_snps": class_counts["low_frequency"],
                    "n_common_snps": class_counts["common"],
                    "n_unknown_maf_snps": class_counts["unknown"],
                    "n_metabolites": len(metabolites),
                    "snps": ";".join(snps),
                    "snp_mafs": ";".join(
                        f"{snp}={snp_maf[snp]:.8g}"
                        if snp_maf[snp] is not None
                        else f"{snp}=NA"
                        for snp in snps
                    ),
                    "snp_maf_classes": ";".join(
                        f"{snp}={snp_maf_class[snp]}" for snp in snps
                    ),
                    "metabolites": ";".join(metabolites),
                }
            )
        for row in region_edges:
            edge = dict(row)
            edge["gim_id"] = node_to_gim[f"S:{row['snp_id']}"]
            edges.append(edge)

    return {
        "matrix_out": normalised_rows,
        "edges": edges,
        "members": members,
        "gim_summary": summary,
    }
