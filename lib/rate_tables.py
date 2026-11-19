"""
rate_tables.py - reader for the pairwise and annotated substitution-rate tables.

Every script reads its input through `read_rate_table()` and works in the
manuscript's notation: dN, dS and dN_dS. Tables written with the legacy Ka, Ks
and Ka_Ks column names are mapped on load, so a script never needs to know which
convention the file on disk uses, and no table has to be rewritten to be read.

Matching is case- and separator-insensitive, so Ka_Ks, KaKs, ka/ks and KA_KS all
resolve to dN_dS. A table carrying two spellings of the same quantity raises
rather than silently choosing between columns that may disagree.

`to_publication()` applies the column names used in the supplementary tables --
gene_ascaris rather than gene1, distance_to_nearest_CBR_bp rather than
min_cbr_dist1. It is applied only when a table is exported, so working files keep
whatever names the pipeline gave them.
"""
from __future__ import annotations

import pandas as pd

__all__ = ["ALIASES", "normalize_columns", "read_rate_table",
           "PUBLICATION_COLUMNS", "to_publication"]


#: Legacy column name -> canonical name.
ALIASES = {
    "ka": "dN",
    "ks": "dS",
    "kaks": "dN_dS",
    "ka_ks": "dN_dS",
    "kaks_ratio": "dN_dS",
    "omega": "dN_dS",
    "dn": "dN",
    "ds": "dS",
    "dnds": "dN_dS",
    "dn_ds": "dN_dS",
}


def _key(name: str) -> str:
    return name.strip().lower().replace("/", "_").replace("-", "_").replace(" ", "_")


def normalize_columns(df: pd.DataFrame, source: str = "<dataframe>") -> pd.DataFrame:
    """Rename rate columns to dN / dS / dN_dS, leaving everything else alone.

    Raises if a table carries both spellings of the same quantity, rather than
    choosing between two columns that may not agree.
    """
    rename, seen = {}, {}
    for col in df.columns:
        canon = ALIASES.get(_key(col))
        if canon is None:
            continue
        if canon in seen:
            raise ValueError(
                f"{source}: columns {seen[canon]!r} and {col!r} both map to "
                f"{canon!r}. Resolve the duplicate before loading.")
        seen[canon] = col
        if col != canon:
            rename[col] = canon
    return df.rename(columns=rename) if rename else df


def read_rate_table(path: str, sep: str = "\t", require=("dS",), **kw) -> pd.DataFrame:
    """Read a pairwise or annotated rate table with canonical column names.

    require : columns that must be present after normalization. Reporting the
              file path and its actual columns here is more informative than a
              KeyError raised further down the call stack.
    """
    df = pd.read_csv(path, sep=sep, **kw)
    df = normalize_columns(df, source=path)
    missing = [c for c in require if c not in df.columns]
    if missing:
        raise KeyError(
            f"{path}: missing required column(s) {missing}. "
            f"Columns present: {list(df.columns)}")
    return df


#: Internal name -> the name used in a supplementary table.
PUBLICATION_COLUMNS = {
    "gene1": "gene_ascaris",
    "gene2": "gene_parascaris",
    "chrom1": "chromosome_ascaris",
    "chrom2": "chromosome_parascaris",
    "region1": "region_ascaris",
    "region2": "region_parascaris",
    "consensus_region": "region_consensus",
    "tox_region": "region_toxocara",
    "partner_region": "region_partner_annotation",
    "min_cbr_dist1": "distance_to_nearest_CBR_bp",
    "cbr_proximal": "CBR_proximal_100kb",
    "eliminated_type": "eliminated_subtype",
    "dN": "dN",
    "dS": "dS",
    "dN_dS": "dN_dS",
    "gc3": "GC3",
    "gc": "GC",
    "peak_stage": "peak_expression_stage",
}


def to_publication(df: pd.DataFrame, order=None, extra=None) -> pd.DataFrame:
    """Return a copy with publication column names, optionally reordered.

    order : sequence of internal names fixing the column order; anything not
            listed is appended in its existing order.
    extra : additional internal -> published name mappings for one table.
    """
    mapping = dict(PUBLICATION_COLUMNS)
    if extra:
        mapping.update(extra)
    out = normalize_columns(df).copy()
    if order:
        first = [c for c in order if c in out.columns]
        out = out[first + [c for c in out.columns if c not in first]]
    return out.rename(columns={k: v for k, v in mapping.items()
                               if k in out.columns})
