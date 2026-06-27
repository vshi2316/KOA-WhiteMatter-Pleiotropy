#!/usr/bin/env python3

from pathlib import Path
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

PROJECT_DIR = Path(os.environ.get("KOA_WM_PROJECT_DIR", "."))
DATA_DIR = Path(os.environ.get("KOA_WM_DATA_DIR", PROJECT_DIR / "data"))
RESULTS_DIR = Path(os.environ.get("KOA_WM_RESULTS_DIR", PROJECT_DIR / "results"))
TABLES_DIR = RESULTS_DIR / "tables"
FIGURES_DIR = RESULTS_DIR / "figures"
TABLES_DIR.mkdir(parents=True, exist_ok=True)
FIGURES_DIR.mkdir(parents=True, exist_ok=True)

MHC = ("6", 25119106, 33854733)
REGION_8P23 = ("8", 7200000, 12500000)


def read_table(path):
    return pd.read_csv(path, sep=None, engine="python")


def normalize(df):
    cols = {c.upper(): c for c in df.columns}
    out = df.rename(columns={
        cols.get("RSID", cols.get("SNP", "SNP")): "SNP",
        cols.get("CHR", "CHR"): "CHR",
        cols.get("BP", cols.get("POS", "BP")): "BP",
        cols.get("A1", "A1"): "A1",
        cols.get("A2", "A2"): "A2",
        cols.get("Z", "Z"): "Z",
        cols.get("P", cols.get("PVAL", "P")): "P",
    })
    out["CHR"] = out["CHR"].astype(str).str.replace("chr", "", regex=False)
    out["BP"] = pd.to_numeric(out["BP"], errors="coerce")
    out["P"] = pd.to_numeric(out["P"], errors="coerce").clip(lower=1e-300, upper=1)
    return out.dropna(subset=["SNP", "CHR", "BP", "P"])


def exclude_regions(df):
    in_mhc = (df["CHR"] == MHC[0]) & df["BP"].between(MHC[1], MHC[2])
    in_8p23 = (df["CHR"] == REGION_8P23[0]) & df["BP"].between(REGION_8P23[1], REGION_8P23[2])
    return df.loc[~(in_mhc | in_8p23)].copy()


def empirical_condfdr(p_primary, p_secondary):
    order = np.argsort(p_secondary)
    p1 = p_primary[order]
    rank_primary = pd.Series(p1).rank(method="max").to_numpy()
    rank_secondary = np.arange(1, len(p1) + 1)
    cdf_primary = rank_primary / len(p1)
    cdf_cond = rank_primary / rank_secondary
    out_sorted = np.minimum(1.0, cdf_primary / np.maximum(cdf_cond, 1e-12))
    out = np.empty_like(out_sorted)
    out[order] = out_sorted
    return out


def condqq(df, primary, secondary, label, out_png):
    plt.figure(figsize=(5, 5))
    for cut, lab in [(1, "all"), (0.1, "P2<=0.1"), (0.01, "P2<=0.01"), (0.001, "P2<=0.001")]:
        sub = df if cut == 1 else df[df[secondary] <= cut]
        if len(sub) < 50:
            continue
        obs = -np.log10(np.sort(sub[primary].clip(lower=1e-300).values))
        exp = -np.log10(np.arange(1, len(obs) + 1) / (len(obs) + 1))
        plt.plot(exp, obs, lw=1, label=lab)
    lim = max(plt.xlim()[1], plt.ylim()[1])
    plt.plot([0, lim], [0, lim], color="grey", ls="--", lw=0.8)
    plt.xlabel("Expected -log10(P)")
    plt.ylabel("Observed -log10(P)")
    plt.title(label)
    plt.legend(frameon=False, fontsize=8)
    plt.tight_layout()
    plt.savefig(out_png, dpi=300)
    plt.close()


manifest = read_table(os.environ.get("KOA_WM_CONDFDR_PAIR_MANIFEST", DATA_DIR / "condfdr_pair_manifest.tsv"))
cond_tables = []
conj_tables = []

for _, row in manifest.iterrows():
    pair_id = str(row["pair_id"])
    wm_label = str(row["wm_label"])
    koa = normalize(read_table(row["koa_file"]))
    wm = normalize(read_table(row["wm_file"]))
    m = koa.merge(wm, on="SNP", suffixes=("_KOA", "_WM"))
    m["CHR"] = m.get("CHR_KOA", m.get("CHR_WM")).astype(str)
    m["BP"] = pd.to_numeric(m.get("BP_KOA", m.get("BP_WM")), errors="coerce")
    m = exclude_regions(m)
    m["condFDR_KOA_given_WM"] = empirical_condfdr(m["P_KOA"].to_numpy(), m["P_WM"].to_numpy())
    m["condFDR_WM_given_KOA"] = empirical_condfdr(m["P_WM"].to_numpy(), m["P_KOA"].to_numpy())
    m["conjFDR"] = np.maximum(m["condFDR_KOA_given_WM"], m["condFDR_WM_given_KOA"])
    m["pair_id"] = pair_id
    m["wm_label"] = wm_label
    cond_tables.append(m[(m["condFDR_KOA_given_WM"] < 0.01) | (m["condFDR_WM_given_KOA"] < 0.01)].copy())
    conj_tables.append(m[m["conjFDR"] < 0.05].copy())
    condqq(m, "P_KOA", "P_WM", f"KOA conditioned on {wm_label}", FIGURES_DIR / f"Supplementary_Figure_condQQ_{pair_id}_KOA_given_WM.png")
    condqq(m, "P_WM", "P_KOA", f"{wm_label} conditioned on KOA", FIGURES_DIR / f"Supplementary_Figure_condQQ_{pair_id}_WM_given_KOA.png")

cond = pd.concat(cond_tables, ignore_index=True) if cond_tables else pd.DataFrame()
conj = pd.concat(conj_tables, ignore_index=True) if conj_tables else pd.DataFrame()
cond.to_csv(TABLES_DIR / "Supplementary_Table_6_7_condFDR_outputs.tsv", sep="\t", index=False)
conj.to_csv(TABLES_DIR / "conjFDR_SNP_level_enrichment_records.tsv", sep="\t", index=False)
conj.groupby(["pair_id", "wm_label"], as_index=False).size().rename(columns={"size": "n_conjFDR_records"}).to_csv(
    TABLES_DIR / "Figure_2A_F_condFDR_conjFDR_count_source.tsv", sep="\t", index=False
)
print("Completed KOA-WM condFDR/conjFDR enrichment analysis.")
