#!/usr/bin/env python3
"""
plot_variant_summary.py
------------------------
Companion script for food_allergy_pipeline_v3.sh.

NOTE FOR REVIEWERS: this file was referenced by the original v2.0 pipeline
(`python3 "$(dirname "$0")/plot_variant_summary.py" ...`) but was never
included in the submitted bundle, so the summarize_sample() step would have
failed (silently, via `|| warn ...`) on every single run. This is a minimal,
functional implementation provided so the pipeline is actually complete and
runnable end-to-end; extend it as needed for your dissertation figures.

Produces, per sample, in --outdir:
  <sample>_variant_counts_by_gene.png   - STAT6 vs FCER1A variant counts
  <sample>_snp_indel_breakdown.png      - SNP vs INDEL counts (from bcftools stats)
  <sample>_quality_distribution.png     - QUAL score histogram, both genes combined
"""

import argparse
import re
import sys
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


def parse_args():
    p = argparse.ArgumentParser(description="Plot per-sample variant summary figures.")
    p.add_argument("--sample", required=True)
    p.add_argument("--stats", required=True, help="bcftools stats output file")
    p.add_argument("--stat6-vcf", required=True)
    p.add_argument("--fcer1a-vcf", required=True)
    p.add_argument("--outdir", required=True)
    return p.parse_args()


def count_records(vcf_path):
    n = 0
    try:
        with open(vcf_path) as fh:
            for line in fh:
                if not line.startswith("#"):
                    n += 1
    except FileNotFoundError:
        pass
    return n


def qual_values(vcf_path):
    quals = []
    try:
        with open(vcf_path) as fh:
            for line in fh:
                if line.startswith("#"):
                    continue
                fields = line.rstrip("\n").split("\t")
                if len(fields) >= 6:
                    try:
                        quals.append(float(fields[5]))
                    except ValueError:
                        pass
    except FileNotFoundError:
        pass
    return quals


def parse_snp_indel(stats_path):
    snp = indel = 0
    try:
        with open(stats_path) as fh:
            for line in fh:
                if line.startswith("SN"):
                    if "number of SNPs:" in line:
                        snp = int(line.strip().split("\t")[-1])
                    elif "number of indels:" in line:
                        indel = int(line.strip().split("\t")[-1])
    except FileNotFoundError:
        pass
    return snp, indel


def main():
    args = parse_args()

    stat6_n = count_records(args.stat6_vcf)
    fcer1a_n = count_records(args.fcer1a_vcf)

    # --- Figure 1: variant counts by gene ---
    fig, ax = plt.subplots(figsize=(5, 4))
    genes = ["STAT6", "FCER1A"]
    counts = [stat6_n, fcer1a_n]
    ax.bar(genes, counts, color=["#4C72B0", "#DD8452"])
    ax.set_ylabel("PASS variant count")
    ax.set_title(f"{args.sample}: variants per target gene")
    for i, c in enumerate(counts):
        ax.text(i, c, str(c), ha="center", va="bottom")
    fig.tight_layout()
    fig.savefig(f"{args.outdir}/{args.sample}_variant_counts_by_gene.png", dpi=150)
    plt.close(fig)

    # --- Figure 2: SNP vs INDEL breakdown (from bcftools stats) ---
    snp, indel = parse_snp_indel(args.stats)
    fig, ax = plt.subplots(figsize=(5, 4))
    ax.bar(["SNPs", "Indels"], [snp, indel], color=["#55A868", "#C44E52"])
    ax.set_ylabel("Count")
    ax.set_title(f"{args.sample}: SNP vs Indel (target regions)")
    fig.tight_layout()
    fig.savefig(f"{args.outdir}/{args.sample}_snp_indel_breakdown.png", dpi=150)
    plt.close(fig)

    # --- Figure 3: QUAL distribution, both genes combined ---
    quals = qual_values(args.stat6_vcf) + qual_values(args.fcer1a_vcf)
    fig, ax = plt.subplots(figsize=(5, 4))
    if quals:
        ax.hist(quals, bins=20, color="#8172B2", edgecolor="black")
    ax.set_xlabel("QUAL")
    ax.set_ylabel("Number of variants")
    ax.set_title(f"{args.sample}: QUAL score distribution")
    fig.tight_layout()
    fig.savefig(f"{args.outdir}/{args.sample}_quality_distribution.png", dpi=150)
    plt.close(fig)

    print(f"[{args.sample}] Plots written to {args.outdir}")


if __name__ == "__main__":
    sys.exit(main())
