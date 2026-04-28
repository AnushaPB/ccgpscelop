#!/usr/bin/env python3

from cyvcf2 import VCF
import pandas as pd
from collections import Counter
import os

"""
Compute nucleotide diversity (π) in three adjacent windows (left, center, right)
around each SNP and compare the center window to its flanking regions. The focal
SNP itself is excluded from the center window calculation to avoid biasing the result due 
to the SNP being variable by definition.
"""

# -----------------------------
# input files
# -----------------------------
VCF_FILE = "../../data/ccgp_data/58-Sceloporus_complete_coords_annotated_final_samples.vcf.gz"
SNP_FILE = "outputs/bio1ndvi_gea_genes_nonsyn_ids.txt"
OUT_FILE = "outputs/bio1ndvi_gea_genes_nonsyn_flanking_pi.tsv"

WINDOW_SIZE = 50000

# -----------------------------
# helper functions
# -----------------------------
def parse_snp_id(s):
    """
    Parse SNP IDs like:
    chr1_65730_C_G
    into chrom, pos, ref, alt
    """
    parts = s.strip().split("_")
    if len(parts) < 4:
        raise ValueError(f"Could not parse SNP ID: {s}")
    chrom = parts[0]
    pos = int(parts[1])
    ref = parts[2]
    alt = parts[3]
    return chrom, pos, ref, alt


def site_pi_from_genotypes(genotypes):
    """
    Calculate per-site nucleotide diversity (average pairwise differences)
    from diploid genotypes.

    cyvcf2 genotypes look like:
    [[0, 0, False], [0, 1, False], [1, 1, False], ...]

    Missing alleles are coded as -1.
    Works for multiallelic sites too.
    Returns None if fewer than 2 chromosomes are called.
    """
    allele_counts = Counter()

    for gt in genotypes:
        a1, a2 = gt[0], gt[1]
        if a1 != -1:
            allele_counts[a1] += 1
        if a2 != -1:
            allele_counts[a2] += 1

    n_chr = sum(allele_counts.values())
    if n_chr < 2:
        return None

    same_pairs = sum(c * (c - 1) / 2 for c in allele_counts.values())
    total_pairs = n_chr * (n_chr - 1) / 2
    pi_site = 1 - (same_pairs / total_pairs)

    return pi_site


def window_pi(vcf, chrom, start, end, exclude_pos=None):
    """
    Sum site pi across all variant sites in the window and divide by window length
    to get per-base pi.

    Parameters
    ----------
    vcf : cyvcf2.VCF
    chrom : str
    start : int
    end : int
    exclude_pos : int or None
        If provided, exclude that exact genomic position from the calculation.
        Useful for removing the focal SNP from the center window.

    Notes
    -----
    This assumes sites absent from the VCF contribute 0, which is standard if your
    VCF represents callable variation and you are comfortable treating other sites
    in the window as invariant/callable.
    """
    if end < start:
        return {
            "pi": None,
            "pi_sum": 0.0,
            "window_length": 0,
            "n_variant_sites": 0,
            "n_sites_used": 0,
        }

    region = f"{chrom}:{start}-{end}"

    pi_sum = 0.0
    n_variant_sites = 0
    n_sites_used = 0

    try:
        for var in vcf(region):
            if exclude_pos is not None and var.POS == exclude_pos:
                continue

            n_variant_sites += 1
            pi_site = site_pi_from_genotypes(var.genotypes)
            if pi_site is not None:
                pi_sum += pi_site
                n_sites_used += 1
    except Exception as e:
        print(f"Warning: failed on region {region}: {e}")

    window_len = end - start + 1
    pi_per_base = pi_sum / window_len if window_len > 0 else None

    return {
        "pi": pi_per_base,
        "pi_sum": pi_sum,
        "window_length": window_len,
        "n_variant_sites": n_variant_sites,
        "n_sites_used": n_sites_used,
    }


def mean_ignore_none(values):
    vals = [x for x in values if x is not None]
    if not vals:
        return None
    return sum(vals) / len(vals)


# -----------------------------
# main
# -----------------------------
def main():
    if not os.path.exists(VCF_FILE):
        raise FileNotFoundError(f"VCF not found: {VCF_FILE}")
    if not os.path.exists(SNP_FILE):
        raise FileNotFoundError(f"SNP list not found: {SNP_FILE}")

    vcf = VCF(VCF_FILE)
    results = []

    with open(SNP_FILE) as f:
        snps = [line.strip() for line in f if line.strip()]

    for snp_id in snps:
        chrom, pos, ref, alt = parse_snp_id(snp_id)

        # Define adjacent non-overlapping windows:
        # left   = [pos-WINDOW_SIZE, pos-1]
        # center = [pos, pos+WINDOW_SIZE-1]
        # right  = [pos+WINDOW_SIZE, pos+2*WINDOW_SIZE-1]
        #
        # This makes all three windows the same nominal size.
        left_start = max(1, pos - WINDOW_SIZE)
        left_end = pos - 1

        center_start = pos
        center_end = pos + WINDOW_SIZE - 1

        right_start = pos + WINDOW_SIZE
        right_end = pos + (2 * WINDOW_SIZE) - 1

        left_stats = window_pi(vcf, chrom, left_start, left_end)
        center_stats = window_pi(vcf, chrom, center_start, center_end, exclude_pos=pos)
        right_stats = window_pi(vcf, chrom, right_start, right_end)

        flank_mean_pi = mean_ignore_none([left_stats["pi"], right_stats["pi"]])

        if center_stats["pi"] is not None and flank_mean_pi is not None:
            delta_pi = center_stats["pi"] - flank_mean_pi
            pi_ratio = center_stats["pi"] / flank_mean_pi if flank_mean_pi != 0 else None
        else:
            delta_pi = None
            pi_ratio = None

        results.append({
            "snp_id": snp_id,
            "chrom": chrom,
            "pos": pos,
            "ref": ref,
            "alt": alt,

            "left_start": left_start,
            "left_end": left_end,
            "left_window_size": max(0, left_end - left_start + 1),
            "left_pi": left_stats["pi"],
            "left_pi_sum": left_stats["pi_sum"],
            "left_n_variant_sites": left_stats["n_variant_sites"],
            "left_n_sites_used": left_stats["n_sites_used"],

            "center_start": center_start,
            "center_end": center_end,
            "center_window_size": center_end - center_start + 1,
            "center_pi": center_stats["pi"],
            "center_pi_sum": center_stats["pi_sum"],
            "center_n_variant_sites": center_stats["n_variant_sites"],
            "center_n_sites_used": center_stats["n_sites_used"],

            "right_start": right_start,
            "right_end": right_end,
            "right_window_size": right_end - right_start + 1,
            "right_pi": right_stats["pi"],
            "right_pi_sum": right_stats["pi_sum"],
            "right_n_variant_sites": right_stats["n_variant_sites"],
            "right_n_sites_used": right_stats["n_sites_used"],

            "flank_mean_pi": flank_mean_pi,
            "delta_pi": delta_pi,
            "pi_ratio": pi_ratio,
        })

    df = pd.DataFrame(results)
    df.to_csv(OUT_FILE, sep="\t", index=False)

    print(f"Wrote: {OUT_FILE}")
    print(df.head())

    # Quick summary
    if "delta_pi" in df.columns:
        print("\nSummary of delta_pi (center - mean flanks):")
        print(df["delta_pi"].describe())

        n_lower = (df["delta_pi"] < 0).sum()
        n_higher = (df["delta_pi"] > 0).sum()
        n_equal = (df["delta_pi"] == 0).sum()

        print(f"\nCenter lower than flanks:  {n_lower}")
        print(f"Center higher than flanks: {n_higher}")
        print(f"Center equal to flanks:    {n_equal}")


if __name__ == "__main__":
    main()