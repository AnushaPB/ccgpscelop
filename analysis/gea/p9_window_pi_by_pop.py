#!/usr/bin/env python3

#!/usr/bin/env python3

from cyvcf2 import VCF
import pandas as pd
from collections import Counter
import os
import re

# -----------------------------
# input files
# -----------------------------
VCF_FILE = "../../data/ccgp_data/58-Sceloporus_complete_coords_annotated_final_samples.vcf.gz"
SNP_FILE = "outputs/bio1ndvi_gea_genes_nonsyn_ids.txt"
POP_DIR = "../../analysis/admixture/outputs"
OUT_FILE = "outputs/bio1ndvi_gea_genes_nonsyn_flanking_pi_by_pop.tsv"

WINDOW_SIZE = 50000

"""
Compute nucleotide diversity (π) in three adjacent windows (left, center, right)
around each SNP within each population. The focal SNP itself is excluded from
the center window calculation.
"""

# -----------------------------
# helper functions
# -----------------------------
def parse_snp_id(s):
    parts = s.strip().split("_")
    if len(parts) < 4:
        raise ValueError(f"Could not parse SNP ID: {s}")
    chrom = parts[0]
    pos = int(parts[1])
    ref = parts[2]
    alt = parts[3]
    return chrom, pos, ref, alt


def load_population_files(pop_dir):
    """
    Read only files named exactly:
        k8_pop1.txt, k8_pop2.txt, ..., k8_pop8.txt

    Ignores files such as:
        k8_pop4_distinguished_individuals.txt
    """
    pattern = re.compile(r"^k8_pop([1-8])\.txt$")
    populations = {}

    for fname in sorted(os.listdir(pop_dir)):
        match = pattern.match(fname)
        if not match:
            continue

        pop_name = os.path.splitext(fname)[0]
        filepath = os.path.join(pop_dir, fname)

        with open(filepath) as f:
            samples = [line.strip() for line in f if line.strip()]

        populations[pop_name] = samples

    if not populations:
        raise FileNotFoundError(
            f"No population files matching exactly k8_pop1.txt ... k8_pop8.txt found in {pop_dir}"
        )

    return populations


def map_population_samples_to_vcf_indices(vcf_samples, populations):
    sample_to_index = {sample: i for i, sample in enumerate(vcf_samples)}

    pop_to_indices = {}
    pop_to_missing = {}

    for pop_name, samples in populations.items():
        indices = []
        missing = []

        for s in samples:
            if s in sample_to_index:
                indices.append(sample_to_index[s])
            else:
                missing.append(s)

        pop_to_indices[pop_name] = indices
        pop_to_missing[pop_name] = missing

    return pop_to_indices, pop_to_missing


def site_pi_from_genotypes_subset(genotypes, sample_indices):
    """
    Calculate per-site nucleotide diversity from diploid genotypes,
    using only the samples in sample_indices.
    """
    allele_counts = Counter()

    for idx in sample_indices:
        gt = genotypes[idx]
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


def window_pi_by_population(vcf, chrom, start, end, sample_indices, exclude_pos=None):
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
            pi_site = site_pi_from_genotypes_subset(var.genotypes, sample_indices)
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
    if not os.path.exists(POP_DIR):
        raise FileNotFoundError(f"Population directory not found: {POP_DIR}")

    vcf = VCF(VCF_FILE)
    vcf_samples = list(vcf.samples)

    populations = load_population_files(POP_DIR)
    pop_to_indices, pop_to_missing = map_population_samples_to_vcf_indices(vcf_samples, populations)

    print("\nPopulation sample matching summary:")
    for pop_name in sorted(populations):
        n_total = len(populations[pop_name])
        n_found = len(pop_to_indices[pop_name])
        n_missing = len(pop_to_missing[pop_name])
        print(f"{pop_name}: found {n_found}/{n_total} samples in VCF")
        if n_missing > 0:
            print(f"  Missing examples: {pop_to_missing[pop_name][:5]}")

    with open(SNP_FILE) as f:
        snps = [line.strip() for line in f if line.strip()]

    results = []

    for snp_id in snps:
        chrom, pos, ref, alt = parse_snp_id(snp_id)

        left_start = max(1, pos - WINDOW_SIZE)
        left_end = pos - 1

        center_start = pos
        center_end = pos + WINDOW_SIZE - 1

        right_start = pos + WINDOW_SIZE
        right_end = pos + (2 * WINDOW_SIZE) - 1

        for pop_name, sample_indices in pop_to_indices.items():
            if len(sample_indices) == 0:
                continue

            left_stats = window_pi_by_population(
                vcf, chrom, left_start, left_end, sample_indices
            )
            center_stats = window_pi_by_population(
                vcf, chrom, center_start, center_end, sample_indices, exclude_pos=pos
            )
            right_stats = window_pi_by_population(
                vcf, chrom, right_start, right_end, sample_indices
            )

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
                "population": pop_name,
                "n_samples_in_pop": len(sample_indices),

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

    print(f"\nWrote: {OUT_FILE}")
    print(df.head())

    print("\nSummary by population:")
    summary = (
        df.groupby("population", dropna=False)["delta_pi"]
        .agg(["count", "mean", "median", "std"])
        .reset_index()
    )
    print(summary.to_string(index=False))


if __name__ == "__main__":
    main()