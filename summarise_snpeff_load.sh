#!/usr/bin/env bash

set -euo pipefail

usage() {
  cat <<'EOF'
Usage:
  summarize_snpeff_load.sh --vcf FILE --groups FILE [options]

Required arguments:
  --vcf FILE           Annotated VCF/VCF.GZ from snpEff
  --groups FILE        Two-column sample group file: SAMPLE<TAB>GROUP
                       GROUP must be ingroup or outgroup

Optional arguments:
  --out FILE           Output TSV
                       Default: input with .load_summary.tsv suffix
  --python BIN         Python binary to use
                       Default: python3 from PATH
  --jackknife-bins N   Compute delete-a-group jackknife uncertainty using N gene bins
                       Adds SE and 95% CI columns for normalised metrics
                       Default: 0 (disabled)
  --allow-missing-outgroup
                       Allow missing outgroup genotypes; all called outgroups must still be 0/0
  --force              Overwrite output if it already exists
  --skip-if-exists     Exit cleanly if output already exists
  -h, --help           Show this help

Output columns:
  sample
  modifier_realised      modifier_masked
  synonymous_realised    synonymous_masked
  nonsynonymous_realised nonsynonymous_masked
  lof_realised           lof_masked
  synonymous_total
  lof_realised_over_synonymous_total
  lof_masked_over_synonymous_total
  nonsynonymous_realised_over_synonymous_total
  nonsynonymous_masked_over_synonymous_total
  other_realised         other_masked

Notes:
  - Sites are counted only if all outgroup individuals are fixed for the reference allele.
  - "Realised" means the derived allele is homozygous in the ingroup individual.
  - "Masked" means the derived allele is present in a heterozygous state.
  - For multiallelic sites, each derived ALT allele carried by a sample is considered separately.
  - Variant class follows your previous definitions exactly:
    loss_of_function = stop_gained, frameshift_variant, splice_acceptor_variant, splice_donor_variant
    nonsynonymous = missense_variant
    synonymous = synonymous_variant, start_retained_variant, stop_retained_variant
    modifier = all remaining MODIFIER-impact annotations
  - If jackknife bins are requested, genes are ordered by first appearance in the VCF and split
    into roughly equal-sized bins. SE and 95% CI are reported for normalised metrics only,
    using estimate +/- 1.96 * SE.

Example groups file:
  Polar_bear_UAmap      outgroup
  American_black_bear_UAmap outgroup
  SW_01_UAmap           ingroup
  SW_02_UAmap           ingroup

Example:
  bash summarize_snpeff_load.sh \
    --vcf /path/sample.anno.vcf.gz \
    --groups /path/groups.tsv
EOF
}

die() {
  printf 'Error: %s\n' "$1" >&2
  exit 1
}

require_file() {
  local path="$1"
  [[ -f "$path" ]] || die "File not found: $path"
}

VCF=""
GROUPS_FILE=""
OUT=""
PYTHON_BIN=""
JACKKNIFE_BINS=0
ALLOW_MISSING_OUTGROUP=0
FORCE=0
SKIP_IF_EXISTS=0

while [[ $# -gt 0 ]]; do
  case "$1" in
    --vcf)
      VCF="${2:-}"
      shift 2
      ;;
    --groups)
      GROUPS_FILE="${2:-}"
      shift 2
      ;;
    --out)
      OUT="${2:-}"
      shift 2
      ;;
    --python)
      PYTHON_BIN="${2:-}"
      shift 2
      ;;
    --jackknife-bins)
      JACKKNIFE_BINS="${2:-}"
      shift 2
      ;;
    --allow-missing-outgroup)
      ALLOW_MISSING_OUTGROUP=1
      shift
      ;;
    --force)
      FORCE=1
      shift
      ;;
    --skip-if-exists)
      SKIP_IF_EXISTS=1
      shift
      ;;
    -h|--help)
      usage
      exit 0
      ;;
    *)
      die "Unknown argument: $1"
      ;;
  esac
done

[[ -n "$VCF" ]] || die "--vcf is required"
[[ -n "$GROUPS_FILE" ]] || die "--groups is required"
[[ "$JACKKNIFE_BINS" =~ ^[0-9]+$ ]] || die "--jackknife-bins must be a non-negative integer"

require_file "$VCF"
require_file "$GROUPS_FILE"

if [[ -z "$PYTHON_BIN" ]]; then
  PYTHON_BIN="$(command -v python3 || true)"
fi
[[ -n "$PYTHON_BIN" ]] || die "python3 not found in PATH and --python was not provided"
require_file "$PYTHON_BIN"

if [[ -z "$OUT" ]]; then
  case "$VCF" in
    *.vcf.gz)
      OUT="${VCF%.vcf.gz}.load_summary.tsv"
      ;;
    *.vcf)
      OUT="${VCF%.vcf}.load_summary.tsv"
      ;;
    *)
      OUT="${VCF}.load_summary.tsv"
      ;;
  esac
fi

if [[ -f "$OUT" ]]; then
  if [[ "$SKIP_IF_EXISTS" -eq 1 && "$FORCE" -ne 1 ]]; then
    printf 'Output already exists, skipping: %s\n' "$OUT"
    exit 0
  fi
  if [[ "$FORCE" -ne 1 ]]; then
    die "Output already exists: $OUT (use --force or --skip-if-exists)"
  fi
fi

printf 'Summarising realised and masked load from snpEff VCF\n'
printf '  vcf    = %s\n' "$VCF"
printf '  groups = %s\n' "$GROUPS_FILE"
printf '  out    = %s\n' "$OUT"
printf '  bins   = %s\n' "$JACKKNIFE_BINS"

"$PYTHON_BIN" - "$VCF" "$GROUPS_FILE" "$OUT" "$ALLOW_MISSING_OUTGROUP" "$JACKKNIFE_BINS" <<'PYCODE'
import csv
import gzip
import math
import sys
from collections import defaultdict

vcf_path, groups_path, out_path, allow_missing_flag, jackknife_bins_arg = sys.argv[1:]
allow_missing_outgroup = allow_missing_flag == "1"
requested_jackknife_bins = int(jackknife_bins_arg)

LOF_EFFECTS = {
    "frameshift_variant",
    "stop_gained",
    "splice_acceptor_variant",
    "splice_donor_variant",
}

NONSYNONYMOUS_EFFECTS = {
    "missense_variant",
}

SYNONYMOUS_EFFECTS = {
    "synonymous_variant",
    "start_retained_variant",
    "stop_retained_variant",
}

IMPACT_RANK = {
    "HIGH": 4,
    "MODERATE": 3,
    "LOW": 2,
    "MODIFIER": 1,
}

CATEGORY_ORDER = ["modifier", "synonymous", "nonsynonymous", "lof", "other"]
OUTPUT_COLUMNS = [
    "sample",
    "modifier_realised", "modifier_masked",
    "synonymous_realised", "synonymous_masked",
    "nonsynonymous_realised", "nonsynonymous_masked",
    "lof_realised", "lof_masked",
    "synonymous_total",
    "lof_realised_over_synonymous_total",
    "lof_masked_over_synonymous_total",
    "nonsynonymous_realised_over_synonymous_total",
    "nonsynonymous_masked_over_synonymous_total",
    "other_realised", "other_masked",
]

JACKKNIFE_METRICS = [
    "lof_realised_over_synonymous_total",
    "lof_masked_over_synonymous_total",
    "nonsynonymous_realised_over_synonymous_total",
    "nonsynonymous_masked_over_synonymous_total",
]

for metric in JACKKNIFE_METRICS:
    OUTPUT_COLUMNS.extend([
        f"{metric}_jackknife_bins",
        f"{metric}_jackknife_se",
        f"{metric}_jackknife_ci_lower",
        f"{metric}_jackknife_ci_upper",
    ])


def open_text(path):
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "rt")


def normalize_group(value):
    v = value.strip().lower()
    if v in {"ingroup", "in"}:
        return "ingroup"
    if v in {"outgroup", "out"}:
        return "outgroup"
    return None


def parse_groups(path):
    ingroup = []
    outgroup = []
    seen = set()
    with open(path, "rt") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split()
            if len(parts) < 2:
                continue
            sample, group = parts[0], normalize_group(parts[1])
            if sample.lower() == "sample" and parts[1].strip().lower() == "group":
                continue
            if group is None:
                raise SystemExit(f"Invalid group '{parts[1]}' in {path}. Use ingroup or outgroup.")
            if sample in seen:
                raise SystemExit(f"Sample appears more than once in groups file: {sample}")
            seen.add(sample)
            if group == "ingroup":
                ingroup.append(sample)
            else:
                outgroup.append(sample)
    if not ingroup:
        raise SystemExit("No ingroup samples were found in the groups file.")
    if not outgroup:
        raise SystemExit("No outgroup samples were found in the groups file.")
    return ingroup, outgroup


def parse_info_ann(info_field):
    for entry in info_field.split(";"):
        if entry.startswith("ANN="):
            value = entry[4:]
            return [] if value == "." else value.split(",")
    return []


def classify_ann_entries(entries):
    best_category = "other"
    best_rank = (-1, -1)
    best_gene = None
    for ann in entries:
        fields = ann.split("|")
        effects_field = fields[1] if len(fields) > 1 else ""
        impact_field = fields[2] if len(fields) > 2 else ""
        gene_name = fields[3] if len(fields) > 3 else ""
        gene_id = fields[4] if len(fields) > 4 else ""
        effects = set(filter(None, effects_field.split("&")))
        impact_score = IMPACT_RANK.get(impact_field, 0)

        if effects & LOF_EFFECTS:
            category = "lof"
            priority = 5
        elif effects & NONSYNONYMOUS_EFFECTS:
            category = "nonsynonymous"
            priority = 4
        elif effects & SYNONYMOUS_EFFECTS:
            category = "synonymous"
            priority = 3
        elif impact_field == "MODIFIER":
            category = "modifier"
            priority = 2
        else:
            category = "other"
            priority = 1

        rank = (priority, impact_score)
        if rank > best_rank:
            best_rank = rank
            best_category = category
            best_gene = gene_id or gene_name or None

    return best_category, best_gene


def alt_details_from_ann(info_field, alts, chrom, pos):
    ann_entries = parse_info_ann(info_field)
    by_alt = defaultdict(list)
    for ann in ann_entries:
        allele = ann.split("|", 1)[0]
        by_alt[allele].append(ann)

    details = {}
    for alt in alts:
        # Some records may lack a matching ANN entry for a particular ALT.
        category, gene = classify_ann_entries(by_alt.get(alt, []))
        gene_key = gene or f"{chrom}:{pos}:{alt}"
        details[alt] = {"category": category, "gene": gene_key}
    return details


def genotype_alleles(sample_field):
    gt = sample_field.split(":", 1)[0]
    if gt in {".", "./.", ".|."}:
        return None
    sep = "/" if "/" in gt else "|"
    alleles = gt.split(sep)
    if any(a == "." for a in alleles):
        return None
    return alleles


ingroup_names, outgroup_names = parse_groups(groups_path)
counts = {
    sample: {f"{category}_{state}": 0 for category in CATEGORY_ORDER for state in ("realised", "masked")}
    for sample in ingroup_names
}
gene_counts = defaultdict(lambda: {
    sample: {f"{category}_{state}": 0 for category in CATEGORY_ORDER for state in ("realised", "masked")}
    for sample in ingroup_names
})
gene_order = []
seen_genes = set()

header_samples = None
ingroup_idx = {}
outgroup_idx = {}

with open_text(vcf_path) as handle:
    for raw_line in handle:
        if raw_line.startswith("##"):
            continue
        if raw_line.startswith("#CHROM"):
            fields = raw_line.rstrip("\n").split("\t")
            header_samples = fields[9:]
            sample_to_idx = {sample: idx for idx, sample in enumerate(header_samples)}

            missing = [s for s in ingroup_names + outgroup_names if s not in sample_to_idx]
            if missing:
                raise SystemExit("Samples in groups file not found in VCF header: " + ", ".join(missing))

            ingroup_idx = {sample: sample_to_idx[sample] for sample in ingroup_names}
            outgroup_idx = {sample: sample_to_idx[sample] for sample in outgroup_names}
            continue

        if raw_line.startswith("#"):
            continue

        fields = raw_line.rstrip("\n").split("\t")
        if len(fields) < 10:
            continue

        chrom = fields[0]
        pos = fields[1]
        ref = fields[3]
        alts = fields[4].split(",")
        info = fields[7]
        sample_fields = fields[9:]

        outgroup_ok = True
        called_outgroups = 0
        for sample, idx in outgroup_idx.items():
            alleles = genotype_alleles(sample_fields[idx])
            if alleles is None:
                if allow_missing_outgroup:
                    continue
                outgroup_ok = False
                break
            called_outgroups += 1
            if any(allele != "0" for allele in alleles):
                outgroup_ok = False
                break

        if not outgroup_ok:
            continue
        if allow_missing_outgroup and called_outgroups == 0:
            continue

        alt_details = alt_details_from_ann(info, alts, chrom, pos)

        for sample, idx in ingroup_idx.items():
            alleles = genotype_alleles(sample_fields[idx])
            if alleles is None:
                continue

            ploidy = len(alleles)
            allele_counts = defaultdict(int)
            for allele in alleles:
                allele_counts[allele] += 1

            for allele_code, dosage in allele_counts.items():
                if allele_code == "0":
                    continue
                try:
                    alt_index = int(allele_code) - 1
                except ValueError:
                    continue
                if alt_index < 0 or alt_index >= len(alts):
                    continue

                alt = alts[alt_index]
                detail = alt_details.get(alt, {"category": "other", "gene": f"{chrom}:{pos}:{alt}"})
                category = detail["category"]
                gene_key = detail["gene"]
                if gene_key not in seen_genes:
                    seen_genes.add(gene_key)
                    gene_order.append(gene_key)
                state = "realised" if dosage == ploidy else "masked"
                counts[sample][f"{category}_{state}"] += 1
                gene_counts[gene_key][sample][f"{category}_{state}"] += 1


def synonymous_total(count_dict):
    return count_dict["synonymous_realised"] + 0.5 * count_dict["synonymous_masked"]


def build_ratio_fields(count_dict):
    syn_total = synonymous_total(count_dict)
    row = {"synonymous_total": syn_total}
    if syn_total > 0:
        row["lof_realised_over_synonymous_total"] = count_dict["lof_realised"] / syn_total
        row["lof_masked_over_synonymous_total"] = count_dict["lof_masked"] / syn_total
        row["nonsynonymous_realised_over_synonymous_total"] = count_dict["nonsynonymous_realised"] / syn_total
        row["nonsynonymous_masked_over_synonymous_total"] = count_dict["nonsynonymous_masked"] / syn_total
    else:
        row["lof_realised_over_synonymous_total"] = "NA"
        row["lof_masked_over_synonymous_total"] = "NA"
        row["nonsynonymous_realised_over_synonymous_total"] = "NA"
        row["nonsynonymous_masked_over_synonymous_total"] = "NA"
    return row


def zero_count_dict():
    return {f"{category}_{state}": 0 for category in CATEGORY_ORDER for state in ("realised", "masked")}


def add_counts_in_place(target, source):
    for key, value in source.items():
        target[key] += value


def subtract_counts(base, minus):
    result = dict(base)
    for key, value in minus.items():
        result[key] -= value
    return result


def format_output_value(value):
    if isinstance(value, float):
        return f"{value:.6f}"
    return value


jackknife_results = {
    sample: {metric: None for metric in JACKKNIFE_METRICS}
    for sample in ingroup_names
}

jackknife_bins_used = 0
if requested_jackknife_bins > 0 and gene_order:
    jackknife_bins_used = min(requested_jackknife_bins, len(gene_order))
    bins = [[] for _ in range(jackknife_bins_used)]
    for idx, gene in enumerate(gene_order):
        bins[idx % jackknife_bins_used].append(gene)

    bin_counts = []
    for genes in bins:
        per_sample = {sample: zero_count_dict() for sample in ingroup_names}
        for gene in genes:
            gene_sample_counts = gene_counts[gene]
            for sample in ingroup_names:
                add_counts_in_place(per_sample[sample], gene_sample_counts[sample])
        bin_counts.append(per_sample)

    for sample in ingroup_names:
        theta_hat = build_ratio_fields(counts[sample])
        metric_values = {metric: [] for metric in JACKKNIFE_METRICS}
        for bin_idx in range(jackknife_bins_used):
            reduced_counts = subtract_counts(counts[sample], bin_counts[bin_idx][sample])
            reduced_ratios = build_ratio_fields(reduced_counts)
            for metric in JACKKNIFE_METRICS:
                value = reduced_ratios[metric]
                if value == "NA":
                    metric_values[metric] = None
                    continue
                if metric_values[metric] is not None:
                    metric_values[metric].append(value)

        for metric in JACKKNIFE_METRICS:
            values = metric_values[metric]
            theta = theta_hat[metric]
            if theta == "NA" or not values or len(values) < 2:
                continue
            mean_theta = sum(values) / len(values)
            variance = ((len(values) - 1) / len(values)) * sum((x - mean_theta) ** 2 for x in values)
            se = math.sqrt(max(variance, 0.0))
            jackknife_results[sample][metric] = {
                "bins": len(values),
                "se": se,
                "ci_lower": theta - 1.96 * se,
                "ci_upper": theta + 1.96 * se,
            }

with open(out_path, "wt", newline="") as out_handle:
    writer = csv.DictWriter(out_handle, fieldnames=OUTPUT_COLUMNS, delimiter="\t")
    writer.writeheader()
    for sample in ingroup_names:
        row = {"sample": sample}
        row.update(counts[sample])
        row.update(build_ratio_fields(counts[sample]))
        for metric in JACKKNIFE_METRICS:
            result = jackknife_results[sample][metric]
            if result is None:
                row[f"{metric}_jackknife_bins"] = jackknife_bins_used if jackknife_bins_used > 0 else "NA"
                row[f"{metric}_jackknife_se"] = "NA"
                row[f"{metric}_jackknife_ci_lower"] = "NA"
                row[f"{metric}_jackknife_ci_upper"] = "NA"
            else:
                row[f"{metric}_jackknife_bins"] = result["bins"]
                row[f"{metric}_jackknife_se"] = result["se"]
                row[f"{metric}_jackknife_ci_lower"] = result["ci_lower"]
                row[f"{metric}_jackknife_ci_upper"] = result["ci_upper"]
        formatted_row = {key: format_output_value(value) for key, value in row.items()}
        writer.writerow(formatted_row)
PYCODE

printf 'Finished: %s\n' "$OUT"
