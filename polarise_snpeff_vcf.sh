#!/usr/bin/env bash
#SBATCH --job-name=snpeff_polarise
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=24:00:00
#SBATCH --output=slurm-%j.out

set -euo pipefail

usage() {
  cat <<'EOF'
Usage:
  polarise_snpeff_vcf.sh --vcf FILE --groups FILE [options]

Required arguments:
  --vcf FILE           Input VCF/VCF.GZ to polarise before snpEff
  --groups FILE        Two-column sample group file: SAMPLE<TAB>GROUP
                       GROUP must be ingroup or outgroup

Optional arguments:
  --out FILE           Output polarised VCF.GZ
                       Default: input with .polarised.vcf.gz suffix
  --python BIN         Python binary to use
                       Default: python3 from PATH
  --bgzip BIN          bgzip binary to use
                       Default: bgzip from PATH
  --allow-missing-outgroup
                       Allow missing outgroup genotypes; all called outgroups
                       must still agree on the same homozygous allele
  --force              Overwrite output if it already exists
  --skip-if-exists     Exit cleanly if output already exists
  -h, --help           Show this help

What it does:
  1. Reads a mixed VCF and keeps biallelic SNPs only
  2. Uses outgroup samples to infer the ancestral allele at each site
  3. Keeps only sites where all called outgroups agree on the same homozygous allele
  4. If the ancestral allele is ALT, swaps REF and ALT and flips all genotypes
  5. Writes a simplified load-only VCF with FORMAT=GT

Important notes:
  - This script is meant for the load workflow only
  - Indels and multiallelic sites are skipped automatically
  - It drops allele-dependent FORMAT and INFO fields after polarisation, because
    values such as AD, PL, AC, and AF are no longer trustworthy after flipping
  - The output should be re-annotated with snpEff before running the load summary

Example:
  sbatch polarise_snpeff_vcf.sh \
    --vcf /path/input.vcf.gz \
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
BGZIP_BIN=""
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
    --bgzip)
      BGZIP_BIN="${2:-}"
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

require_file "$VCF"
require_file "$GROUPS_FILE"

if [[ -z "$PYTHON_BIN" ]]; then
  PYTHON_BIN="$(command -v python3 || true)"
fi
[[ -n "$PYTHON_BIN" ]] || die "python3 not found in PATH and --python was not provided"
require_file "$PYTHON_BIN"

if [[ -z "$BGZIP_BIN" ]]; then
  BGZIP_BIN="$(command -v bgzip || true)"
fi
[[ -n "$BGZIP_BIN" ]] || die "bgzip not found in PATH and --bgzip was not provided"
require_file "$BGZIP_BIN"

if [[ -z "$OUT" ]]; then
  case "$VCF" in
    *.vcf.gz)
      OUT="${VCF%.vcf.gz}.polarised.vcf.gz"
      ;;
    *.vcf)
      OUT="${VCF%.vcf}.polarised.vcf.gz"
      ;;
    *)
      OUT="${VCF}.polarised.vcf.gz"
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

TMP_VCF="${OUT%.gz}.tmp.vcf"
TMP_OUT="${OUT}.tmp"
rm -f "$TMP_VCF" "$TMP_OUT"

printf 'Polarising load VCF using outgroups\n'
printf '  vcf    = %s\n' "$VCF"
printf '  groups = %s\n' "$GROUPS_FILE"
printf '  out    = %s\n' "$OUT"

"$PYTHON_BIN" - "$VCF" "$GROUPS_FILE" "$ALLOW_MISSING_OUTGROUP" > "$TMP_VCF" <<'PYCODE'
import gzip
import sys
from collections import Counter

vcf_path, groups_path, allow_missing_flag = sys.argv[1:]
allow_missing_outgroup = allow_missing_flag == "1"


def die(message):
    print(f"Error: {message}", file=sys.stderr)
    sys.exit(1)


def parse_groups(path):
    groups = {}
    with open(path, "r", encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, start=1):
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            fields = line.split()
            if len(fields) != 2:
                die(f"Group file line {line_number} does not have exactly two columns")
            sample, group = fields
            if group not in {"ingroup", "outgroup"}:
                die(f"Group file line {line_number} has invalid group '{group}'")
            groups[sample] = group
    return groups


def open_vcf(path):
    if path.endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8")
    return open(path, "r", encoding="utf-8")


def parse_gt(sample_field):
    gt = sample_field.split(":", 1)[0]
    if gt in {".", "./.", ".|."}:
        return None
    if "/" in gt:
        sep = "/"
    elif "|" in gt:
        sep = "|"
    else:
        return "unsupported"
    alleles = gt.split(sep)
    if len(alleles) != 2:
        return "unsupported"
    if any(allele == "." for allele in alleles):
        return None
    if any(allele not in {"0", "1"} for allele in alleles):
        return "unsupported"
    return tuple(int(allele) for allele in alleles), sep


def flip_gt(parsed_gt, flip):
    if parsed_gt is None:
        return "./."
    (a1, a2), sep = parsed_gt
    if flip:
        a1 = 1 - a1
        a2 = 1 - a2
    return f"{a1}{sep}{a2}"


groups = parse_groups(groups_path)
outgroup_names = {sample for sample, group in groups.items() if group == "outgroup"}
if not outgroup_names:
    die("No outgroup samples were found in the groups file")

stats = Counter()
header_written = False
gt_header_present = False

with open_vcf(vcf_path) as handle:
    for raw_line in handle:
        if raw_line.startswith("##"):
            if raw_line.startswith("##FORMAT=<ID=GT,"):
                gt_header_present = True
            print(raw_line, end="")
            continue

        if raw_line.startswith("#CHROM"):
            fields = raw_line.rstrip("\n").split("\t")
            if len(fields) < 10:
                die("VCF does not contain sample columns")
            sample_names = fields[9:]
            sample_to_index = {sample: idx for idx, sample in enumerate(sample_names)}

            missing_outgroups = sorted(outgroup_names - set(sample_to_index))
            if missing_outgroups:
                die(
                    "Outgroup samples missing from VCF header: "
                    + ", ".join(missing_outgroups)
                )

            if not gt_header_present:
                print('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">')
            print('##polarise_snpeff_vcf=Ancestral allele set from agreeing homozygous outgroups; output reduced to GT-only for load analyses')
            print(raw_line, end="")
            outgroup_indices = [sample_to_index[name] for name in sample_names if name in outgroup_names]
            header_written = True
            continue

        if not header_written:
            die("VCF header was not found before data rows")

        stats["total_sites"] += 1

        fields = raw_line.rstrip("\n").split("\t")
        if len(fields) < 10:
            stats["skipped_malformed"] += 1
            continue

        chrom, pos, site_id, ref, alt, qual, filt = fields[:7]
        sample_fields = fields[9:]

        if "," in alt:
            stats["skipped_multiallelic"] += 1
            continue

        if len(ref) != 1 or len(alt) != 1:
            stats["skipped_non_snp"] += 1
            continue

        parsed_sample_gts = []
        unsupported_site = False
        for sample_field in sample_fields:
            parsed = parse_gt(sample_field)
            if parsed == "unsupported":
                unsupported_site = True
                break
            parsed_sample_gts.append(parsed)

        if unsupported_site:
            stats["skipped_unsupported_gt"] += 1
            continue

        outgroup_called_alleles = []
        skip_reason = None

        for idx in outgroup_indices:
            parsed_gt = parsed_sample_gts[idx]
            if parsed_gt is None:
                if allow_missing_outgroup:
                    continue
                skip_reason = "skipped_missing_outgroup"
                break

            alleles, _sep = parsed_gt
            if alleles[0] != alleles[1]:
                skip_reason = "skipped_heterozygous_outgroup"
                break
            outgroup_called_alleles.append(alleles[0])

        if skip_reason is not None:
            stats[skip_reason] += 1
            continue

        if not outgroup_called_alleles:
            stats["skipped_no_called_outgroup"] += 1
            continue

        if len(set(outgroup_called_alleles)) != 1:
            stats["skipped_discordant_outgroup"] += 1
            continue

        ancestral_allele = outgroup_called_alleles[0]
        flip = ancestral_allele == 1

        out_ref = alt if flip else ref
        out_alt = ref if flip else alt
        out_samples = [flip_gt(parsed_gt, flip) for parsed_gt in parsed_sample_gts]

        out_fields = [
            chrom,
            pos,
            site_id,
            out_ref,
            out_alt,
            qual,
            filt,
            ".",
            "GT",
            *out_samples,
        ]
        print("\t".join(out_fields))

        if flip:
            stats["kept_flipped_sites"] += 1
        else:
            stats["kept_ref_ancestral_sites"] += 1
        stats["kept_sites"] += 1

print("Polarisation summary", file=sys.stderr)
for key in [
    "total_sites",
    "kept_sites",
    "kept_ref_ancestral_sites",
    "kept_flipped_sites",
    "skipped_multiallelic",
    "skipped_non_snp",
    "skipped_unsupported_gt",
    "skipped_missing_outgroup",
    "skipped_no_called_outgroup",
    "skipped_heterozygous_outgroup",
    "skipped_discordant_outgroup",
    "skipped_malformed",
]:
    print(f"  {key} = {stats.get(key, 0)}", file=sys.stderr)
PYCODE

"$BGZIP_BIN" -c "$TMP_VCF" > "$TMP_OUT"
mv "$TMP_OUT" "$OUT"
rm -f "$TMP_VCF"

printf 'Finished: %s\n' "$OUT"
