#!/usr/bin/env bash

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
  --tabix BIN          tabix binary to use for chromosome batching
                       Default: tabix from PATH
  -c, --chromosomes-per-batch N
                       Process N chromosomes at a time using tabix regions
                       Requires a bgzipped VCF; index is created automatically if missing
                       If N > 0 and --batch-index is not set, the script processes
                       all chromosome batches internally in parallel
  --batch-index N      Zero-based chromosome batch index to process
                       Optional: mainly useful for advanced manual splitting
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
  6. Optionally splits chromosomes into internal batches and runs them in parallel
     using the CPUs allocated to the job

Important notes:
  - This script is meant for the load workflow only
  - Indels and multiallelic sites are skipped automatically
  - It drops allele-dependent FORMAT and INFO fields after polarisation, because
    values such as AD, PL, AC, and AF are no longer trustworthy after flipping
  - The output should be re-annotated with snpEff before running the load summary
  - In chromosome-batch mode, a tabix index is created automatically if missing
  - When running internal chromosome batches, the script uses SLURM_CPUS_PER_TASK
    if available, otherwise it runs one batch at a time

Example:
  sbatch polarise_snpeff_vcf.sh \
    --vcf /path/input.vcf.gz \
    --groups /path/groups.tsv

  sbatch --cpus-per-task=10 polarise_snpeff_vcf.sh \
    --vcf /path/input.vcf.gz \
    --groups /path/groups.tsv \
    -c 10
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
TABIX_BIN=""
CHROMS_PER_BATCH=0
BATCH_INDEX=""
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
    --tabix)
      TABIX_BIN="${2:-}"
      shift 2
      ;;
    -c|--chromosomes-per-batch)
      CHROMS_PER_BATCH="${2:-}"
      shift 2
      ;;
    --batch-index)
      BATCH_INDEX="${2:-}"
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
[[ "$CHROMS_PER_BATCH" =~ ^[0-9]+$ ]] || die "--chromosomes-per-batch must be a non-negative integer"
[[ -z "$BATCH_INDEX" || "$BATCH_INDEX" =~ ^[0-9]+$ ]] || die "--batch-index must be a non-negative integer"

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

REGION_MODE=0
INTERNAL_BATCH_MODE=0
if [[ "$CHROMS_PER_BATCH" -gt 0 ]]; then
  REGION_MODE=1
  if [[ -z "$TABIX_BIN" ]]; then
    TABIX_BIN="$(command -v tabix || true)"
  fi
  [[ -n "$TABIX_BIN" ]] || die "tabix not found in PATH and --tabix was not provided"
  require_file "$TABIX_BIN"
  [[ "$VCF" == *.gz ]] || die "--chromosomes-per-batch requires a bgzipped VCF"
  if [[ ! -f "${VCF}.tbi" && ! -f "${VCF}.csi" ]]; then
    printf 'No tabix index found, creating one for: %s\n' "$VCF"
    "$TABIX_BIN" -f -p vcf "$VCF"
  fi

  if [[ -z "$BATCH_INDEX" && -n "${SLURM_ARRAY_TASK_ID:-}" ]]; then
    BATCH_INDEX="${SLURM_ARRAY_TASK_ID}"
  fi
fi

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

SELECTED_CHROMS=()
TOTAL_CHROMS=0
TOTAL_BATCHES=0
BATCH_SUFFIX=""
MAX_PARALLEL_JOBS="${SLURM_CPUS_PER_TASK:-1}"
[[ "$MAX_PARALLEL_JOBS" =~ ^[0-9]+$ ]] || MAX_PARALLEL_JOBS=1
(( MAX_PARALLEL_JOBS >= 1 )) || MAX_PARALLEL_JOBS=1

if [[ "$REGION_MODE" -eq 1 ]]; then
  mapfile -t ALL_CHROMS < <("$TABIX_BIN" -l "$VCF")
  TOTAL_CHROMS="${#ALL_CHROMS[@]}"
  [[ "$TOTAL_CHROMS" -gt 0 ]] || die "No chromosomes were found in the tabix index"

  TOTAL_BATCHES=$(( (TOTAL_CHROMS + CHROMS_PER_BATCH - 1) / CHROMS_PER_BATCH ))
  if [[ -n "$BATCH_INDEX" ]]; then
    (( BATCH_INDEX < TOTAL_BATCHES )) || die "--batch-index ${BATCH_INDEX} is out of range for ${TOTAL_BATCHES} batches"

    START_INDEX=$(( BATCH_INDEX * CHROMS_PER_BATCH ))
    END_INDEX=$(( START_INDEX + CHROMS_PER_BATCH ))
    if (( END_INDEX > TOTAL_CHROMS )); then
      END_INDEX=$TOTAL_CHROMS
    fi

    for (( i=START_INDEX; i<END_INDEX; i++ )); do
      SELECTED_CHROMS+=("${ALL_CHROMS[$i]}")
    done

    BATCH_SUFFIX=".batch$(printf '%03d' "$BATCH_INDEX")"
    case "$OUT" in
      *.vcf.gz)
        OUT="${OUT%.vcf.gz}${BATCH_SUFFIX}.vcf.gz"
        ;;
      *.vcf)
        OUT="${OUT%.vcf}${BATCH_SUFFIX}.vcf.gz"
        ;;
      *)
        OUT="${OUT}${BATCH_SUFFIX}.vcf.gz"
        ;;
    esac
  else
    INTERNAL_BATCH_MODE=1
  fi
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
if [[ "$REGION_MODE" -eq 1 && "$INTERNAL_BATCH_MODE" -eq 0 ]]; then
  printf '  batch  = %s of %s\n' "$((BATCH_INDEX + 1))" "$TOTAL_BATCHES"
  printf '  chroms = %s\n' "$(IFS=,; printf '%s' "${SELECTED_CHROMS[*]}")"
elif [[ "$INTERNAL_BATCH_MODE" -eq 1 ]]; then
  printf '  batches = %s\n' "$TOTAL_BATCHES"
  printf '  cpus    = %s\n' "$MAX_PARALLEL_JOBS"
fi

run_python_to_output() {
local vcf_input="$1"
local out_vcf="$2"
"$PYTHON_BIN" - "$vcf_input" "$GROUPS_FILE" "$ALLOW_MISSING_OUTGROUP" > "$out_vcf" <<'PYCODE'
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
}

process_batch() {
  local batch_index="$1"
  shift
  local batch_out="$1"
  shift
  local batch_log="$1"
  shift
  local -a batch_chroms=("$@")

  {
    printf 'Batch %s chromosomes: %s\n' "$batch_index" "$(IFS=,; printf '%s' "${batch_chroms[*]}")"
    "$TABIX_BIN" -h "$VCF" "${batch_chroms[@]}" | run_python_to_output "-" "$batch_out"
  } >"$batch_log" 2>&1
}

if [[ "$INTERNAL_BATCH_MODE" -eq 1 ]]; then
  TMP_DIR="$(mktemp -d)"
  trap 'rm -rf "$TMP_DIR"' EXIT

  pids=()
  batch_outputs=()
  batch_logs=()

  for (( batch_i=0; batch_i<TOTAL_BATCHES; batch_i++ )); do
    start_i=$(( batch_i * CHROMS_PER_BATCH ))
    end_i=$(( start_i + CHROMS_PER_BATCH ))
    if (( end_i > TOTAL_CHROMS )); then
      end_i=$TOTAL_CHROMS
    fi

    batch_chroms=()
    for (( j=start_i; j<end_i; j++ )); do
      batch_chroms+=("${ALL_CHROMS[$j]}")
    done

    batch_out="${TMP_DIR}/batch$(printf '%03d' "$batch_i").vcf"
    batch_log="${TMP_DIR}/batch$(printf '%03d' "$batch_i").log"
    batch_outputs+=("$batch_out")
    batch_logs+=("$batch_log")

    process_batch "$batch_i" "$batch_out" "$batch_log" "${batch_chroms[@]}" &
    pids+=("$!")

    while (( ${#pids[@]} >= MAX_PARALLEL_JOBS )); do
      wait "${pids[0]}"
      pids=("${pids[@]:1}")
    done
  done

  for pid in "${pids[@]}"; do
    wait "$pid"
  done

  : > "$TMP_VCF"
  first_batch=1
  for batch_file in "${batch_outputs[@]}"; do
    if [[ "$first_batch" -eq 1 ]]; then
      cat "$batch_file" >> "$TMP_VCF"
      first_batch=0
    else
      awk '!/^#/' "$batch_file" >> "$TMP_VCF"
    fi
  done

  for batch_log in "${batch_logs[@]}"; do
    cat "$batch_log" >&2
  done
elif [[ "$REGION_MODE" -eq 1 ]]; then
  "$TABIX_BIN" -h "$VCF" "${SELECTED_CHROMS[@]}" | run_python_to_output "-" "$TMP_VCF"
else
  run_python_to_output "$VCF" "$TMP_VCF"
fi

"$BGZIP_BIN" -c "$TMP_VCF" > "$TMP_OUT"
mv "$TMP_OUT" "$OUT"
rm -f "$TMP_VCF"

printf 'Finished: %s\n' "$OUT"
