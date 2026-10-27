#!/usr/bin/env bash

set -euo pipefail

usage() {
  cat <<'EOF'
Usage:
  prepare_snpeff_ancestral_fasta.sh --fasta FILE --patches FILE --polarised-vcf FILE --out FILE [options]

Required arguments:
  --fasta FILE           Original reference FASTA used for variant calling
  --patches FILE         Original-orientation patch VCF from polarise_snpeff_vcf.sh
  --polarised-vcf FILE   Companion polarised VCF; every REF will be checked
  --out FILE             Output ancestral FASTA for a new snpEff database

Optional arguments:
  --bcftools BIN         bcftools binary (default: from PATH)
  --force                Replace an existing output FASTA
  --skip-if-exists       Skip if the output FASTA already exists
  -h, --help             Show this help

The patch VCF must have original-genome REF and ancestral ALT. Do not apply
the polarised VCF as a patch: its REF is already ancestral. This script checks
patch REF against the original FASTA, applies the patch with bcftools consensus,
then checks every polarised VCF REF against the new FASTA. No sites are held
in memory by this shell script.
EOF
}

die() { printf 'Error: %s\n' "$1" >&2; exit 1; }

FASTA=""
PATCHES=""
POLARISED_VCF=""
OUT=""
BCFTOOLS_BIN=""
FORCE=0
SKIP_IF_EXISTS=0

while [[ $# -gt 0 ]]; do
  case "$1" in
    --fasta) FASTA="${2:-}"; shift 2 ;;
    --patches) PATCHES="${2:-}"; shift 2 ;;
    --polarised-vcf) POLARISED_VCF="${2:-}"; shift 2 ;;
    --out) OUT="${2:-}"; shift 2 ;;
    --bcftools) BCFTOOLS_BIN="${2:-}"; shift 2 ;;
    --force) FORCE=1; shift ;;
    --skip-if-exists) SKIP_IF_EXISTS=1; shift ;;
    -h|--help) usage; exit 0 ;;
    *) die "Unknown argument: $1" ;;
  esac
done

[[ -n "$FASTA" && -n "$PATCHES" && -n "$POLARISED_VCF" && -n "$OUT" ]] || die "--fasta, --patches, --polarised-vcf and --out are required"
for file in "$FASTA" "$PATCHES" "$POLARISED_VCF"; do
  [[ -f "$file" ]] || die "File not found: $file"
done
[[ "$OUT" != "$FASTA" ]] || die "Output FASTA must differ from original FASTA"
[[ "$OUT" != *.gz ]] || die "Output FASTA must be uncompressed for snpEff"
if [[ -z "$BCFTOOLS_BIN" ]]; then
  BCFTOOLS_BIN="$(command -v bcftools || true)"
fi
[[ -n "$BCFTOOLS_BIN" && -x "$BCFTOOLS_BIN" ]] || die "bcftools not found"

if [[ -e "$OUT" ]]; then
  if (( SKIP_IF_EXISTS && ! FORCE )); then
    printf 'Checking existing ancestral FASTA before skipping...\n'
    "$BCFTOOLS_BIN" norm --check-ref e --fasta-ref "$FASTA" -Ou "$PATCHES" > /dev/null
    "$BCFTOOLS_BIN" norm --check-ref e --fasta-ref "$OUT" -Ou "$POLARISED_VCF" > /dev/null
    printf 'Validated existing ancestral FASTA, skipping: %s\n' "$OUT"
    exit 0
  fi
  (( FORCE )) || die "Output exists: $OUT (use --force or --skip-if-exists)"
fi

mkdir -p "$(dirname "$OUT")"
TMP_OUT="${OUT}.tmp.$$"
trap 'rm -f "$TMP_OUT" "${TMP_OUT}.fai"' EXIT

printf 'Checking patch REF against original FASTA...\n'
"$BCFTOOLS_BIN" norm --check-ref e --fasta-ref "$FASTA" -Ou "$PATCHES" > /dev/null

if [[ ! -f "${PATCHES}.csi" && ! -f "${PATCHES}.tbi" ]]; then
  printf 'Indexing patch VCF for bcftools consensus...\n'
  "$BCFTOOLS_BIN" index --csi "$PATCHES"
fi

printf 'Applying ancestral patches to original FASTA...\n'
"$BCFTOOLS_BIN" consensus --fasta-ref "$FASTA" "$PATCHES" > "$TMP_OUT"
[[ -s "$TMP_OUT" ]] || die "Consensus produced an empty FASTA"

printf 'Checking every polarised VCF REF against ancestral FASTA...\n'
"$BCFTOOLS_BIN" norm --check-ref e --fasta-ref "$TMP_OUT" -Ou "$POLARISED_VCF" > /dev/null

mv "$TMP_OUT" "$OUT"
if [[ -f "${TMP_OUT}.fai" ]]; then
  mv "${TMP_OUT}.fai" "${OUT}.fai"
else
  rm -f "${OUT}.fai"
fi
printf 'Finished ancestral FASTA: %s\n' "$OUT"
