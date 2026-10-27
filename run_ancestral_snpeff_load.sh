#!/usr/bin/env bash

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

usage() {
  cat <<'EOF'
Usage:
  run_ancestral_snpeff_load.sh --vcf FILE --groups FILE --fasta FILE --gff FILE \
    --outdir DIR --db-id NEW_ID --jar FILE [options]

Required:
  --vcf FILE           Unannotated input VCF/VCF.GZ
  --groups FILE        Ingroup/outgroup sample table
  --fasta FILE         Original FASTA used for calling, or its contig subset
  --gff FILE           Original GFF/GFF3 annotation
  --outdir DIR         snpEff root directory (contains snpEff.config and data/)
  --db-id NEW_ID       New database ID for this ancestral FASTA
  --jar FILE           snpEff.jar

Optional:
  --workdir DIR        VCF, ancestral FASTA and table outputs (default: VCF directory)
  --prefix NAME        Output file prefix (default: input VCF basename)
  --config-template FILE
                       Config template if OUTDIR/snpEff.config does not yet exist
  -c, --chromosomes-per-batch N
                       Internal polarisation batch size (default: 0, serial)
  --jackknife-bins N   Gene bins for load uncertainty (default: 100)
  --allow-missing-outgroup
                       Use the same relaxed outgroup rule in both steps
  -h, --help           Show this help

This runs: polarise VCF -> make ancestral FASTA -> build NEW snpEff database
-> annotate against that database -> summarise load. It will not overwrite
existing outputs. Use a fresh prefix and database ID for a changed input or
ancestral-state definition. No Slurm directives are embedded in this script.
EOF
}

die() { printf 'Error: %s\n' "$1" >&2; exit 1; }

VCF=""
GROUPS_FILE=""
FASTA=""
GFF=""
OUTDIR=""
DB_ID=""
JAR=""
WORKDIR=""
PREFIX=""
CONFIG_TEMPLATE=""
CHROMS_PER_BATCH=0
JACKKNIFE_BINS=100
ALLOW_MISSING=0

while [[ $# -gt 0 ]]; do
  case "$1" in
    --vcf) VCF="${2:-}"; shift 2 ;;
    --groups) GROUPS_FILE="${2:-}"; shift 2 ;;
    --fasta) FASTA="${2:-}"; shift 2 ;;
    --gff) GFF="${2:-}"; shift 2 ;;
    --outdir) OUTDIR="${2:-}"; shift 2 ;;
    --db-id) DB_ID="${2:-}"; shift 2 ;;
    --jar) JAR="${2:-}"; shift 2 ;;
    --workdir) WORKDIR="${2:-}"; shift 2 ;;
    --prefix) PREFIX="${2:-}"; shift 2 ;;
    --config-template) CONFIG_TEMPLATE="${2:-}"; shift 2 ;;
    -c|--chromosomes-per-batch) CHROMS_PER_BATCH="${2:-}"; shift 2 ;;
    --jackknife-bins) JACKKNIFE_BINS="${2:-}"; shift 2 ;;
    --allow-missing-outgroup) ALLOW_MISSING=1; shift ;;
    -h|--help) usage; exit 0 ;;
    *) die "Unknown argument: $1" ;;
  esac
done

[[ -n "$VCF" ]] || die "--vcf is required"
[[ -n "$GROUPS_FILE" ]] || die "--groups is required"
[[ -n "$FASTA" ]] || die "--fasta is required"
[[ -n "$GFF" ]] || die "--gff is required"
[[ -n "$OUTDIR" ]] || die "--outdir is required"
[[ -n "$DB_ID" ]] || die "--db-id is required"
[[ -n "$JAR" ]] || die "--jar is required"
for file in "$VCF" "$GROUPS_FILE" "$FASTA" "$GFF" "$JAR"; do
  [[ -f "$file" ]] || die "File not found: $file"
done
[[ "$CHROMS_PER_BATCH" =~ ^[0-9]+$ ]] || die "-c must be a non-negative integer"
[[ "$JACKKNIFE_BINS" =~ ^[0-9]+$ ]] || die "--jackknife-bins must be a non-negative integer"
[[ "$DB_ID" =~ ^[A-Za-z0-9_-]+$ ]] || die "--db-id may contain letters, digits, underscores and hyphens only"

if [[ -z "$WORKDIR" ]]; then
  WORKDIR="$(dirname "$VCF")"
fi
if [[ -z "$PREFIX" ]]; then
  PREFIX="$(basename "$VCF")"
  PREFIX="${PREFIX%.gz}"
  PREFIX="${PREFIX%.vcf}"
fi
[[ "$PREFIX" =~ ^[A-Za-z0-9_.-]+$ ]] || die "--prefix must be a plain filename (no slashes or spaces)"

mkdir -p "$WORKDIR" "$OUTDIR"
POLARISED="${WORKDIR}/${PREFIX}.polarised.vcf.gz"
PATCHES="${WORKDIR}/${PREFIX}.polarised.ancestral_patches.vcf.gz"
ANCESTRAL_FASTA="${WORKDIR}/${PREFIX}.ancestral.fa"
ANNOTATED="${WORKDIR}/${PREFIX}.polarised.ancestral_ref.anno.vcf.gz"
SUMMARY="${WORKDIR}/${PREFIX}.polarised.ancestral_ref.anno.load_summary.tsv"

polarise_args=(--vcf "$VCF" --groups "$GROUPS_FILE" --out "$POLARISED" --ancestral-patches "$PATCHES" -c "$CHROMS_PER_BATCH")
summary_args=(--vcf "$ANNOTATED" --groups "$GROUPS_FILE" --out "$SUMMARY" --jackknife-bins "$JACKKNIFE_BINS")
if (( ALLOW_MISSING )); then
  polarise_args+=(--allow-missing-outgroup)
  summary_args+=(--allow-missing-outgroup)
fi

printf 'Step 1/5: polarise VCF and record ancestral patches\n'
bash "$SCRIPT_DIR/polarise_snpeff_vcf.sh" "${polarise_args[@]}"

printf 'Step 2/5: make and validate ancestral FASTA\n'
bash "$SCRIPT_DIR/prepare_snpeff_ancestral_fasta.sh" \
  --fasta "$FASTA" --patches "$PATCHES" --polarised-vcf "$POLARISED" --out "$ANCESTRAL_FASTA"

printf 'Step 3/5: build a new snpEff database\n'
build_args=(--fasta "$ANCESTRAL_FASTA" --gff "$GFF" --outdir "$OUTDIR" --db-id "$DB_ID" --jar "$JAR" --build)
if [[ -n "$CONFIG_TEMPLATE" ]]; then
  build_args+=(--config-template "$CONFIG_TEMPLATE")
fi
bash "$SCRIPT_DIR/build_snpeff_db.sh" "${build_args[@]}"

printf 'Step 4/5: annotate only against the matching database\n'
bash "$SCRIPT_DIR/annotate_snpeff_vcf.sh" \
  --vcf "$POLARISED" --db-id "$DB_ID" --jar "$JAR" \
  --config "$OUTDIR/snpEff.config" --out "$ANNOTATED"

printf 'Step 5/5: summarise load\n'
bash "$SCRIPT_DIR/summarise_snpeff_load.sh" "${summary_args[@]}"

printf 'Finished: %s\n' "$SUMMARY"
