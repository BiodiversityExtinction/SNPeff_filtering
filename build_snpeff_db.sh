#!/usr/bin/env bash

set -euo pipefail

usage() {
  cat <<'EOF'
Usage:
  build_snpeff_db.sh --fasta FILE --gff FILE --outdir DIR --db-id ID [options]

Required arguments:
  --fasta FILE         Reference FASTA file
  --gff FILE           Annotation GFF/GFF3 file (.gz accepted)
  --outdir DIR         snpEff root output directory
  --db-id ID           Database ID used under data/ID and in the build command

Optional arguments:
  --genome-name NAME   Display name written to snpEff.config (default: same as --db-id)
  --config-template FILE
                       Template snpEff.config to copy if outdir/snpEff.config does not exist
  --config FILE        Explicit config file path to update/use (default: OUTDIR/snpEff.config)
  --java BIN           Java binary for the build step (default: java from PATH)
  --jar FILE           snpEff jar for the build step
  --build              Run snpEff build after preparing the files
  --skip-build-if-done Skip snpEff build if a database marker already exists
  --force              Overwrite existing genes.gff and re-create symlink
  -h, --help           Show this help

What it does:
  1. Creates OUTDIR/data/DB_ID
  2. Creates sequences.fa as a symlink to the input FASTA
  3. Filters the GFF so only sequence IDs found in the FASTA are kept
  4. Writes genes.gff in the snpEff database folder
  5. Adds "DB_ID.genome : GENOME_NAME" to snpEff.config if missing
  6. Optionally runs snpEff build
  7. Reuses existing files when possible unless --force is given

Example:
  build_snpeff_db.sh \
    --fasta /path/genome.fna \
    --gff /path/genes.gff.gz \
    --outdir /path/SNPeff \
    --db-id brown_bear \
    --genome-name Brown_bear \
    --config-template /path/snpEff.config \
    --java /path/java \
    --jar /path/snpEff.jar \
    --build
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

FASTA=""
GFF=""
OUTDIR=""
DB_ID=""
GENOME_NAME=""
CONFIG_TEMPLATE=""
CONFIG=""
JAVA_BIN=""
SNPEFF_JAR=""
RUN_BUILD=0
SKIP_BUILD_IF_DONE=0
FORCE=0

while [[ $# -gt 0 ]]; do
  case "$1" in
    --fasta)
      FASTA="${2:-}"
      shift 2
      ;;
    --gff)
      GFF="${2:-}"
      shift 2
      ;;
    --outdir)
      OUTDIR="${2:-}"
      shift 2
      ;;
    --db-id)
      DB_ID="${2:-}"
      shift 2
      ;;
    --genome-name)
      GENOME_NAME="${2:-}"
      shift 2
      ;;
    --config-template)
      CONFIG_TEMPLATE="${2:-}"
      shift 2
      ;;
    --config)
      CONFIG="${2:-}"
      shift 2
      ;;
    --java)
      JAVA_BIN="${2:-}"
      shift 2
      ;;
    --jar)
      SNPEFF_JAR="${2:-}"
      shift 2
      ;;
    --build)
      RUN_BUILD=1
      shift
      ;;
    --skip-build-if-done)
      SKIP_BUILD_IF_DONE=1
      shift
      ;;
    --force)
      FORCE=1
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

[[ -n "$FASTA" ]] || die "--fasta is required"
[[ -n "$GFF" ]] || die "--gff is required"
[[ -n "$OUTDIR" ]] || die "--outdir is required"
[[ -n "$DB_ID" ]] || die "--db-id is required"

require_file "$FASTA"
require_file "$GFF"

if [[ -z "$GENOME_NAME" ]]; then
  GENOME_NAME="$DB_ID"
fi

if [[ -z "$CONFIG" ]]; then
  CONFIG="${OUTDIR}/snpEff.config"
fi

DATA_DIR="${OUTDIR}/data/${DB_ID}"
SEQUENCES_LINK="${DATA_DIR}/sequences.fa"
GENES_GFF="${DATA_DIR}/genes.gff"
SEQIDS_FILE="${DATA_DIR}/seqids.txt"
TMP_GFF="${DATA_DIR}/genes.all.gff"
SNPEFF_BIN_FILE="${DATA_DIR}/snpEffectPredictor.bin"

mkdir -p "$DATA_DIR"

if [[ ! -f "$CONFIG" ]]; then
  [[ -n "$CONFIG_TEMPLATE" ]] || die "Config not found at $CONFIG and --config-template was not provided"
  require_file "$CONFIG_TEMPLATE"
  cp "$CONFIG_TEMPLATE" "$CONFIG"
fi

if grep -q "^${DB_ID}\.genome[[:space:]]*:" "$CONFIG"; then
  :
else
  printf '%s.genome : %s\n' "$DB_ID" "$GENOME_NAME" >> "$CONFIG"
fi

if [[ -L "$SEQUENCES_LINK" ]] && [[ "$(readlink "$SEQUENCES_LINK")" == "$FASTA" ]] && [[ "$FORCE" -ne 1 ]]; then
  printf 'Reusing existing sequences.fa symlink: %s\n' "$SEQUENCES_LINK"
else
  if [[ -e "$SEQUENCES_LINK" && "$FORCE" -ne 1 ]]; then
    die "${SEQUENCES_LINK} already exists and does not point to ${FASTA}. Use --force to replace it."
  fi
  ln -sfn "$FASTA" "$SEQUENCES_LINK"
fi

if [[ -f "$GENES_GFF" && "$FORCE" -ne 1 ]]; then
  printf 'Reusing existing genes.gff: %s\n' "$GENES_GFF"
else
  awk '/^>/ {sub(/^>/, "", $1); print $1}' "$SEQUENCES_LINK" > "$SEQIDS_FILE"

  case "$GFF" in
    *.gz)
      gzip -cd "$GFF" > "$TMP_GFF"
      ;;
    *)
      cp "$GFF" "$TMP_GFF"
      ;;
  esac

  awk 'BEGIN{FS=OFS="\t"} NR==FNR {keep[$1]=1; next} /^#/ {print; next} ($1 in keep) {print}' \
    "$SEQIDS_FILE" "$TMP_GFF" > "$GENES_GFF"

  rm -f "$SEQIDS_FILE" "$TMP_GFF"
fi

printf 'Prepared snpEff database files in %s\n' "$DATA_DIR"
printf '  sequences.fa -> %s\n' "$FASTA"
printf '  genes.gff    <- %s\n' "$GFF"
printf '  config       =  %s\n' "$CONFIG"

if [[ "$RUN_BUILD" -eq 1 ]]; then
  if [[ -z "$JAVA_BIN" ]]; then
    JAVA_BIN="$(command -v java || true)"
  fi

  [[ -n "$JAVA_BIN" ]] || die "--build was requested but Java was not found in PATH and --java was not provided"
  [[ -n "$SNPEFF_JAR" ]] || die "--build was requested but --jar was not provided"
  require_file "$JAVA_BIN"
  require_file "$SNPEFF_JAR"

  if [[ "$SKIP_BUILD_IF_DONE" -eq 1 && -f "$SNPEFF_BIN_FILE" && "$FORCE" -ne 1 ]]; then
    printf 'Existing snpEff build found, skipping: %s\n' "$SNPEFF_BIN_FILE"
    exit 0
  fi

  (
    cd "$OUTDIR"
    "$JAVA_BIN" -jar "$SNPEFF_JAR" \
      build -gff3 -v "$DB_ID" -c "$CONFIG" \
      -noCheckCds -noCheckProtein -quiet
  )

  printf 'snpEff build completed for %s\n' "$DB_ID"
fi
