#!/usr/bin/env bash
#SBATCH --job-name=snpeff_annotate
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=48:00:00
#SBATCH --output=slurm-%j.out

set -euo pipefail

usage() {
  cat <<'EOF'
Usage:
  annotate_snpeff_vcf.sh --vcf FILE --db-id ID [options]

Required arguments:
  --vcf FILE           Input VCF/VCF.GZ file to annotate
  --db-id ID           snpEff database ID (for example: brown_bear)

Optional arguments:
  --out FILE           Output annotated VCF.GZ
                       Default: input with .anno.vcf.gz suffix
  --jar FILE           Path to snpEff.jar
                       Default: none; must be provided unless you edit the script
  --config FILE        snpEff config file
                       Default: none; must be provided unless you edit the script
  --java BIN           Java binary to use
                       Default: java from PATH
  --bgzip BIN          bgzip binary to use
                       Default: bgzip from PATH
  --force              Overwrite output if it already exists
  --skip-if-exists     Exit cleanly if output already exists
  -h, --help           Show this help

What it does:
  1. Checks that the input VCF, snpEff jar, and config file exist
  2. Verifies that the requested database is present in the config
  3. Runs snpEff annotation
  4. Compresses the output with bgzip

Example:
  sbatch annotate_snpeff_vcf.sh \
    --vcf /path/sample.vcf.gz \
    --db-id brown_bear \
    --jar /path/snpEff.jar \
    --config /path/snpEff.config
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
DB_ID=""
OUT=""
SNPEFF_JAR=""
CONFIG=""
JAVA_BIN=""
BGZIP_BIN=""
FORCE=0
SKIP_IF_EXISTS=0

while [[ $# -gt 0 ]]; do
  case "$1" in
    --vcf)
      VCF="${2:-}"
      shift 2
      ;;
    --db-id)
      DB_ID="${2:-}"
      shift 2
      ;;
    --out)
      OUT="${2:-}"
      shift 2
      ;;
    --jar)
      SNPEFF_JAR="${2:-}"
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
    --bgzip)
      BGZIP_BIN="${2:-}"
      shift 2
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
[[ -n "$DB_ID" ]] || die "--db-id is required"

require_file "$VCF"
[[ -n "$SNPEFF_JAR" ]] || die "--jar is required unless you set a default in the script"
[[ -n "$CONFIG" ]] || die "--config is required unless you set a default in the script"

require_file "$CONFIG"
require_file "$SNPEFF_JAR"

if [[ -z "$JAVA_BIN" ]]; then
  JAVA_BIN="$(command -v java || true)"
fi
[[ -n "$JAVA_BIN" ]] || die "Java not found in PATH and --java was not provided"
require_file "$JAVA_BIN"

if [[ -z "$BGZIP_BIN" ]]; then
  BGZIP_BIN="$(command -v bgzip || true)"
fi
[[ -n "$BGZIP_BIN" ]] || die "bgzip not found in PATH and --bgzip was not provided"
require_file "$BGZIP_BIN"

if ! grep -q "^${DB_ID}\.genome[[:space:]]*:" "$CONFIG"; then
  die "Database ID '${DB_ID}' not found in config: $CONFIG"
fi

if [[ -z "$OUT" ]]; then
  case "$VCF" in
    *.vcf.gz)
      OUT="${VCF%.vcf.gz}.anno.vcf.gz"
      ;;
    *.vcf)
      OUT="${VCF%.vcf}.anno.vcf.gz"
      ;;
    *)
      OUT="${VCF}.anno.vcf.gz"
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

TMP_OUT="${OUT}.tmp"
rm -f "$TMP_OUT"

printf 'Annotating VCF with snpEff\n'
printf '  input  = %s\n' "$VCF"
printf '  output = %s\n' "$OUT"
printf '  db     = %s\n' "$DB_ID"
printf '  config = %s\n' "$CONFIG"

"$JAVA_BIN" -jar "$SNPEFF_JAR" -v "$DB_ID" -c "$CONFIG" "$VCF" | "$BGZIP_BIN" -c > "$TMP_OUT"
mv "$TMP_OUT" "$OUT"

printf 'Finished: %s\n' "$OUT"
