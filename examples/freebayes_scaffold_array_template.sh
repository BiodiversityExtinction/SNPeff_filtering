#!/usr/bin/env bash
# FreeBayes SNP calling in parallel over one BED file per scaffold/region.
# Configure REF, BAMLIST, BEDDIR, and OUTDIR in the environment before sbatch.

#SBATCH --job-name=freebayes
#SBATCH --cpus-per-task=1
#SBATCH --mem=6G
#SBATCH --time=100:00:00
#SBATCH --output=logs/freebayes_%A_%a.out
#SBATCH --error=logs/freebayes_%A_%a.err

set -euo pipefail

ref="${REF:?Set REF to the indexed reference FASTA}"
bamlist="${BAMLIST:?Set BAMLIST to a text file containing one BAM path per line}"
beddir="${BEDDIR:?Set BEDDIR to a directory containing one .bed file per scaffold/region}"
outdir="${OUTDIR:?Set OUTDIR for the per-scaffold VCF output}"
freebayes_bin="${FREEBAYES:-freebayes}"

for file in "$ref" "$bamlist"; do
  [[ -s "$file" ]] || { echo "ERROR: Missing or empty file: $file" >&2; exit 2; }
done
[[ -d "$beddir" ]] || { echo "ERROR: BED directory not found: $beddir" >&2; exit 2; }

mkdir -p "$outdir" logs

mapfile -t beds < <(find "$beddir" -maxdepth 1 -type f -name '*.bed' | sort)
if [[ ${#beds[@]} -eq 0 ]]; then
  echo "ERROR: No .bed files found in $beddir" >&2
  exit 2
fi

task_id="${SLURM_ARRAY_TASK_ID:?Submit this script as a Slurm array}"
if (( task_id < 0 || task_id >= ${#beds[@]} )); then
  echo "ERROR: Array index $task_id outside 0-$((${#beds[@]} - 1))" >&2
  exit 2
fi

bed="${beds[$task_id]}"
name="$(basename "$bed" .bed)"
out="$outdir/${name}.vcf.gz"
tmp="${out}.tmp"
trap 'rm -f "$tmp"' EXIT

echo "BED chunk: $bed"
echo "BAMs: $(grep -cve '^\s*$' "$bamlist")"
echo "Output: $out"

"$freebayes_bin" \
  -f "$ref" \
  -t "$bed" \
  --min-mapping-quality 20 \
  --min-base-quality 20 \
  --min-alternate-fraction 0.25 \
  -L "$bamlist" \
| bcftools norm -f "$ref" -m-any -Ou \
| bcftools view -i 'QUAL>30' -v snps -m2 -M2 -Ou \
| bcftools sort -Oz -o "$tmp"

mv "$tmp" "$out"
tabix -f -p vcf "$out"
bcftools index -n "$out"
trap - EXIT
echo "Finished: $out"
