# SNPeff Load Workflow Notes

This folder contains five core scripts that I use together to build a `snpEff` database, polarise a VCF to an outgroup-defined ancestral state, annotate that VCF, summarise realised and masked genetic load, and then plot the results in a consistent way across projects. An optional FreeBayes example shows one way to generate the input VCF.

This is not meant to be polished software documentation. It is a working note for future analyses, so the main aim is to explain clearly what each script does, what goes in, what comes out, and what assumptions are being made.

## Example: Make An Input VCF

[`examples/freebayes_scaffold_array_template.sh`](./examples/freebayes_scaffold_array_template.sh)
shows one way to call variants before starting the `snpEff` workflow. It is an
example of an upstream calling pipeline, not a required part of the load scripts.
It uses FreeBayes on a BAM list, with one Slurm array task per BED file, then
normalises and filters the calls with BCFtools. The output is one indexed,
compressed VCF per BED region. The example assumes Slurm, FreeBayes, BCFtools
and tabix are available.

Prepare a FASTA reference with an index, indexed BAM files, a text file listing
one BAM path per line, and non-overlapping BED files for the regions to call.
Each BED filename becomes its output VCF basename, so filenames must be unique.
From the directory where you will run `sbatch`, create `logs/` **before**
submission: Slurm opens its output files before the script can create directories.
Set `REF`, `BAMLIST`, `BEDDIR` and `OUTDIR` to your own paths and export them so
the array tasks can read them. If `FREEBAYES` is unset, the script uses
`freebayes` from `PATH`.

For example, from the repository root:

```bash
export REF=/path/reference.fa
export BAMLIST=/path/samples.bamlist
export BEDDIR=/path/scaffold_beds
export OUTDIR=/path/scaffold_vcfs
mkdir -p logs "$OUTDIR"
mapfile -t beds < <(find "$BEDDIR" -maxdepth 1 -type f -name '*.bed' | sort)
# Confirm that the BED list is non-empty before submitting.
sbatch --array="0-$((${#beds[@]} - 1))" examples/freebayes_scaffold_array_template.sh
```

The example requires mapping quality at least 20, base quality at least 20 and
an alternate allele fraction of at least 0.25. It keeps biallelic SNPs with
`QUAL > 30`. These settings are **example choices**, not universal thresholds;
record the values actually used for each study. SNP-only output excludes indels
and therefore frameshift indels from downstream LoF counts. This script does
not remove sites with missing genotypes. If that is required, do it explicitly
after combining the VCFs, and record that filtering step in the methods.

After the array finishes, combine the VCFs in genomic coordinate order, check
that regions do not overlap or duplicate sites, and index the combined VCF.
Use that combined VCF as input to polarisation (if needed) and subsequent
annotation. The exact combination and missing-data filtering steps depend on
how the BED regions and samples were prepared; they are not done by the example
script. Use the same reference FASTA for calling and for the `snpEff` database.

## Scripts In This Folder

- `build_snpeff_db.sh`
- `polarise_snpeff_vcf.sh`
- `annotate_snpeff_vcf.sh`
- `summarise_snpeff_load.sh`
- `plot_snpeff_load.R`

## Overall Logic

The workflow is:

1. Prepare a `snpEff` reference database from a FASTA and GFF.
2. Polarise a VCF so the outgroup-defined ancestral allele becomes `REF`.
3. Annotate the polarised VCF with `snpEff`.
4. Summarise load statistics from the annotated VCF.
5. Plot per-individual load values and pairwise population `Rxy`.

The purpose of keeping these steps separate is:

- the reference database usually only needs to be built once per reference genome
- the polarisation step may be repeated with different outgroup definitions
- annotation may be repeated for multiple polarised VCFs against the same reference
- load summaries may be rerun multiple times with different outgroup definitions, population files, or uncertainty settings
- plots may be rerun with different population files without touching the earlier steps

### When to use the polarisation step

This step is mainly for datasets mapped to an ingroup reference.

- If the VCF was mapped to an ingroup reference, I should usually polarise first.
- If the VCF was mapped to an outgroup reference, I can usually skip the polarisation step and go straight to `snpEff` annotation.

The reason is simple:

- with an ingroup reference, the reference allele is not guaranteed to be ancestral
- with an outgroup reference, the reference allele is often a much better approximation to the ancestral state

So in practice:

- ingroup-mapped dataset: build -> polarise -> annotate -> summarise -> plot
- outgroup-mapped dataset: build -> annotate -> summarise -> plot

## Software Needed

These scripts assume the following are available:

- `bash`
- `java`
- `python3`
- `bgzip`
- `snpEff.jar`
- a `snpEff.config` file

They were written for HPC use and can be submitted with `sbatch`, but they also run directly with `bash`.

## Step 1: Build The SNPeff Database

Script:
- [`build_snpeff_db.sh`](./build_snpeff_db.sh)

### Why this step exists

`snpEff` needs a database directory with a very specific structure. In particular, it expects:

- a FASTA called `sequences.fa`
- a GFF called `genes.gff`
- both placed under `data/<db-id>/`
- an entry for the database in `snpEff.config`

This script prepares those files in a reproducible way.

### What the script does

Given a reference FASTA and a GFF/GFF3 file, it:

1. creates `OUTDIR/data/<db-id>/`
2. creates `sequences.fa` as a symlink to the input FASTA
3. filters the GFF so that only sequence IDs found in the FASTA are retained
4. writes that filtered annotation as `genes.gff`
5. adds `<db-id>.genome : <genome-name>` to `snpEff.config` if the entry is missing
6. optionally runs `snpEff build`

The filtering step matters because `snpEff` is happier when the annotation and FASTA sequence names are perfectly aligned.

### Important input arguments

- `--fasta`: reference FASTA
- `--gff`: annotation GFF or GFF3, gzipped or plain text
- `--outdir`: root `snpEff` directory
- `--db-id`: short database ID used by `snpEff`

Useful optional arguments:

- `--genome-name`: display name in the config
- `--config-template`: template `snpEff.config` to copy if one does not exist yet
- `--jar`: path to `snpEff.jar`
- `--build`: actually run the build after preparing the files
- `--force`: overwrite prepared files
- `--skip-build-if-done`: skip rebuilding if the database already exists

### Example

```bash
bash build_snpeff_db.sh \
  --fasta /path/reference.fna \
  --gff /path/reference.gff.gz \
  --outdir /path/SNPeff \
  --db-id brown_bear \
  --genome-name Brown_bear \
  --config-template /path/snpEff.config \
  --jar /path/snpEff.jar \
  --build
```

### What should exist after this step

At minimum:

```text
SNPeff/
  snpEff.config
  data/
    brown_bear/
      sequences.fa
      genes.gff
```

If `--build` was used successfully, the database directory should also contain binary `snpEff` build files such as:

```text
snpEffectPredictor.bin
sequence.bin
sequence.<contig>.bin
```

## Step 2: Polarise A VCF

Script:
- [`polarise_snpeff_vcf.sh`](./polarise_snpeff_vcf.sh)

### Why this step exists

For load analyses, what I really want is a VCF where:

- `REF` is the ancestral allele
- `ALT` is the derived allele

If I use the raw reference-based VCF, that assumption only holds when the reference genome already carries the ancestral state. That can bias load estimates, especially if the reference genome comes from one of the focal populations.

This script fixes that for a load-only workflow by using outgroups to decide which allele is ancestral and then rewriting the VCF before running `snpEff`.

### What the script does

It:

1. reads a VCF or VCF.GZ
2. keeps biallelic SNPs only
3. checks the outgroup genotypes at each site
4. keeps a site only if all called outgroups agree on the same homozygous allele
5. if outgroups are fixed for `ALT`, swaps `REF` and `ALT`
6. flips all sample genotypes to match the new allele order
7. writes a simplified VCF with `FORMAT=GT` only

The simplification step is deliberate. Once a site has been flipped, many allele-dependent fields such as `AD`, `PL`, `AC`, and `AF` are no longer safe to trust unless they are recalculated. Since I only need this VCF for `snpEff` and the downstream load script, it is cleaner to keep genotypes only.

### Optional parallel batching inside one job

For large VCFs, the polarisation step can be sped up without using a Slurm array. The script can:

1. detect chromosome names automatically from a bgzipped VCF
2. create a tabix index automatically if it is missing
3. split chromosomes into internal batches
4. run those batches in parallel within a single job
5. merge the batch outputs back into one final polarised VCF

This is controlled with:

- `-c N` or `--chromosomes-per-batch N`

The number of simultaneous batch workers is taken from `SLURM_CPUS_PER_TASK` if available, so I can submit one job with multiple CPUs and let the script manage the parallel chromosome batches internally.

For example:

```bash
sbatch --cpus-per-task=10 polarise_snpeff_vcf.sh \
  --vcf /path/input.vcf.gz \
  --groups /path/Pops.txt \
  -c 10
```

In that example:

- `--cpus-per-task=10` means up to 10 chromosome batches can run at once
- `-c 10` means each batch contains 10 chromosomes

These are different settings. One controls parallelism, the other controls batch size.

### Important input arguments

- `--vcf`: unannotated input VCF or VCF.GZ
- `--groups`: two-column sample group file with `ingroup` and `outgroup`

Useful optional arguments:

- `--out`: custom output file name
- `-c N`: process chromosomes in internal batches of `N` at a time
- `--allow-missing-outgroup`: allow missing outgroup calls, as long as all called outgroups still agree on the same homozygous allele
- `--force`: overwrite existing output
- `--skip-if-exists`: do nothing if the output already exists

### Important limitations

- this is meant for load analyses, not as a general-purpose rewritten VCF
- it keeps biallelic SNPs only
- indels and multiallelic sites are skipped automatically
- it will not recover indel-based LoF classes such as `frameshift_variant`
- after polarisation, the VCF should be re-annotated with `snpEff`

So if I use this route, my LoF category becomes effectively SNP-based LoF unless I later extend the same logic to indels.

### Example

```bash
sbatch polarise_snpeff_vcf.sh \
  --vcf /path/input.vcf.gz \
  --groups /path/Pops.txt
```

### Default output

If the input is:

```text
input.vcf.gz
```

the default output will be:

```text
input.polarised.vcf.gz
```

### What to check after this step

- the output VCF should exist and be bgzipped
- `FORMAT` should now be `GT`
- some sites may now have a different `REF` and `ALT` than in the original VCF
- the run summary should report how many sites were kept and how many were flipped

## Step 3: Annotate A VCF

Script:
- [`annotate_snpeff_vcf.sh`](./annotate_snpeff_vcf.sh)

### Why this step exists

Once the database is built, the next step is to annotate variants with predicted effects. This script runs `snpEff` and writes a compressed annotated VCF.

The most important new information added by `snpEff` is the `ANN=` field in the VCF `INFO` column. The load script uses those annotations to decide whether a site is LoF, nonsynonymous, synonymous, modifier, or something else.

### What the script does

It:

1. checks that the input VCF exists
2. checks that `snpEff.jar` and the config file exist
3. checks that the requested database ID is present in the config
4. runs `snpEff`
5. compresses the result with `bgzip`

### Important input arguments

- `--vcf`: input VCF or VCF.GZ
- `--db-id`: the database name built in step 1

Useful optional arguments:

- `--out`: custom output file name
- `--jar`: path to `snpEff.jar`
- `--config`: explicit config file path
- `--force`: overwrite existing output
- `--skip-if-exists`: do nothing if the output already exists

### Example

```bash
sbatch annotate_snpeff_vcf.sh \
  --vcf /path/input.polarised.vcf.gz \
  --db-id brown_bear \
  --jar /path/snpEff.jar \
  --config /path/snpEff.config
```

### Default output

If the input is:

```text
input.polarised.vcf.gz
```

the default output will be:

```text
input.polarised.anno.vcf.gz
```

### What to check after this step

- the output VCF should exist and be bgzipped
- the header should contain an `ANN` definition
- records should contain `ANN=` annotations in the `INFO` field

## Step 4: Summarise Realised And Masked Load

Script:
- [`summarise_snpeff_load.sh`](./summarise_snpeff_load.sh)

### Why this step exists

The annotated VCF still contains site-by-site information. This script turns that into per-individual load summaries.

The main idea is:

- start from the polarised and annotated VCF
- treat `REF` as ancestral and `ALT` as derived
- use the group file to decide which samples are ingroup and which are outgroup
- keep only sites where the outgroups still support the ancestral state
- count how often each ingroup individual carries derived alleles in different annotation classes

### Outgroup assumption

The script uses a strict ancestral-state rule:

- every outgroup individual must be `0/0`

If any outgroup carries an ALT allele, the site is excluded.

This is intentionally conservative. It reduces the number of usable sites, but it makes the polarisation clearer and more consistent.

There is also an optional relaxed mode:

- `--allow-missing-outgroup`

In that case, missing outgroup genotypes are allowed, but all called outgroups still need to be `0/0`.

### Group file format

The script needs a two-column tab-delimited file telling it which samples are ingroup and which are outgroup.

Example:

```tsv
sample	group
Polar_bear_UAmap	outgroup
American_black_bear_UAmap	outgroup
SW_01_UAmap	ingroup
SW_02_UAmap	ingroup
```

Accepted group labels:

- `ingroup`
- `outgroup`

### How sites are classified

The script follows the same strict annotation rules each time.

#### LoF

A site is counted as loss of function if `snpEff` annotates it as one of:

- `stop_gained`
- `frameshift_variant`
- `splice_acceptor_variant`
- `splice_donor_variant`

#### Nonsynonymous

A site is counted as nonsynonymous if it is annotated as:

- `missense_variant`

#### Synonymous

A site is counted as synonymous if it is annotated as one of:

- `synonymous_variant`
- `start_retained_variant`
- `stop_retained_variant`

#### Modifier

If a site is not in one of the categories above but has `MODIFIER` impact, it is counted as modifier.

#### Other

Anything not captured above goes into `other`.

### Realised and masked load

The script separates homozygous and heterozygous derived states:

- `realised` = homozygous derived
- `masked` = heterozygous derived

So, for example:

- a homozygous derived LoF genotype contributes to `lof_realised`
- a heterozygous derived LoF genotype contributes to `lof_masked`

This keeps the biologically useful distinction between load that is exposed in homozygotes and load that is currently masked in heterozygotes.

### Why synonymous_total is used

The script also calculates:

```text
synonymous_total = synonymous_realised + 0.5 * synonymous_masked
```

This is used as a denominator for normalised load measures.

The logic is:

- homozygous derived synonymous sites contribute `1`
- heterozygous derived synonymous sites contribute `0.5`

This gives a simple dosage-like neutral denominator while still keeping realised and masked load separate in the numerator.

### Output columns

The output table contains one row per ingroup sample.

The main raw-count columns are:

- `lof_realised`
- `lof_masked`
- `nonsynonymous_realised`
- `nonsynonymous_masked`
- `synonymous_realised`
- `synonymous_masked`
- `modifier_realised`
- `modifier_masked`

The main normalised columns are:

- `synonymous_total`
- `lof_realised_over_synonymous_total`
- `lof_masked_over_synonymous_total`
- `nonsynonymous_realised_over_synonymous_total`
- `nonsynonymous_masked_over_synonymous_total`

Floating-point values are written to 6 decimal places.

### Separate genotype-specific normalisation table

Each run also writes `input.anno.load_summary.genotype_normalised.tsv`. The
original summary keeps its shared `synonymous_total` denominator. The additional
table uses the synonymous count for the same genotype state as its numerator:

| Column | Calculation |
| --- | --- |
| `lof_realised_over_synonymous_realised` | Homozygous derived LoF / homozygous derived synonymous |
| `lof_masked_over_synonymous_masked` | Heterozygous LoF / heterozygous synonymous |
| `nonsynonymous_realised_over_synonymous_realised` | Homozygous derived missense / homozygous derived synonymous |
| `nonsynonymous_masked_over_synonymous_masked` | Heterozygous missense / heterozygous synonymous |

The table includes `sample` and all six underlying counts, so the denominators
and the ratios can be checked. Homozygous synonymous here means **homozygous
derived**, not homozygous reference. There is no 0.5 weighting in these ratios:
both numerator and denominator count the same genotype state. Variant
classification, outgroup filtering and sample inclusion are identical for the
two tables. Numbers use British spelling and six decimal places.

These additional ratios describe predicted damaging variation relative to
synonymous variation within each genotype state. They help examine whether a
population difference tracks background heterozygosity or homozygosity, but they
do not measure the absolute number of damaging genotypes. For example, 100
heterozygous LoF variants / 10,000 heterozygous synonymous variants and 200 /
20,000 give the same ratio despite different burdens. Use both tables and the
raw counts when interpreting differences; neither normalisation removes all
effects of demography, callability or annotation quality.

With `--jackknife-bins N`, the new table includes SE and 95% CI columns for all
four ratios. The script uses the **same bins** for both tables and recalculates
both numerator and denominator after each bin is omitted. The new intervals
cannot be recovered from the old summary's SEs/CIs alone: rerun the annotated
VCF to calculate them.

A zero denominator gives `NA`. If any leave-one-bin-out denominator is zero, the
full-data ratio is still reported when defined, but its SE and CI are `NA`.
Uncertainty is also `NA` when jackknifing is disabled or fewer than two bins are
available. `_jackknife_bins` records the number of bins used/attempted, so read it
together with the SE/CI fields. Zero observed numerator counts can produce zero
SE and a zero-width interval; this does not establish that the underlying burden
is certainly zero. Normal-approximation CI bounds are not clipped at zero.

### Jackknife uncertainty

If `--jackknife-bins N` is provided, the script also estimates uncertainty.

The logic is:

1. genes are tracked from the `ANN` field
2. genes are ordered by first appearance in the VCF
3. genes are assigned round-robin to approximately equal-sized bins (gene index modulo bin count)
4. one bin is left out at a time
5. the normalised load metrics are recalculated
6. jackknife SE and 95% CI are calculated from those replicate values

Variants assigned to the same gene stay together. If an annotation has no gene
identifier/name, the script uses chromosome, position and ALT as a fallback unit.
These bins are not contiguous genomic windows and do not necessarily preserve
linkage between neighbouring genes. The intervals describe variation across the
chosen bins; they do not include uncertainty from genotype errors or annotation.

That is useful because:

- sites inside the same gene are not independent
- the unit of interest here is closer to genes than to random SNPs across the whole genome

The 95% confidence interval is:

```text
estimate +/- 1.96 * SE
```

So each normalised metric can have:

- `_jackknife_bins`
- `_jackknife_se`
- `_jackknife_ci_lower`
- `_jackknife_ci_upper`

### Important input arguments

- `--vcf`: annotated VCF from step 3
- `--groups`: ingroup/outgroup table

Useful optional arguments:

- `--out`: output file name
- `--genotype-out`: name for the separate genotype-specific table; defaults to
  `--out` with `.tsv` replaced by `.genotype_normalised.tsv`
- `--genotype-only`: write only the genotype-specific table, preserving any existing summary
- `--jackknife-bins N`: add uncertainty using `N` gene bins
- `--allow-missing-outgroup`: allow missing outgroup calls
- `--force`: overwrite selected outputs
- `--skip-if-exists`: do nothing only if all selected outputs already exist;
  a partially existing pair requires `--genotype-only` or `--force`

### Example

```bash
sbatch summarise_snpeff_load.sh \
  --vcf /path/input.polarised.anno.vcf.gz \
  --groups /path/Pops.txt \
  --jackknife-bins 100
```

### Default output

If the input is:

```text
input.polarised.anno.vcf.gz
```

the default output will be:

```text
input.polarised.anno.load_summary.tsv
input.polarised.anno.load_summary.genotype_normalised.tsv
```

To add the new table to a completed analysis with 100 jackknife bins:

```bash
bash summarise_snpeff_load.sh \
  --vcf /path/input.polarised.anno.vcf.gz \
  --groups /path/Pops.txt \
  --jackknife-bins 100 \
  --genotype-only
```

If the original summary used a custom `--out`, supply that same `--out` to derive
the matching companion filename, or specify `--genotype-out` explicitly. The
plotting script currently reads the original shared-denominator summary; the new
table is a separate analysis output.

## Step 5: Plot Individual Load And Population Rxy

Script:
- [`plot_snpeff_load.R`](./plot_snpeff_load.R)

### Why this step exists

After making the load summary table, I usually want two things:

- per-individual plots of realised and masked load
- pairwise `Rxy` comparisons between populations

This script takes the final load summary table and a sample-to-population file, filters to the samples I want to include, and writes both figures and tables.

### Input files

The script needs:

1. a load summary table from `summarise_snpeff_load.sh`
2. a two-column sample/population file

The plotting population file is separate from the ingroup/outgroup file used during the load calculation. That way I can keep the load definition fixed but regroup samples however I want for plotting.

### Plot population file format

Example:

```tsv
DRR250459_UAmap	Japan
DRR320311_UAmap	Japan
NE_01_UAmap	Mainland
NE_02_UAmap	Mainland
SW_01_UAmap	Mainland
SW_02_UAmap	Mainland
```

Only samples present in this file are kept for plotting. If a sample is in the load summary table but not in this file, it is ignored.

### What the script produces

#### Individual load plots

The script makes:

- one LoF plot
- one nonsynonymous plot

For each individual:

- realised load is plotted on the left y-axis
- masked load is plotted on the right y-axis
- points are coloured by population
- confidence intervals are drawn if present in the load table

The four plotted values of interest are:

- `lof_realised_over_synonymous_total`
- `lof_masked_over_synonymous_total`
- `nonsynonymous_realised_over_synonymous_total`
- `nonsynonymous_masked_over_synonymous_total`

#### Rxy output

The script also calculates pairwise population `Rxy` from the same four metrics.

If there are exactly two populations:

- it calculates all pairwise between-population sample ratios
- it writes those values to `*.rxy_pairwise.tsv`
- it plots them as boxplots
- it overlays the population-mean `Rxy` as a point
- it draws a horizontal line at `Rxy = 1`

If there are more than two populations:

- it writes the pairwise population summary table
- it makes a heatmap instead

### Output files

Typical outputs are:

```text
prefix.lof_dual_axis.pdf
prefix.lof_dual_axis.png
prefix.nonsynonymous_dual_axis.pdf
prefix.nonsynonymous_dual_axis.png
prefix.rxy.tsv
prefix.population_means.tsv
prefix.samples_used.tsv
```

If there are two populations:

```text
prefix.rxy_boxplot.pdf
prefix.rxy_boxplot.png
prefix.rxy_pairwise.tsv
```

If there are more than two populations:

```text
prefix.rxy_heatmap.pdf
prefix.rxy_heatmap.png
```

### Example

```bash
Rscript plot_snpeff_load.R \
  --load /path/input.anno.load_summary.tsv \
  --groups /path/Plot_pops.txt \
  --outdir /path/load_plots \
  --prefix modern_bears
```

## End-To-End Example

### 1. Build the database once for a reference

```bash
bash build_snpeff_db.sh \
  --fasta /path/reference.fna \
  --gff /path/reference.gff.gz \
  --outdir /path/SNPeff \
  --db-id brown_bear \
  --config-template /path/snpEff.config \
  --jar /path/snpEff.jar \
  --build
```

### 2. Polarise a VCF if needed

```bash
sbatch polarise_snpeff_vcf.sh \
  --vcf /path/modern_bears.vcf.gz \
  --groups /path/Pops.txt
```

### 3. Annotate the polarised VCF against that database

```bash
sbatch annotate_snpeff_vcf.sh \
  --vcf /path/modern_bears.polarised.vcf.gz \
  --db-id brown_bear \
  --jar /path/snpEff.jar \
  --config /path/SNPeff/snpEff.config
```

### 4. Make a population file

Example:

```tsv
sample	group
Polar_bear_UAmap	outgroup
American_black_bear_UAmap	outgroup
SW_01_UAmap	ingroup
SW_02_UAmap	ingroup
SW_03_UAmap	ingroup
```

### 5. Summarise realised and masked load

```bash
sbatch summarise_snpeff_load.sh \
  --vcf /path/modern_bears.polarised.anno.vcf.gz \
  --groups /path/Pops.txt \
  --jackknife-bins 100
```

### 6. Plot individual load and Rxy

```bash
Rscript plot_snpeff_load.R \
  --load /path/modern_bears.polarised.anno.load_summary.tsv \
  --groups /path/Plot_pops.txt \
  --outdir /path/load_plots \
  --prefix modern_bears
```

## Practical Notes For Future Use

- The database build usually only needs to be done once per reference genome.
- If a script says the output already exists, either use `--force` or provide a new `--out` file name.
- If I want to preserve earlier versions of an output table, it is usually cleaner to give a new `--out` path rather than overwriting.
- The load summary table only includes ingroup samples. Outgroup samples are used for filtering but are not reported as rows.
- If there are no synonymous sites for a sample after filtering, ratio columns become `NA`.
- Negative lower confidence intervals can happen mathematically when the estimate is small and the SE is relatively large.
- The scripts default to `java`, `python3`, and `bgzip` from `PATH`, so I usually only need to specify paths when working on a new system.

## Short Reminder

If I forget the order:

1. build database
2. prepare `Pops.txt`
3. polarise VCF
4. annotate VCF
5. summarise load
6. prepare `Plot_pops.txt`
7. make plots
