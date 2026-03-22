# SNPeff Load Workflow Notes

This folder contains three scripts that I use together to build a `snpEff` database, annotate a VCF, and then summarise realised and masked genetic load in a consistent way across projects.

This is not meant to be polished software documentation. It is a working note for future analyses, so the main aim is to explain clearly what each script does, what goes in, what comes out, and what assumptions are being made.

## Scripts In This Folder

- `build_snpeff_db.sh`
- `annotate_snpeff_vcf.sh`
- `summarize_snpeff_load.sh`

## Overall Logic

The workflow is:

1. Prepare a `snpEff` reference database from a FASTA and GFF.
2. Annotate a VCF with `snpEff`.
3. Polarise variants using an outgroup and summarise load statistics from the annotated VCF.

The purpose of keeping these three steps separate is:

- the reference database usually only needs to be built once per reference genome
- annotation may be repeated for multiple VCFs against the same reference
- load summaries may be rerun multiple times with different outgroup definitions, population files, or uncertainty settings

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

## Step 2: Annotate A VCF

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
- `--snpeff-dir`: directory containing `snpEff.jar`
- `--config`: explicit config file path
- `--force`: overwrite existing output
- `--skip-if-exists`: do nothing if the output already exists

### Example

```bash
sbatch annotate_snpeff_vcf.sh \
  --vcf /path/input.vcf.gz \
  --db-id brown_bear
```

### Default output

If the input is:

```text
input.vcf.gz
```

the default output will be:

```text
input.anno.vcf.gz
```

### What to check after this step

- the output VCF should exist and be bgzipped
- the header should contain an `ANN` definition
- records should contain `ANN=` annotations in the `INFO` field

## Step 3: Summarise Realised And Masked Load

Script:
- [`summarize_snpeff_load.sh`](./summarize_snpeff_load.sh)

### Why this step exists

The annotated VCF still contains site-by-site information. This script turns that into per-individual load summaries.

The main idea is:

- use outgroup individuals to define the ancestral state
- keep only sites where all outgroups are fixed for the reference allele
- treat the ALT allele as derived at those retained sites
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

### Jackknife uncertainty

If `--jackknife-bins N` is provided, the script also estimates uncertainty.

The logic is:

1. genes are tracked from the `ANN` field
2. genes are ordered by first appearance in the VCF
3. those genes are split into approximately equal-sized bins
4. one bin is left out at a time
5. the normalised load metrics are recalculated
6. jackknife SE and 95% CI are calculated from those replicate values

This means the uncertainty is based on groups of genes rather than on individual sites.

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

- `--vcf`: annotated VCF from step 2
- `--groups`: ingroup/outgroup table

Useful optional arguments:

- `--out`: output file name
- `--jackknife-bins N`: add uncertainty using `N` gene bins
- `--allow-missing-outgroup`: allow missing outgroup calls
- `--force`: overwrite existing output
- `--skip-if-exists`: do nothing if the output already exists

### Example

```bash
sbatch summarize_snpeff_load.sh \
  --vcf /path/input.anno.vcf.gz \
  --groups /path/Pops.txt \
  --jackknife-bins 100
```

### Default output

If the input is:

```text
input.anno.vcf.gz
```

the default output will be:

```text
input.anno.load_summary.tsv
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

### 2. Annotate a VCF against that database

```bash
sbatch annotate_snpeff_vcf.sh \
  --vcf /path/modern_bears.vcf.gz \
  --db-id brown_bear
```

### 3. Make a population file

Example:

```tsv
sample	group
Polar_bear_UAmap	outgroup
American_black_bear_UAmap	outgroup
SW_01_UAmap	ingroup
SW_02_UAmap	ingroup
SW_03_UAmap	ingroup
```

### 4. Summarise realised and masked load

```bash
sbatch summarize_snpeff_load.sh \
  --vcf /path/modern_bears.anno.vcf.gz \
  --groups /path/Pops.txt \
  --jackknife-bins 100
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
2. annotate VCF
3. prepare `Pops.txt`
4. summarise load
