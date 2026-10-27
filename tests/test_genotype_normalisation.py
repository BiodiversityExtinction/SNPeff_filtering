"""Small VCF integration tests for the two load normalisations."""

import csv
import math
from pathlib import Path
import subprocess
import tempfile
import unittest


SCRIPT = Path(__file__).resolve().parents[1] / "summarise_snpeff_load.sh"


class GenotypeNormalisationTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)
        self.groups = self.root / "groups.tsv"
        self.groups.write_text("out outgroup\nbear ingroup\nempty ingroup\n")
        self.vcf = self.root / "input.vcf"
        rows = ["##fileformat=VCFv4.2", "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tout\tbear\tempty"]
        self.records = []
        for gene in range(4):
            self.records.extend([(gene, "synonymous_variant", "1/1"), (gene, "synonymous_variant", "0/1")])
        self.records.extend([
            (0, "stop_gained", "1/1"), (1, "stop_gained", "0/1"),
            (2, "missense_variant", "1/1"), (3, "missense_variant", "0/1"),
        ])
        self.header = rows
        self.write_vcf(self.records)
        self.out = self.root / "summary.tsv"
        self.genotype_out = self.root / "summary.genotype_normalised.tsv"

    def write_vcf(self, records):
        rows = list(self.header)
        for pos, (gene, effect, gt) in enumerate(records, 1):
            ann = f"ANN=G|{effect}|HIGH|gene{gene}|id{gene}"
            rows.append(f"chr1\t{pos}\t.\tA\tG\t.\tPASS\t{ann}\tGT\t0/0\t{gt}\t0/0")
        # An outgroup ALT carrier must never contribute to either table.
        rows.append("chr1\t999\t.\tA\tG\t.\tPASS\tANN=G|stop_gained|HIGH|excluded|excluded\tGT\t0/1\t1/1\t0/0")
        self.vcf.write_text("\n".join(rows) + "\n")

    def run_script(self, *args, ok=True):
        result = subprocess.run([
            "bash", str(SCRIPT), "--vcf", str(self.vcf), "--groups", str(self.groups),
            "--out", str(self.out), *args,
        ], capture_output=True, text=True)
        if ok:
            self.assertEqual(result.returncode, 0, result.stderr)
        else:
            self.assertNotEqual(result.returncode, 0)
        return result

    def read_table(self, path):
        with path.open() as handle:
            return list(csv.DictReader(handle, delimiter="\t"))

    def test_ratios_and_independently_calculated_jackknife(self):
        self.run_script("--jackknife-bins", "4")
        old = self.read_table(self.out)[0]
        row, empty = self.read_table(self.genotype_out)
        self.assertEqual(old["synonymous_total"], "6.000000")
        self.assertEqual(old["lof_realised_over_synonymous_total"], "0.166667")
        self.assertNotIn("lof_realised_over_synonymous_realised", old)
        # Each genotype denominator loses one synonymous count per omitted gene.
        replicates = [0, 1 / 3, 1 / 3, 1 / 3]
        mean = sum(replicates) / 4
        se = math.sqrt(3 / 4 * sum((x - mean) ** 2 for x in replicates))
        for category in ("lof", "nonsynonymous"):
            for state in ("realised", "masked"):
                metric = f"{category}_{state}_over_synonymous_{state}"
                self.assertEqual(row[metric], "0.250000")
                self.assertEqual(row[metric + "_jackknife_bins"], "4")
                self.assertEqual(row[metric + "_jackknife_se"], f"{se:.6f}")
                self.assertEqual(row[metric + "_jackknife_ci_lower"], f"{0.25 - 1.96 * se:.6f}")
                self.assertEqual(row[metric + "_jackknife_ci_upper"], f"{0.25 + 1.96 * se:.6f}")
                self.assertEqual(empty[metric], "NA")
                self.assertEqual(empty[metric + "_jackknife_se"], "NA")

    def test_undefined_replicate_does_not_discard_valid_point_estimate(self):
        records = [(g, e, gt) for g, e, gt in self.records if e != "synonymous_variant" or gt != "0/1" or g == 0]
        self.write_vcf(records)
        self.run_script("--jackknife-bins", "4")
        row = self.read_table(self.genotype_out)[0]
        self.assertEqual(row["lof_masked_over_synonymous_masked"], "1.000000")
        self.assertEqual(row["lof_masked_over_synonymous_masked_jackknife_se"], "NA")
        self.assertEqual(row["lof_realised_over_synonymous_realised_jackknife_se"], "0.250000")

    def test_disabled_and_single_bin_uncertainty(self):
        for bins in (0, 1):
            self.run_script("--jackknife-bins", str(bins), "--force")
            row = self.read_table(self.genotype_out)[0]
            self.assertEqual(row["lof_realised_over_synonymous_realised"], "0.250000")
            self.assertEqual(row["lof_realised_over_synonymous_realised_jackknife_se"], "NA")

    def test_companion_only_skip_and_collision_handling(self):
        self.out.write_text("preserve the existing summary\n")
        self.run_script("--skip-if-exists", ok=False)
        self.run_script("--genotype-only", "--jackknife-bins", "4")
        self.assertEqual(self.out.read_text(), "preserve the existing summary\n")
        before = self.genotype_out.read_bytes()
        self.run_script("--skip-if-exists")
        self.assertEqual(self.genotype_out.read_bytes(), before)
        self.run_script("--genotype-only", ok=False)
        self.run_script("--genotype-out", str(self.out), "--force", ok=False)
        self.run_script("--genotype-out", str(self.vcf), "--force", ok=False)


if __name__ == "__main__":
    unittest.main()
