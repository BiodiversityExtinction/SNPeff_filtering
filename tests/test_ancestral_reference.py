"""End-to-end checks for SNP polarisation and the matching snpEff FASTA."""

import gzip
import os
from pathlib import Path
import shutil
import subprocess
import tempfile
import unittest


ROOT = Path(__file__).resolve().parents[1]


@unittest.skipUnless(all(shutil.which(tool) for tool in ("bgzip", "tabix", "bcftools")), "htslib tools required")
class AncestralReferenceTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)
        self.fasta = self.root / "original.fa"
        self.fasta.write_text(">chr1\nACAAA\n>chr2\nAAAAA\n")
        self.groups = self.root / "groups.tsv"
        self.groups.write_text("out1 outgroup\nout2 outgroup\nin ingroup\n")
        self.input_vcf = self.root / "input.vcf.gz"
        raw_vcf = self.root / "input.vcf"
        raw_vcf.write_text(
            "##fileformat=VCFv4.2\n"
            "##contig=<ID=chr1,length=5>\n"
            "##contig=<ID=chr2,length=5>\n"
            "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n"
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tout1\tout2\tin\n"
            "chr1\t1\t.\tA\tG\t.\tPASS\t.\tGT\t1/1\t1/1\t0/0\n"
            "chr1\t2\t.\tC\tT\t.\tPASS\t.\tGT\t0/0\t0/0\t1/1\n"
            "chr2\t1\t.\tA\tG\t.\tPASS\t.\tGT\t1/1\t1/1\t0/1\n"
        )
        with self.input_vcf.open("wb") as out:
            subprocess.run(["bgzip", "-c", str(raw_vcf)], stdout=out, check=True)
        subprocess.run(["tabix", "-p", "vcf", str(self.input_vcf)], check=True)
        self.polarised = self.root / "input.polarised.vcf.gz"
        self.patches = self.root / "input.polarised.ancestral_patches.vcf.gz"
        self.ancestral = self.root / "ancestral.fa"

    def run_script(self, name, *args, ok=True, env=None):
        result = subprocess.run(["bash", str(ROOT / name), *map(str, args)],
                                capture_output=True, text=True, env=env)
        if ok:
            self.assertEqual(result.returncode, 0, result.stderr)
        else:
            self.assertNotEqual(result.returncode, 0)
        return result

    def polarise(self, batch=False):
        args = ["--vcf", self.input_vcf, "--groups", self.groups, "--out", self.polarised]
        env = os.environ.copy()
        if batch:
            args += ["-c", "1"]
            env["SLURM_CPUS_PER_TASK"] = "2"
        self.run_script("polarise_snpeff_vcf.sh", *args, env=env)

    def test_serial_and_batched_patches_make_matching_fasta(self):
        for batch in (False, True):
            with self.subTest(batch=batch):
                self.polarise(batch)
                with gzip.open(self.patches, "rt") as handle:
                    records = [line.split("\t") for line in handle if not line.startswith("#")]
                self.assertEqual([(r[0], r[1], r[3], r[4]) for r in records],
                                 [("chr1", "1", "A", "G"), ("chr2", "1", "A", "G")])
                self.run_script("prepare_snpeff_ancestral_fasta.sh", "--fasta", self.fasta,
                                "--patches", self.patches, "--polarised-vcf", self.polarised,
                                "--out", self.ancestral)
                self.assertEqual(self.ancestral.read_text(), ">chr1\nGCAAA\n>chr2\nGAAAA\n")
                if not batch:
                    self.polarised.unlink()
                    self.patches.unlink()
                    (self.root / "input.polarised.ancestral_patches.vcf.gz.csi").unlink()
                    self.ancestral.unlink()

    def test_mismatch_is_rejected_before_annotation(self):
        self.polarise()
        self.run_script("prepare_snpeff_ancestral_fasta.sh", "--fasta", self.fasta,
                        "--patches", self.patches, "--polarised-vcf", self.polarised,
                        "--out", self.ancestral)
        db_dir = self.root / "data" / "test"
        db_dir.mkdir(parents=True)
        (db_dir / "snpEffectPredictor.bin").touch()
        link = db_dir / "sequences.fa"
        link.symlink_to(self.fasta)
        config = self.root / "snpEff.config"
        config.write_text("test.genome : Test\n")
        jar = self.root / "snpEff.jar"
        jar.touch()
        java = self.root / "fake_java"
        java.write_text("#!/usr/bin/env bash\nfor arg in \"$@\"; do last=\"$arg\"; done\ngzip -cd \"$last\"\n")
        java.chmod(0o755)
        args = ["--vcf", self.polarised, "--db-id", "test", "--jar", jar,
                "--config", config, "--java", java]
        result = self.run_script("annotate_snpeff_vcf.sh", *args, ok=False)
        self.assertIn("VCF REF does not match", result.stderr)
        link.unlink()
        link.symlink_to(self.ancestral)
        self.run_script("annotate_snpeff_vcf.sh", *args)
        self.assertTrue((self.root / "input.polarised.anno.vcf.gz").is_file())

    def test_reference_is_unchanged_when_no_sites_flip(self):
        raw_vcf = self.root / "input.vcf"
        rows = raw_vcf.read_text().splitlines()
        rows = [row for row in rows if row.startswith("#") or row.startswith("chr1\t2\t")]
        raw_vcf.write_text("\n".join(rows) + "\n")
        with self.input_vcf.open("wb") as out:
            subprocess.run(["bgzip", "-c", str(raw_vcf)], stdout=out, check=True)
        self.polarise()
        self.run_script("prepare_snpeff_ancestral_fasta.sh", "--fasta", self.fasta,
                        "--patches", self.patches, "--polarised-vcf", self.polarised,
                        "--out", self.ancestral)
        self.assertEqual(self.ancestral.read_text(), self.fasta.read_text())

    def test_builder_links_relative_fasta_to_absolute_path(self):
        gff = self.root / "genes.gff3"
        gff.write_text("##gff-version 3\nchr1\ttest\tgene\t1\t3\t.\t+\t.\tID=g1\n")
        template = self.root / "template.config"
        template.write_text("data.dir = ./data\n")
        result = subprocess.run(
            ["bash", str(ROOT / "build_snpeff_db.sh"), "--fasta", "original.fa",
             "--gff", "genes.gff3", "--outdir", "db", "--db-id", "relative",
             "--config-template", "template.config"],
            cwd=self.root, capture_output=True, text=True,
        )
        self.assertEqual(result.returncode, 0, result.stderr)
        link = self.root / "db" / "data" / "relative" / "sequences.fa"
        self.assertTrue(link.is_symlink())
        self.assertEqual(link.resolve(), self.fasta.resolve())

    def test_force_regeneration_removes_stale_patch_index(self):
        self.polarise()
        subprocess.run(["bcftools", "index", "--csi", str(self.patches)], check=True)
        index = Path(str(self.patches) + ".csi")
        self.assertTrue(index.is_file())
        self.run_script("polarise_snpeff_vcf.sh", "--vcf", self.input_vcf,
                        "--groups", self.groups, "--out", self.polarised, "--force")
        self.assertFalse(index.exists())

    def test_batching_creates_csi_index_when_input_is_unindexed(self):
        Path(str(self.input_vcf) + ".tbi").unlink()
        self.polarise(batch=True)
        self.assertTrue(Path(str(self.input_vcf) + ".csi").is_file())

    @unittest.skipUnless(os.environ.get("SNPEFF_TEST_JAR") and os.environ.get("SNPEFF_TEST_CONFIG"),
                         "Set SNPEFF_TEST_JAR and SNPEFF_TEST_CONFIG for a real snpEff build")
    def test_real_snpeff_calls_flipped_coding_site_missense(self):
        self.polarise()
        self.run_script("prepare_snpeff_ancestral_fasta.sh", "--fasta", self.fasta,
                        "--patches", self.patches, "--polarised-vcf", self.polarised,
                        "--out", self.ancestral)
        gff = self.root / "genes.gff3"
        gff.write_text(
            "##gff-version 3\n"
            "chr1\ttest\tgene\t1\t3\t.\t+\t.\tID=gene1;Name=gene1\n"
            "chr1\ttest\tmRNA\t1\t3\t.\t+\t.\tID=tx1;Parent=gene1\n"
            "chr1\ttest\texon\t1\t3\t.\t+\t.\tID=ex1;Parent=tx1\n"
            "chr1\ttest\tCDS\t1\t3\t.\t+\t0\tID=cds1;Parent=tx1\n"
        )
        snpeff_root = self.root / "snpEff"
        jar = Path(os.environ["SNPEFF_TEST_JAR"])
        self.run_script("build_snpeff_db.sh", "--fasta", self.ancestral,
                        "--gff", gff, "--outdir", snpeff_root,
                        "--db-id", "test_ancestral", "--config-template",
                        os.environ["SNPEFF_TEST_CONFIG"], "--jar", jar, "--build")
        self.run_script("annotate_snpeff_vcf.sh", "--vcf", self.polarised,
                        "--db-id", "test_ancestral", "--jar", jar,
                        "--config", snpeff_root / "snpEff.config")
        with gzip.open(self.root / "input.polarised.anno.vcf.gz", "rt") as handle:
            first_record = next(line for line in handle if not line.startswith("#"))
        self.assertIn("missense_variant", first_record)
        self.assertNotIn("WARNING_REF_DOES_NOT_MATCH_GENOME", first_record)

    @unittest.skipUnless(os.environ.get("SNPEFF_TEST_JAR") and os.environ.get("SNPEFF_TEST_CONFIG"),
                         "Set SNPEFF_TEST_JAR and SNPEFF_TEST_CONFIG for a real snpEff build")
    def test_runner_completes_load_summary(self):
        gff = self.root / "genes.gff3"
        gff.write_text(
            "##gff-version 3\n"
            "chr1\ttest\tgene\t1\t3\t.\t+\t.\tID=gene1;Name=gene1\n"
            "chr1\ttest\tmRNA\t1\t3\t.\t+\t.\tID=tx1;Parent=gene1\n"
            "chr1\ttest\texon\t1\t3\t.\t+\t.\tID=ex1;Parent=tx1\n"
            "chr1\ttest\tCDS\t1\t3\t.\t+\t0\tID=cds1;Parent=tx1\n"
        )
        self.run_script("run_ancestral_snpeff_load.sh",
                        "--vcf", self.input_vcf, "--groups", self.groups,
                        "--fasta", self.fasta, "--gff", gff,
                        "--outdir", self.root / "snpEff", "--db-id", "test_ancestral",
                        "--jar", os.environ["SNPEFF_TEST_JAR"],
                        "--config-template", os.environ["SNPEFF_TEST_CONFIG"],
                        "--workdir", self.root / "results", "--prefix", "case",
                        "--jackknife-bins", "0")
        summary = self.root / "results" / "case.polarised.ancestral_ref.anno.load_summary.tsv"
        self.assertTrue(summary.is_file())
        self.assertIn("nonsynonymous_realised", summary.read_text())


if __name__ == "__main__":
    unittest.main()
