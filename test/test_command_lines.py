import unittest
from types import SimpleNamespace
from IGenotyper.command_lines.alignments import Align
from IGenotyper.command_lines.snps import Snps


class CommandCapture:
    def __init__(self):
        self.calls = []

    def __call__(self, command, output):
        self.calls.append((command, output))


class CommandLineTests(unittest.TestCase):
    def setUp(self):
        self.files = SimpleNamespace(ref="reference.fasta")
        self.cpu = SimpleNamespace(threads=8)

    def test_minimap2_hifi_command(self):
        align = Align(self.files, self.cpu, "sample")
        capture = CommandCapture()
        align.run_command = capture
        align.map_reads_with_minimap2("reads.fasta", "mapped.bam", self.files.ref)
        self.assertEqual(capture.calls[0][1], "mapped.bam.bai")
        self.assertIn("minimap2 -t 8 -a -x map-hifi", capture.calls[0][0])
        self.assertIn("| samtools sort -@ 8 -o mapped.bam -", capture.calls[0][0])

    def test_minimap2_streams_bam_as_fasta(self):
        align = Align(self.files, self.cpu, "sample")
        capture = CommandCapture()
        align.run_command = capture
        align.map_reads_with_minimap2(
            "subreads.bam", "mapped.bam", self.files.ref, "map-pb"
        )
        self.assertIn("samtools fasta -@ 8 subreads.bam |", capture.calls[0][0])
        self.assertIn("-x map-pb reference.fasta -", capture.calls[0][0])

    def test_assembly_uses_asm20_only(self):
        self.files.assembly_fasta = "contigs.fasta"
        self.files.assembly_to_ref = "contigs.bam"
        align = Align(self.files, self.cpu, "sample")
        capture = CommandCapture()
        align.run_command = capture
        align.map_assembly()
        self.assertIn("-x asm20", capture.calls[0][0])
        self.assertEqual(capture.calls[0][0].count("-x asm20"), 1)

    def test_whatshap_runs_without_conda_environment_switch(self):
        snps = Snps(self.files, self.cpu, "sample 1")
        capture = CommandCapture()
        snps.run_command = capture
        snps.snp_candidates("reads.bam", "candidates.vcf")
        command = capture.calls[0][0]
        self.assertIn("whatshap find_snv_candidates", command)
        self.assertIn("--pacbio", command)
        self.assertNotIn("conda", command)


if __name__ == "__main__":
    unittest.main()
