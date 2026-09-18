import tempfile
import unittest
import sys
import types
from pathlib import Path

from IGenotyper.common.reference import validate_bed

sys.modules.setdefault("vcf", types.ModuleType("vcf"))
from IGenotyper.common.vcffn import Variant


class FakeGenotype:
    def __getitem__(self, key):
        if key == "GT":
            return "0|1"
        raise KeyError(key)


class FakeRecord:
    REF = "A"
    ALT = ["G"]

    def genotype(self, sample):
        return FakeGenotype()


class ReferenceAndPhasingTests(unittest.TestCase):
    def test_phased_alleles_are_reusable(self):
        variant = Variant()
        variant.from_record_allele_bases(FakeRecord(), "sample")
        self.assertEqual(tuple(variant.allele_bases), ("A", "G"))
        self.assertEqual(tuple(variant.allele_bases), ("A", "G"))

    def test_bed_validation_rejects_unknown_contig(self):
        with tempfile.TemporaryDirectory() as directory:
            bed = Path(directory) / "bad.bed"
            bed.write_text("missing\t0\t1\n")
            with self.assertRaisesRegex(ValueError, "unknown contig"):
                validate_bed(str(bed), {"igh": 10})


if __name__ == "__main__":
    unittest.main()
