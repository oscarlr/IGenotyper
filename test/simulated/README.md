# Simulated integration dataset

This deterministic diploid dataset contains a 2 kb `igh` reference and twelve
perfect HiFi-like reads: six from each haplotype. Haplotype 2 carries four
heterozygous SNVs recorded in `truth.vcf`.

Regenerate and test it with:

```bash
python test/simulated/generate.py
test/simulated/run.sh
```

The integration test maps reads with minimap2, sorts and indexes them with
SAMtools, finds and genotypes candidates with WhatsHap, phases the calls, and
checks all four positions against the truth VCF.
