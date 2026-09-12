# MVA Hackathon 2026, Track 1: Phenotype-Anchored, ClinVar-First Prioritization (samwell)

Pipeline and methods for the Rare Disease, Real Kid MVA Hackathon 2026 (Track 1, variant
prediction).

**Finding:** a biallelic *BUB1B* compound heterozygote, p.Leu737* (in ClinVar as pathogenic for
MVA syndrome 1) together with p.Asn1002Lys (a kinase-domain missense absent from gnomAD). A robust depth and
allele-fraction scan of the same WGS also finds a low-level allelic-imbalance signal (mosaic gains
of chr20 and chr22 in roughly 5% of blood cells), the direction MVA predicts. Full details are in
[report/methods_report.md](report/methods_report.md).

**Contents:** `pipeline.py` (the full reproducible analysis), `predictions/` (the submission CSV),
and `report/` (the methods report). No patient data is stored in this repository. Access to the
challenge dataset is gated by the organizers (WCG IRB #20252010). The pipeline uses public
references only (ClinVar, gnomAD streamed by region, UCSC RefSeq), and patient genotypes never
leave the analysis machine.

License: CC BY 4.0, per the hackathon rules.
