# Track 1 Methods Report: Phenotype-Anchored, ClinVar-First Variant Prioritization

**Team:** samwell · **Case:** MVA Hackathon 2026 proband (WGS_EX2312012) · **Build:** GRCh38

## Summary of findings

**Primary candidate (rank 1): a biallelic *BUB1B* compound heterozygote.**

| Allele | GRCh38 | HGVS (NM_001211.6) | Evidence |
|---|---|---|---|
| 1 | chr15:40209701 T>G | c.2210T>G, p.Leu737* (nonsense) | ClinVar Pathogenic/Likely pathogenic for Mosaic variegated aneuploidy syndrome 1 (VCV000533901, multiple submitters, no conflicts); rs759242053; gnomAD AF about 3e-5; PASS, GT 0/1, AD 21,25 |
| 2 | chr15:40220612 T>G | c.3006T>G, p.Asn1002Lys (missense, exon 23) | Absent from gnomAD v4.1 genomes; sits in the BUBR1 kinase domain; SIFT 0.01 (deleterious); PolyPhen-2 0.997 (probably damaging); PASS, GT 0/1, AD 15,13 |

Biallelic loss of *BUB1B* function causes MVA syndrome 1 (MIM 257300), and the proband's phenotype
matches it point for point: rhabdomyosarcoma (the malignancy most specific to *BUB1B*-related MVA),
growth restriction and short stature, nephrocalcinosis, failure to thrive, prematurity, and a
parental history of recurrent miscarriage. The genotype we observe, one truncating allele together
with a rarer missense allele on the other copy, is the configuration most often reported in
*BUB1B*-related MVA (Matsuura et al. described the pattern of one null allele plus one
reduced-function allele).

**The same file also shows the aneuploidy itself.** A per-chromosome depth and allele-fraction scan
of the WGS finds chromosome-level gains consistent with mosaic trisomies in blood: chr21 mean depth
55.8x against a 43.1x autosomal baseline (a 29% excess), with chr22 up 21%, chr20 up 16%, and chr16
up 16%. Each of these chromosomes also shows a wider spread of heterozygous allele fractions
(population standard deviation 0.116 to 0.120, against a baseline near 0.092) and a slightly
depressed median alternate fraction. GC bias does not explain the pattern, because the GC-rich
chr19 sits at baseline depth. So the variant-level diagnosis and the chromosomal phenotype are both
visible in one dataset, independently of each other.

**Alternative second alleles (ranks 2 and 3), kept as hedges:**
- chr15:40192892 C>T, an intronic *BUB1B* variant of uncertain significance (rs185599777,
  VCV001676498, gnomAD AF 0.21%, intron 8, about 3.6 kb from the exon 9 acceptor). It would fit a
  low-expression allele model, but its population frequency (including one gnomAD homozygote) and
  its depth in the intron make it a weaker candidate than p.Asn1002Lys.
- chr15:40216470 A>G, deep in intron 20 and absent from gnomAD.

**Secondary and incidental findings (ranks 5 to 7, flagged as secondary):** heterozygous
carrier-grade ClinVar pathogenic variants in *LZTR1* (chr22:20996720 C>G), *PRSS1*
(chr7:142750561 C>T), and *HK1* (chr10:69315762 A>G). None has a second allele anywhere in the
genome and none matches the presenting phenotype. We list them because the rules invite secondary
findings.

## Method

1. **Phenotype anchoring.** The HPO cluster (HP:0002859 rhabdomyosarcoma, HP:0000121
   nephrocalcinosis, HP:0004322, HP:0001508, HP:0001622, HP:0001518, HP:0200067 parental recurrent
   miscarriage) is the published clinical signature of the mosaic variegated aneuploidy syndromes,
   so the four known MVA genes (*BUB1B*, *CEP57*, *TRIP13*, *CENATAC*) formed the prior panel, with
   rhabdomyosarcoma specifically favoring *BUB1B*.
2. **Targeted scan.** We extracted every PASS variant in the four gene bodies from the provided
   GATK VCF with tabix and intersected them with a locally downloaded copy of the ClinVar GRCh38
   VCF. This step alone surfaced the pathogenic *BUB1B* nonsense and the intronic VUS.
3. **Second-allele search.** We annotated every heterozygous *BUB1B* gene-body variant with a
   population frequency by streaming the relevant slices of gnomAD v4.1 genomes over remote tabix,
   so only region coordinates ever left the analysis machine, never patient genotypes. We then
   mapped each variant onto the NM_001211.6 exon model (UCSC RefSeq) and computed codon-level
   consequences locally. The kinase-domain missense p.Asn1002Lys emerged as the only coding
   heterozygote absent from gnomAD.
4. **Genome-wide exclusion.** We intersected every PASS variant in the genome with all 346,571
   ClinVar pathogenic or likely pathogenic records. There were 7 hits in total. Six are
   heterozygous carrier findings for recessive or unrelated conditions, which leaves *BUB1B* as the
   only gene with both phenotype concordance and a plausible biallelic genotype. No competing
   diagnosis survives this pass.
5. **Aneuploidy scan.** A single streaming pass over all PASS variants collected per-chromosome
   mean depth and heterozygous allele-fraction distributions (DP between 15 and 80), flagging
   chromosomes that deviate in depth ratio, median allele fraction, or allele-fraction spread. This
   recovered the mosaic chromosomal gains described above.
6. **Ranking.** The ClinVar-anchored compound heterozygote received epcr 0.90, with the alternates
   and incidentals ranked beneath it at decreasing epcr.

## Limitations

- **Phase is inferred, not proven.** No parental samples are provided, and the two variants sit
  10.9 kb apart, which is beyond short-read phasing range. We infer the trans configuration from
  allele frequencies (a cis haplotype carrying both a 3e-5 allele and a gnomAD-absent allele is far
  less likely than trans) and from the clinical diagnosis itself.
- p.Asn1002Lys is unreported as far as we can tell (absent from ClinVar and gnomAD), so its
  pathogenicity rests on in-silico consensus, its domain location, its rarity, and the need for a
  second allele under a recessive model with a confirmed phenotype.
- The aneuploidy analysis reflects the blood compartment at the time of sampling, nothing more.

## Reproducibility

All code in this repository (`pipeline.py`) reruns the full analysis from the provided VCF plus
public references (ClinVar, gnomAD, UCSC RefSeq). No patient data is redistributed here. The run
takes about 10 minutes on a laptop, with no GPU and no external service ever receiving patient
genotypes.
