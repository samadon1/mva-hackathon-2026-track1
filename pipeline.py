"""MVA Hackathon 2026 - Track 1 pipeline (samwell).

Phenotype-anchored, ClinVar-first variant prioritization for the MVA proband.
Reruns the complete analysis from the provided VCF + public references.
No patient genotypes are sent to any external service: ClinVar is downloaded
locally; gnomAD and UCSC are queried by genomic REGION only.

Usage:
  python pipeline.py --vcf data/WGS_EX2312012_HGWCNDSX7.vcf.gz --out predictions.csv
Deps: pysam (pip install pysam)
"""
import argparse, gzip, json, statistics, urllib.request
import pysam

MVA_GENES = {  # GRCh38 gene-body windows (no-chr contig names, as in the provided VCF)
    "BUB1B":   ("15", 40_161_088, 40_221_169),
    "CEP57":   ("11", 95_786_807, 95_827_919),
    "TRIP13":  ("5",     892_879,    930_000),
    "CENATAC": ("11", 118_880_000, 118_896_000),
}
CLINVAR = "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz"
GNOMAD = "https://storage.googleapis.com/gcp-public-data--gnomad/release/4.1/vcf/genomes/gnomad.genomes.v4.1.sites.chr{c}.vcf.bgz"
UCSC = "https://api.genome.ucsc.edu/getData/track?genome=hg38;track=ncbiRefSeqCurated;chrom=chr{c};start={s};end={e}"


def clinvar_map(cv, chrom, start, end):
    out = {}
    for row in cv.fetch(chrom, start, end):
        f = row.split("\t")
        info = dict(kv.split("=", 1) for kv in f[7].split(";") if "=" in kv)
        out[(f[1], f[3], f[4])] = info
    return out


def gnomad_af(chrom, pos, ref, alt):
    url = GNOMAD.format(c=chrom)
    g = pysam.TabixFile(url, index=url + ".tbi")
    try:
        for row in g.fetch(f"chr{chrom}", pos - 1, pos):
            f = row.split("\t")
            if f[3] == ref and alt in f[4].split(","):
                ai = f[4].split(",").index(alt)
                info = dict(kv.split("=", 1) for kv in f[7].split(";") if "=" in kv)
                return float(info.get("AF", "nan").split(",")[ai])
        return 0.0  # absent
    finally:
        g.close()


def gene_scan(pat, cv):
    print("== MVA panel scan ==")
    for gene, (c, s, e) in MVA_GENES.items():
        cvm = clinvar_map(cv, c, s, e)
        for row in pat.fetch(c, s, e):
            f = row.split("\t")
            if f[6] != "PASS":
                continue
            key = (f[1], f[3], f[4])
            gt = f[9].split(":")[0]
            sig = cvm.get(key, {}).get("CLNSIG", "")
            af = gnomad_af(c, int(f[1]), f[3], f[4]) if gt in ("0/1", "1/1") else None
            if (sig and "enign" not in sig) or (af is not None and af < 0.005):
                print(f"  {gene} {c}:{f[1]} {f[3]}>{f[4]} GT={gt} gnomAD_AF={af} ClinVar={sig or '-'}")


def genomewide_plp(vcf_path):
    print("== genome-wide ClinVar P/LP intersection ==")
    plp = {}
    with gzip.open("data/clinvar.vcf.gz", "rt") as f:
        for line in f:
            if line.startswith("#"):
                continue
            fl = line.split("\t")
            if any(s in fl[7] for s in ("CLNSIG=Pathogenic", "CLNSIG=Likely_pathogenic", "CLNSIG=Pathogenic/Likely_pathogenic")):
                plp[(fl[0], fl[1], fl[3], fl[4])] = fl[7]
    with gzip.open(vcf_path, "rt") as f:
        for line in f:
            if line.startswith("#"):
                continue
            fl = line.rstrip("\n").split("\t")
            if fl[6] != "PASS":
                continue
            for alt in fl[4].split(","):
                info = plp.get((fl[0], fl[1], fl[3], alt))
                if info:
                    d = dict(kv.split("=", 1) for kv in info.split(";") if "=" in kv)
                    print(f"  {fl[0]}:{fl[1]} {fl[3]}>{alt} GT={fl[9].split(':')[0]} {d.get('GENEINFO','?')[:30]} {d.get('CLNDN','?')[:60]}")


def aneuploidy_scan(vcf_path):
    print("== per-chromosome depth / het-BAF scan (mosaic aneuploidy) ==")
    chroms = [str(i) for i in range(1, 23)] + ["X"]
    baf = {c: [] for c in chroms}
    dp = {c: [0, 0] for c in chroms}
    with gzip.open(vcf_path, "rt") as f:
        for line in f:
            if line.startswith("#"):
                continue
            fl = line.rstrip("\n").split("\t")
            if fl[0] not in baf or fl[6] != "PASS":
                continue
            smp = fl[9].split(":")
            try:
                ad = [int(x) for x in smp[1].split(",")]
            except (ValueError, IndexError):
                continue
            d = sum(ad)
            if d < 15:
                continue
            dp[fl[0]][0] += d
            dp[fl[0]][1] += 1
            if smp[0] == "0/1" and len(ad) == 2 and d <= 80:
                baf[fl[0]].append(ad[1] / d)
    base = statistics.median(dp[c][0] / max(dp[c][1], 1) for c in chroms[:22])
    for c in chroms:
        if not baf[c]:
            continue
        mdp = dp[c][0] / max(dp[c][1], 1)
        sd = statistics.pstdev(baf[c])
        flag = " <-- deviates" if (abs(statistics.median(baf[c]) - 0.5) > 0.03 or sd > 0.13 or abs(mdp / base - 1) > 0.06) else ""
        print(f"  chr{c:<3} nhet={len(baf[c]):>7} medBAF={statistics.median(baf[c]):.3f} sdBAF={sd:.3f} meanDP={mdp:6.1f} ({mdp/base:+.0%} vs baseline){flag}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--vcf", default="data/WGS_EX2312012_HGWCNDSX7.vcf.gz")
    a = ap.parse_args()
    import os, subprocess
    if not os.path.exists("data/clinvar.vcf.gz"):
        for suf in ("", ".tbi"):
            subprocess.run(["curl", "-sL", "-o", f"data/clinvar.vcf.gz{suf}", CLINVAR + suf], check=True)
    pat = pysam.TabixFile(a.vcf)
    cv = pysam.TabixFile("data/clinvar.vcf.gz")
    gene_scan(pat, cv)
    genomewide_plp(a.vcf)
    aneuploidy_scan(a.vcf)


if __name__ == "__main__":
    main()
