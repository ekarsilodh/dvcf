# VCF Trio Analyzer – Required Files & How to Run

## Dependencies
1. bedtools
2. python

---

## To install bedtools, use:

1. for macOS:

```bash
brew install bedtools
```

2. for linux:

```bash
sudo apt install bedtools
```

---

## Required Files
1. **Input VCF file (trio VCF)**
2. `clinvar_SV_clean.bed`
3. `ClinGen_recurrent_CNV-hg38.bed`
4. `ClinGen_haploinsufficiency_gene_GRCh38.bed`
5. `ClinGen_triplosensitivity_gene_GRCh38.bed`
6. `hg38_exons.bed`
7. `hg38_genes.bed`

---

## How to Run

```bash
python vcf_analyzer.py --vcf /path/to/vcf/file --out /path/to/output/directory
```

**Note:**
- `--vcf` is **required**
- `--out` is optional (defaults to `out/`)

---

## Example

```bash
python vcf_analyzer.py --vcf DellyVariation.vcf
```
