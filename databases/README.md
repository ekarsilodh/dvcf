# Databases / Annotation Files

This folder is intentionally left **empty**.  
It acts as a placeholder to show users **which BED files they must provide**
to run the `vcf_analyzer.py` pipeline.

## Required BED Files

The pipeline expects the following annotation files, all in **GRCh38 (hg38)**
coordinates:

### 1. `hg38_genes.bed`
- Gene intervals for GRCh38
- Typically downloaded or generated from:
  - Ensembl Biomart
  - UCSC Table Browser
  - GENCODE (recommended)

### 2. `hg38_exons.bed`
- Exon intervals for GRCh38
- Can be derived from:
  - Ensembl GTF → converted to BED
  - GENCODE comprehensive annotation

### 3. `ClinGen_haploinsufficiency_gene_GRCh38.bed`
- ClinGen curated **haploinsufficient (HI)** genes
- Download from ClinGen Genome Dosage Map

### 4. `ClinGen_triplosensitivity_gene_GRCh38.bed`
- ClinGen curated **triplosensitive (TS)** genes

### 5. `ClinGen recurrent CNV-hg38.bed`
- ClinGen recurrent CNV regions
- Used for SV interpretation based on known pathogenic hotspots

### 6. `clinvar_SV_clean.bed`
- ClinVar **pathogenic / likely pathogenic** structural variants
- Must be filtered and converted to BED
- Suggested workflow:
  1. Download ClinVar VCF (GRCh38)
  2. Filter for `CLNSIG=Pathogenic, Likely_pathogenic`
  3. Convert to BED using:
     ```bash
     bcftools query -f '%CHROM\t%POS\t%INFO/END\t%ID\n' file.vcf > clinvar_SV_clean.bed
     ```

---

## Folder Usage

Place your BED files here when running the pipeline locally:

databases/
├── hg38_genes.bed
├── hg38_exons.bed
├── ClinGen_haploinsufficiency_gene_GRCh38.bed
├── ClinGen_triplosensitivity_gene_GRCh38.bed
├── ClinGen recurrent CNV-hg38.bed
└── clinvar_SV_clean.bed

