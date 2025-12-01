# 🧬 Delly Trio Structural Variant Explorer

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.17781941.svg)](https://doi.org/10.5281/zenodo.17781941)

An interactive Streamlit dashboard for trio-based structural variant (SV) analysis, Mendelian validation, de novo detection, and clinical annotation using DELLY-generated VCF files.

Visit the app: 🔗 **https://dellyvcftrio.streamlit.app/**

> **Note:** Due to Streamlit Cloud limitations, **annotation steps (BEDTools-based)** cannot run in the online app.  
> For the full, annotation-enabled pipeline, visit the complete repository:  
> 🔗 **https://github.com/ekarsilodh/Delly-SV-Trios**

---

## 🚀 Overview
The **Delly Trio SV Explorer** is a visual, talk-ready dashboard for exploring DELLY multi-sample SV VCFs.  
It provides:

- **Automatic trio role inference**
- **Sex calling using chrX heterozygosity**
- **Mendelian violation analysis**
- **De novo structural variant detection**
- **SV-type and chromosome-level visualizations**
- **Exportable tables and interactive filters**

Ideal for genomics teaching, diagnostics, rapid SV interpretation, and trio-based variant exploration.

---

## 🎯 Key Features
### 👨‍👩‍👧 Trio Inference
- Finds best child–mother–father assignment  
- Computes violations across all permutations  
- Flags invalid or unrelated trios  

### 🧬 Sex Calling
- chrX heterozygosity rate–based inference  
- Displays genotype composition of X-linked SVs  

### 💥 De Novo SV Discovery
- Identifies variants where both parents are `0/0` and child is `!= 0/0`  
- Filterable by chromosome & SVTYPE  
- Downloadable tables  

### 📊 Visual Analytics
- SVTYPE distribution with counts  
- Chromosome-level density plots
  
### 🧪 Sample VCF Included
Test the app without uploading your own data.

---

## 📦 Installation
Clone the repository:

```bash
git clone https://github.com/ekarsilodh/dvcf.git
cd dvcf
```

Install dependencies:

```bash
pip install -r requirements.txt
```

---

## ▶️ Running Locally
```bash
streamlit run src/app_local.py
```

Open in your browser at:  
**http://localhost:8501/**

---

## ⚠️ Streamlit Cloud Limitations
Streamlit Cloud does **not** allow installation of `bedtools` or system binaries.

Therefore:
- **Annotation**
- **ClinGen gene mapping**
- **Pathogenicity tiering**

cannot run online.

To use the **full version** of the project with annotation support:

👉 **Visit the complete repository here:**  
https://github.com/ekarsilodh/Delly-SV-Trios

---

## 🗂 Directory Structure
```
project/
│
├── src/
│   ├── app.py
│   ├── app_local.py
│   ├── vcf_analyzer.py
│   ├── vcf_analyzer_local.py
│   ├── plot.py
│   └──data/
│      └── sample_trio.vcf
│
├── assets/
│   ├── logo.png
│   ├── trio_banner.png
│   ├── Pipeline.png
│   └── overview_illustration.png
│
├── databases/ (local use only)
│   ├── hg38_genes.bed
│   ├── hg38_exons.bed
│   ├── ClinGen_haploinsufficiency_gene_GRCh38.bed
│   ├── ClinGen_triplosensitivity_gene_GRCh38.bed
│   └── ClinGen_recurrent_CNV_GRCh38.bed
│
└── requirements.txt
```

---

## 🎨 Streamlit Theme (config.toml)
```toml
[theme]
base="dark"
primaryColor="#ff4d4d"
backgroundColor="#1a1b26"
secondaryBackgroundColor="#24283b"
textColor="#c0caf5"
```

---

## 📥 Usage Workflow
1. Upload your DELLY-generated multi-sample VCF  
2. Let the app infer trio roles  
3. Inspect SV distributions  
4. Identify de novo events  
5. Download filtered tables  
6. (Local only) Run annotation and pathogenicity tiering  

---

## 👨‍💻 Author
**Ekarsi Lodh**  
MSc Bioinformatics
College of Medicine and Health
University of Birmingham  

---

## 📜 License
MIT License — free to use, modify, and extend.

---

## If you find this tool useful and have used it in your work, please ⭐ star the repository and cite the following:
E. Lodh, “Streamlit-Based Interactive Trio VCF Analyzer for Structural Variant Interpretation”. Zenodo, Dec. 01, 2025. doi: [10.5281/zenodo.17781941](https://doi.org/10.5281/zenodo.17781941).
