#!/usr/bin/env python3
# -*- coding: utf-8 -*-

#Author: Ekarsi Lodh

"""
VCF Trio Analyzer (Python + pybedtools)
--------------------------------

Pipeline:
1. Python VCF parsing:
    - Reads VCF
    - Extracts CHROM..SVLEN, and GT for the LAST THREE samples
    - Fast parsing

2. Python:
    - Infers trio roles (child / parent1 / parent2)
    - Determines ORIGINAL COLUMN numbers (exact like $16, $17, $18 etc)
    - Runs Python de-novo finder 
    - Summarizes:

       a) De-novo SV calls
       b) All trio calls (counts SVTYPE, etc.)

3. bedtools (via pybedtools):
    - Annotates SVs against:
       * hg38 gene bodies
       * hg38 exons
       * ClinGen HI genes
       * ClinGen TS genes
       * ClinGen recurrent CNVs

Outputs:
    - trio_table.tsv              : Parsed trio table with GTs and SV annotations.
    - svtype_counts.tsv           : Global count of SVTYPEs.
    - de_novo_calls.txt / .tsv    : De novo candidate SVs (0/0 parents, child non-ref).
    - clinically_prioritised_SVs.tsv : Child's SVs passing size/region filters.
    - sv_for_bedtools.bed         : BED intervals used for all bedtools annotations.
    - bedtools_hg38_genes_raw.tsv : Raw SV × gene overlaps.
    - hg38_exons_raw.tsv / _counts.tsv : Raw + counts for exon overlaps.
    - clingen_hi_raw.tsv / _counts.tsv : Raw + counts for ClinGen HI genes.
    - clingen_ts_raw.tsv / _counts.tsv : Raw + counts for ClinGen TS genes.
    - clingen_recurrent_cnv_raw.tsv / _counts.tsv : Raw + counts for recurrent CNVs.
"""

import argparse
import tempfile
import pandas as pd
import os
import sys
from itertools import permutations
from typing import Optional, Dict


############################################
# CONSTANTS / PATHS FOR BEDTOOLS ANNOTATION
############################################

HG38_GENES_BED = "../databases/hg38_genes.bed"
HG38_EXONS_BED = "../databases/hg38_exons.bed"

CLINGEN_HI_GENES_BED = "../databases/ClinGen_haploinsufficiency_gene_GRCh38.bed"
CLINGEN_TS_GENES_BED = "../databases/ClinGen_triplosensitivity_gene_GRCh38.bed"
CLINGEN_RECURRENT_CNV_BED = "../databases/ClinGen_recurrent_CNV_GRCh38.bed"

CLINVAR_SV_BED = "../databases/clinvar_SV_clean.bed"  # Pathogenic / Likely pathogenic SVs (GRCh38) only

# 1. ARGUMENT PARSING

def parse_args():
    """
    Parse command-line arguments.

    Returns a Namespace with:
      vcf: path to input VCF
      out: output directory (default: "out")
    """
    p = argparse.ArgumentParser()
    p.add_argument("--vcf", required=True)
    p.add_argument("--out", default="../out")
    return p.parse_args()


############################################
# HELPER: GET SAMPLE NAMES FROM VCF HEADER
############################################

def get_vcf_samples(vcf_path):
    """
    Read the VCF header, return list of sample names (columns 10+).
    """
    with open(vcf_path) as f:
        for line in f:
            if line.startswith("#CHROM"):
                fields = line.strip().split("\t")
                return fields[9:]
    return []


def choose_trio_samples(vcf_path):
    """
    Interactively ask user which 3 samples form the trio.

    Returns a list [sample1, sample2, sample3].
    """
    all_samples = get_vcf_samples(vcf_path)

    print("\n===== SAMPLES IN VCF =====")
    for i, s in enumerate(all_samples, start=1):
        print(f"{i}. {s}")

    print("\nEnter the names of the 3 samples you want to analyse,")
    print("exactly as they appear above.")

    chosen = []
    while len(chosen) < 3:
        name = input(f"Sample {len(chosen)+1} name: ").strip()
        if name not in all_samples:
            print(f"  → '{name}' is not a valid sample name in this VCF. Please try again.")
            continue
        if name in chosen:
            print("  → You already selected that sample; please choose a different one.")
            continue
        chosen.append(name)

    return chosen


# 2. AWK SCRIPT FOR FIRST PASS
# (kept for reference, not used any more; actual parsing is now done in Python)

AWK_PARSE = r'''
function extract_gt(sampleField,   a,b) {
    n = split(sampleField, a, ":");
    # return genotype only
    split(a[1], b, ":");
    return a[1];
}

BEGIN { OFS="\t"; }

/^##/ { next }

/^#CHROM/ {
    for (i=1; i<=NF; i++) H[i]=$i;
    n=NF;

    # Find the three sample columns matching the requested names.
    # The variables sample1, sample2 and sample3 are passed in from Python via -v.
    s1 = s2 = s3 = -1;
    for (i=10; i<=NF; i++) {
        if (H[i] == sample1)      s1 = i;
        else if (H[i] == sample2) s2 = i;
        else if (H[i] == sample3) s3 = i;
    }

    # Fallback: if any of the requested names were not found,
    # default to using the last three sample columns.
    if (s1 < 0 || s2 < 0 || s3 < 0) {
        s1 = n-2; s2 = n-1; s3 = n;
    }

    print "CHROM","POS","REF","ALT","QUAL","FILTER","INFO",
          H[s1]"_GT",H[s2]"_GT",H[s3]"_GT", "SVTYPE", "END", "SVLEN";
    next;
}


/^[^#]/ && $7 == "PASS" {
    # parse INFO
    SVTYPE=""; ENDV=""; SVLEN="";
    split($8, infoA, ";");
    for (i in infoA){
        split(infoA[i], kv, "=");
        if (kv[1]=="SVTYPE") SVTYPE=kv[2];
        else if (kv[1]=="END") ENDV=kv[2];
        else if (kv[1]=="SVLEN") SVLEN=kv[2];
    }

    gt1 = extract_gt($(s1));
    gt2 = extract_gt($(s2));
    gt3 = extract_gt($(s3));

    print $1,$2,$4,$5,$6,$7,$8,gt1,gt2,gt3,SVTYPE,ENDV,SVLEN;
}
'''


# 2b. Python VCF parser (replacing AWK first pass)

def parse_vcf_trio_table(vcf_path: str, trio_samples, pass_only: bool = True) -> pd.DataFrame:
    """
    Pure-Python VCF parser that reproduces the output of the AWK first pass.

    It writes a table with columns:
      CHROM, POS, REF, ALT, QUAL, FILTER, INFO,
      <sample1>_GT, <sample2>_GT, <sample3>_GT, SVTYPE, END, SVLEN

    where sample1/2/3 are the chosen trio sample names. By default it
    keeps only records with FILTER == "PASS", matching the AWK logic.
    """
    trio_samples = list(trio_samples)
    if len(trio_samples) != 3:
        raise ValueError("trio_samples must contain exactly 3 sample names.")

    rows = []
    header_parsed = False
    s1 = s2 = s3 = None
    sample_header = None

    with open(vcf_path) as f:
        for line in f:
            if line.startswith("##"):
                continue

            if line.startswith("#CHROM"):
                # Header line with all column names including FORMAT + samples
                sample_header = line.strip().split("\t")
                n = len(sample_header)

                # map requested sample names to column indices (0-based)
                name_to_idx = {name: i for i, name in enumerate(sample_header)}

                for s in trio_samples:
                    if s not in name_to_idx:
                        # we'll fall back to last 3 samples below
                        break
                else:
                    # all requested names present
                    s1 = name_to_idx[trio_samples[0]]
                    s2 = name_to_idx[trio_samples[1]]
                    s3 = name_to_idx[trio_samples[2]]

                # fallback to last three columns if any name not found
                if s1 is None or s2 is None or s3 is None:
                    s1, s2, s3 = n - 3, n - 2, n - 1

                # Build header row exactly as AWK did
                h_s1 = sample_header[s1]
                h_s2 = sample_header[s2]
                h_s3 = sample_header[s3]
                header = [
                    "CHROM",
                    "POS",
                    "REF",
                    "ALT",
                    "QUAL",
                    "FILTER",
                    "INFO",
                    f"{h_s1}_GT",
                    f"{h_s2}_GT",
                    f"{h_s3}_GT",
                    "SVTYPE",
                    "END",
                    "SVLEN",
                ]
                rows.append(header)
                header_parsed = True
                continue

            if not header_parsed:
                # safety: skip anything before #CHROM
                continue

            if line.startswith("#"):
                continue

            parts = line.strip().split("\t")
            if len(parts) < 10:
                continue

            # FILTER == PASS only (AWK: /^[^#]/ && $7 == "PASS")
            if pass_only and parts[6] != "PASS":
                continue

            chrom = parts[0]
            pos = parts[1]
            ref = parts[3]
            alt = parts[4]
            qual = parts[5]
            flt = parts[6]
            info = parts[7]

            # parse INFO for SVTYPE, END, SVLEN
            svtype = ""
            endv = ""
            svlen = ""
            for field in info.split(";"):
                if "=" in field:
                    key, val = field.split("=", 1)
                else:
                    key, val = field, ""
                if key == "SVTYPE":
                    svtype = val
                elif key == "END":
                    endv = val
                elif key == "SVLEN":
                    svlen = val

            def extract_gt(sample_field: str) -> str:
                if not sample_field:
                    return "./."
                # genotype is the first colon-separated item
                return sample_field.split(":", 1)[0]

            gt1 = extract_gt(parts[s1])
            gt2 = extract_gt(parts[s2])
            gt3 = extract_gt(parts[s3])

            rows.append(
                [
                    chrom,
                    pos,
                    ref,
                    alt,
                    qual,
                    flt,
                    info,
                    gt1,
                    gt2,
                    gt3,
                    svtype,
                    endv,
                    svlen,
                ]
            )

    if not rows:
        # no variants
        return pd.DataFrame(
            columns=[
                "CHROM",
                "POS",
                "REF",
                "ALT",
                "QUAL",
                "FILTER",
                "INFO",
                f"{trio_samples[0]}_GT",
                f"{trio_samples[1]}_GT",
                f"{trio_samples[2]}_GT",
                "SVTYPE",
                "END",
                "SVLEN",
            ]
        )

    # first row is header
    columns = rows[0]
    data_rows = rows[1:]
    return pd.DataFrame(data_rows, columns=columns)


# 3. GENOTYPE CODING

def gt_to_code(gt: Optional[str]) -> Optional[int]:
    """Convert GT to simple 0/1/2 code."""
    if gt in (None, ".", "./.", ".|."):
        return None
    if gt in ("0/0", "0|0", "0"):
        return 0
    if gt in ("0/1", "1/0", "0|1", "1|0", "1/2", "2/1", "0/2", "2/0"):
        return 1
    if gt in ("1/1", "1|1", "1"):
        return 2
    return None


# 4. TRIO ROLE INFERENCE by MENDELIAN VIOLATION

def mendelian_violation(child, p1, p2):
    """
    - Parent1 genotype = 0/0  → encoded as 0
    - Parent2 genotype = 0/0  → encoded as 0
    - Child genotype  != 0/0  → encoded as anything non-zero

    Returns True iff this pattern is observed, False otherwise.
    """
    # ignore if any genotype is missing
    if child is None or p1 is None or p2 is None:
        return False

    # both parents 0/0, child not 0/0
    return (p1 == 0 and p2 == 0 and child != 0)


def infer_trio(df, samples):
    """
    Given three samples, try all child/parent assignments and
    pick the one with the fewest Mendelian violations.
    """
    best = None
    roles = None
    trio_stats = []

    for child, p1, p2 in permutations(samples, 3):
        vio = 0
        for _, row in df.iterrows():
            c = gt_to_code(row[f"{child}_GT"])
            m = gt_to_code(row[f"{p1}_GT"])
            f = gt_to_code(row[f"{p2}_GT"])
            if mendelian_violation(c, m, f):
                vio += 1

        trio_stats.append({
            "child": child, 
            "parent1": p1, 
            "parent2": p2, 
            "violations": vio
        })

        if best is None or vio < best:
            best = vio
            roles = {"child": child, "parent1": p1, "parent2": p2}
    # sanity check — if even the best assignment has > 4 violations,
    # this trio probably isn't a valid family.
    if best is not None and best > 4:
        return None, trio_stats, best

    return roles, trio_stats, best


# 5. SEX INFERENCE FROM ChrX HET

def infer_sex_from_chrX(
    df: pd.DataFrame,
    samples,
    het_threshold: float = 0.1,
    min_sites: int = 5,
):
    """
    Sex inference based on heterozygosity on chrX:

      - For each sample:
          H = (# het genotypes on X) / (# non-missing genotypes on X)
      - If H < het_threshold → "M"
      - Else → "F"
    Returns:
      sex_calls: {sample: "M"/"F"/"unknown"}
      het_rates: {sample: float or None}
    """

    sex_calls: Dict[str, str] = {s: "unknown" for s in samples}
    het_rates: Dict[str, Optional[float]] = {s: None for s in samples}

    # Normalise CHROM (strip "chr") and select X only
    chrom_series = df["CHROM"].astype(str)
    chrom_nochr = chrom_series.str.replace("^chr", "", case=False, regex=True)
    mask_x = chrom_nochr == "X"
    df_x = df[mask_x]

    if df_x.empty:
        # no chrX sites
        return sex_calls, het_rates

    for s in samples:
        col = f"{s}_GT"
        if col not in df_x.columns:
            continue

        col_series = df_x[col].fillna(".")

        non_missing_mask = ~col_series.isin([".", "./.", ".|."])
        n_non_missing = int(non_missing_mask.sum())
        if n_non_missing < min_sites:
            het_rates[s] = None
            continue

        het_mask = col_series.isin(["0/1", "1/0", "0|1", "1|0"])
        H = float(het_mask.sum()) / float(n_non_missing)

        het_rates[s] = H
        if H < het_threshold:
            sex_calls[s] = "M"
        else:
            sex_calls[s] = "F"

    return sex_calls, het_rates

# 6. DE NOVO DETECTION (Python implementation)

def run_awk_denovo(vcf, outpath, outpath1, parent1_col, parent2_col, child_col):
    """
    Python implementation of the original AWK-based de novo detector.

    Parents must be 0/0, child must be non-ref and not missing.
    Writes two identical tab-separated files (outpath, outpath1) with columns:
      CHROM, START_POS, END_POS, ID, REF, ALT, SV_TYPE, LENGTH_bp,
      Paired_end_PE, Split_end_SR, Mother_GT, Father_GT, Child_GT
    """
    # parent*_col and child_col are 1-based VCF column numbers.
    header_written = False
    rows = []

    with open(vcf) as f:
        for line in f:
            if line.startswith("#"):
                continue

            fields = line.rstrip("\n").split("\t")
            if len(fields) < max(parent1_col, parent2_col, child_col):
                continue

            # Sample fields (1-based -> 0-based indices)
            p1_field = fields[parent1_col - 1]
            p2_field = fields[parent2_col - 1]
            c_field = fields[child_col - 1]

            # De novo genotype pattern:
            # parents 0/0, child non-0/0 and not missing ("./.")
            if p1_field.startswith("0/0") and p2_field.startswith("0/0"):
                if c_field.startswith("0/0") or c_field.startswith("./."):
                    continue
            else:
                continue

            # positions and basic fields
            chrom = fields[0]
            pos_str = fields[1]
            vid = fields[2]
            ref = fields[3]
            alt = fields[4]
            info = fields[7]

            try:
                pos = int(pos_str)
            except ValueError:
                continue

            # Parse INFO for SVTYPE, END, PE, SR
            svtype = ""
            endpos = pos
            pe = 0
            sr = 0

            for field in info.split(";"):
                if "=" in field:
                    key, val = field.split("=", 1)
                else:
                    key, val = field, ""
                if key == "SVTYPE":
                    svtype = val
                elif key == "END":
                    try:
                        endpos = int(val)
                    except ValueError:
                        endpos = pos
                elif key == "PE":
                    try:
                        pe = int(val)
                    except ValueError:
                        pe = 0
                elif key == "SR":
                    try:
                        sr = int(val)
                    except ValueError:
                        sr = 0

            length_bp = endpos - pos

            def extract_gt(sample_field: str) -> str:
                if not sample_field:
                    return "./."
                return sample_field.split(":", 1)[0]

            mother_gt = extract_gt(p1_field)
            father_gt = extract_gt(p2_field)
            child_gt = extract_gt(c_field)

            if not header_written:
                header = [
                    "CHROM",
                    "START_POS",
                    "END_POS",
                    "ID",
                    "REF",
                    "ALT",
                    "SV_TYPE",
                    "LENGTH_bp",
                    "Paired_end_PE",
                    "Split_end_SR",
                    "Mother_GT",
                    "Father_GT",
                    "Child_GT",
                ]
                rows.append("\t".join(header))
                header_written = True

            row = [
                chrom,
                str(pos),
                str(endpos),
                vid,
                ref,
                alt,
                svtype,
                str(length_bp),
                str(pe),
                str(sr),
                mother_gt,
                father_gt,
                child_gt,
            ]
            rows.append("\t".join(row))

    # Write out the two identical files
    content = "\n".join(rows) + ("\n" if rows else "")
    with open(outpath, "w") as outfh:
        outfh.write(content)
    with open(outpath1, "w") as outfh:
        outfh.write(content)


# 7. SV LENGTH, SIZE PRIORITY

def add_sv_length_and_size_category(df: pd.DataFrame) -> pd.DataFrame:
    """
    Adds:
      - POS_INT, END_INT, SVLEN_INT (numeric)
      - SVLEN_BP      : final length used (abs SVLEN or END - POS)
      - SIZE_PRIORITY : 'high', 'moderate', 'low', or 'unknown'
    This is applied to ALL SVs, but is most meaningful for DEL/DUP.
    """
    # numeric POS
    df["POS_INT"] = pd.to_numeric(df["POS"], errors="coerce")

    # numeric END (from parsed TSV)
    if "END" in df.columns:
        df["END_INT"] = pd.to_numeric(df["END"], errors="coerce")
    else:
        df["END_INT"] = df["POS_INT"]

    # numeric SVLEN
    if "SVLEN" in df.columns:
        df["SVLEN_INT"] = pd.to_numeric(df["SVLEN"], errors="coerce")
    else:
        df["SVLEN_INT"] = None

    # final SVLEN_BP to use
    def _choose_length(row):
        if pd.notna(row.get("SVLEN_INT")):
            try:
                return abs(int(row["SVLEN_INT"]))
            except Exception:
                pass
        # fallback: END - POS
        try:
            if pd.notna(row.get("END_INT")) and pd.notna(row.get("POS_INT")):
                return abs(int(row["END_INT"]) - int(row["POS_INT"]))
        except Exception:
            pass
        return None

    df["SVLEN_BP"] = df.apply(_choose_length, axis=1)

    # categorize size
    def _size_priority(L):
        if L is None or pd.isna(L):
            return "unknown"
        L = int(L)
        if L >= 500000:
            return "high"
        if L >= 100000:
            return "moderate"
        if L >= 10000:
            return "low"
        return "unknown"

    df["SIZE_PRIORITY"] = df["SVLEN_BP"].map(_size_priority)

    return df


# 8. INHERITANCE ANNOTATION FOR ALL VARIANTS

def annotate_inheritance(df: pd.DataFrame, child: str, p1: str, p2: str) -> pd.DataFrame:
    """
    Label inheritance for ALL variants based on trio genotypes:

      - de_novo
      - inherited_p1
      - inherited_p2
      - inherited_both
      - reference_or_absent_child
      - unknown

    This does NOT filter anything; it just adds an INHERITANCE column.
    """
    child_col = f"{child}_GT"
    p1_col = f"{p1}_GT"
    p2_col = f"{p2}_GT"

    def classify(row):
        c = gt_to_code(row.get(child_col))
        m = gt_to_code(row.get(p1_col))
        f = gt_to_code(row.get(p2_col))

        # Missing
        if c is None or m is None or f is None:
            return "unknown"

        # Child reference or absent
        if c == 0:
            return "reference_or_absent_child"

        # both parents 0 => de novo
        if m == 0 and f == 0:
            return "de_novo"

        # child alt, mother alt, father ref
        if m > 0 and f == 0:
            return "inherited_p1"

        # child alt, father alt, mother ref
        if f > 0 and m == 0:
            return "inherited_p2"

        # both parents alt, child alt
        if m > 0 and f > 0:
            return "inherited_both"

        return "unknown"

    df = df.copy()
    df["INHERITANCE"] = df.apply(classify, axis=1)
    return df


############################################
# MAIN SCRIPT
############################################

def main():
    args = parse_args()
    os.makedirs(args.out, exist_ok=True)
    
    print("[INFO] Reading vcf...")

    # ---------------------------------------------
    # ASK USER WHICH 3 SAMPLES TO ANALYSE (by name)
    # ---------------------------------------------
    trio_samples = choose_trio_samples(args.vcf)

    print("\n[INFO] Selected samples:", ", ".join(trio_samples))
    print("[INFO] Starting initial parsing using pure Python...")

    # ------------------------------
    # RUN PYTHON-BASED FIRST PASS
    # ------------------------------
    parsed_tsv = os.path.join(args.out, "trio_table.tsv")

    print("[INFO] Running primary Python parsing to extract trio genotype table...")

    parse_vcf_trio_table(args.vcf, trio_samples).to_csv(parsed_tsv, sep="\t", index=False)


    print("[INFO] Genotype table created:", parsed_tsv)

    # ------------------------------
    # LOAD PARSED TSV
    # ------------------------------
    print("[INFO] Loading parsed trio table...")

    df = pd.read_csv(parsed_tsv, sep="\t", dtype=str)

    # extract sample names from columns ending with _GT
    samples = [c.replace("_GT", "") for c in df.columns if c.endswith("_GT")]

    # ------------------------------
    # INFER TRIO ROLES + TRIO STATS
    # ------------------------------
    print("[INFO] Trio table loaded. Inferring trio roles...")

    roles, trio_stats = infer_trio(df, samples)
    print("[INFO] Determining who is child, father, mother...")

    child = roles["child"]
    p1 = roles["parent1"]
    p2 = roles["parent2"]

    min_vio = min(ts["violations"] for ts in trio_stats)

    print("\n===== TRIO ROLE SEARCH (Mendelian violations) =====")
    for ts in trio_stats:
        print(
            f"Assuming child={ts['child']}, "
            f"parent1={ts['parent1']}, parent2={ts['parent2']} "
            f"→ violations={ts['violations']}"
        )
    print(
        f"\nSelected trio = child={child}, parent1={p1}, parent2={p2}\n "
        f"because this assignment had the smallest number of Mendelian\n "
        f"violations ({min_vio}) among all permutations."
    )
    
    print("\n===== TRIO ROLES =====")
    print(f"Child    = {child}")
    print(f"Parent1  = {p1}")
    print(f"Parent2  = {p2}")

    # ---------------------------------------
    # INFER SEX FROM chrX  (child / parents)
    # ---------------------------------------
    print("\n[INFO] Performing sex inference from chrX heterozygosity...")

    sex_calls, x_het_rates = infer_sex_from_chrX(df, samples)

    child_sex = sex_calls.get(child, "unknown")
    p1_sex = sex_calls.get(p1, "unknown")
    p2_sex = sex_calls.get(p2, "unknown")

    father = None
    mother = None

    parent_sexes = {p1: p1_sex, p2: p2_sex}
    male_parents = [s for s, sx in parent_sexes.items() if sx == "M"]
    female_parents = [s for s, sx in parent_sexes.items() if sx == "F"]

    if male_parents:
        father = male_parents[0]
    if female_parents:
        mother = female_parents[0]

    print("\n===== SEX INFERENCE (chrX) =====")
    for s in samples:
        hr = x_het_rates.get(s, None)
        hr_str = "NA" if hr is None else f"{hr:.4f}"
        print(f"{s}: sex={sex_calls.get(s, 'unknown')}, X_het_rate={hr_str}")
    print(f"\nChild sex : {child} → {child_sex}")
    print(f"Father    : {father if father else 'undetermined'}")
    print(f"Mother    : {mother if mother else 'undetermined'}")

    # ------------------------------
    # FIND ORIGINAL COLUMN NUMBERS
    # ------------------------------
    header_samples = []
    with open(args.vcf) as f:
        for line in f:
            if line.startswith("#CHROM"):
                header_samples = line.strip().split("\t")
                break

    sample_to_col = {}
    for idx, name in enumerate(header_samples):
        if idx >= 9:
            sample_to_col[name] = idx + 1 

    child_col = sample_to_col[child]
    p1_col = sample_to_col[father]
    p2_col = sample_to_col[mother]

    print("\n===== VCF COLUMN NUMBERS =====")
    print(f"Child col   = {child_col}")
    print(f"Father col = {p1_col}")
    print(f"Mother col = {p2_col}")
    
    # ------------------------------
    # RUN AWK DE-NOVO DETECTOR (now Python)
    # ------------------------------
    print("\n[INFO] Running de novo variant discovery ...")

    denovo_file_txt = os.path.join(args.out, "de_novo_calls.txt")
    denovo_file_tsv = os.path.join(args.out, "de_novo_calls.tsv")
    run_awk_denovo(args.vcf, denovo_file_txt, denovo_file_tsv, p1_col, p2_col, child_col)
    print("[INFO] De novo variant detection completed.")
    
    # ------------------------------------------------
    # ANNOTATE SV LENGTH, SIZE PRIORITY & INHERITANCE
    # ------------------------------------------------
    print("[INFO] Getting SV Type Counts ...")
    print("[INFO] Annotating variants with ClinVar & ClinGen datasets...")

    df = add_sv_length_and_size_category(df)
    df = annotate_inheritance(df, child=child, p1=p1, p2=p2)
    
    # child non-ref subset (these are the variants actually present in child)
    child_gt_col = f"{child}_GT"
    alt_gts = {"0/1", "1/0", "0|1", "1|0", "1/1", "1|1", "1"}
    df_child_alt = df[df[child_gt_col].isin(alt_gts)].copy()
    
    # ------------------------------
    # COUNT ALL SVTYPES
    # ------------------------------
    sv_counts = df["SVTYPE"].value_counts().reset_index()
    sv_counts.columns = ["SVTYPE", "COUNT"]
    sv_counts.to_csv(os.path.join(args.out, "svtype_counts.tsv"), sep="\t", index=False)
    
    print("\n===== SVTYPE COUNTS =====")
    print(str(sv_counts))
    

    # Inheritance category counts
    inheritance_counts = df["INHERITANCE"].value_counts(dropna=False).reset_index()
    inheritance_counts.columns = ["INHERITANCE", "COUNT"]
    
    inh_tsv = os.path.join(args.out, "inheritance_counts.tsv")
    inheritance_counts.to_csv(inh_tsv, sep="\t", index=False)

if __name__ == "__main__":
    main()
