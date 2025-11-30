# app.py
import os
import tempfile
import subprocess
from typing import Dict, List, Optional, Tuple

import pandas as pd
import streamlit as st
from PIL import Image

import vcf_analyzer_local as vca 
import plotly.express as px

# Nice qualitative palettes (Seaborn-like via Plotly)
PALETTES = {
    "set2": px.colors.qualitative.Set2,
    "set3": px.colors.qualitative.Set3,
    "pastel": px.colors.qualitative.Pastel,
    "dark": px.colors.qualitative.Dark24,
    "bold": px.colors.qualitative.Bold,
}


def plot_count_bar(df, x, y, title="", palette="set2"):
    """
    Categorical bar chart with:
      - colourful bars (palette by category)
      - count labels on top
      - dark theme
    """
    colors = PALETTES.get(palette, px.colors.qualitative.Set2)

    fig = px.bar(
        df,
        x=x,
        y=y,
        text=y,
        color=x,  # colour by category
        color_discrete_sequence=colors,
        title=title,
        template="plotly_dark",
    )

    fig.update_traces(
        textposition="outside",
        hovertemplate=f"{x}: %{{x}}<br>{y}: %{{y}}<extra></extra>",
        marker=dict(line=dict(width=1, color="rgba(255,255,255,0.4)")),
    )

    fig.update_layout(
        font=dict(size=14),
        margin=dict(t=60, b=40, l=40, r=20),
        xaxis=dict(title="", tickangle=0),
        yaxis=dict(title="Count"),
        height=450,
        legend_title_text="",
    )

    return fig


def chrom_order_key(chrom: str) -> int:
    """
    Map chromosome label to a sortable integer:
    1..23 -> 1..23
    X -> 24
    Y -> 25
    MT/M -> 26
    Everything else -> 999
    Works with or without 'chr' prefix.
    """
    if chrom is None:
        return 999
    c = str(chrom).strip()
    c = c.replace("chr", "").replace("CHR", "")
    c = c.upper()

    # numeric chromosomes
    if c.isdigit():
        return int(c)

    if c == "X":
        return 24
    if c == "Y":
        return 25
    if c in {"MT", "M"}:
        return 26

    return 999


# -------------------------
# Page configuration
# -------------------------
st.set_page_config(
    page_title="VCF Trio SV Analyzer",
    page_icon="🧬",
    layout="wide",
)

# -------------------------
# Custom styling (Tokyo-night vibe)
# -------------------------
st.markdown(
    """
    <style>
    /* Global background */
    .stApp {
        background: radial-gradient(circle at 0% 0%, #1a1b27 0, #15161e 40%, #101018 100%);
        color: #e5e9f0;
    }

    .block-container {
        padding-top: 1.5rem;
    }

    [data-testid="stMetricValue"] {
        color: #c0caf5;
        font-weight: 700;
    }

    [data-testid="stMetricLabel"] {
        color: #9aa5ce;
    }

    .stDataFrame table {
        background-color: #1a1b26 !important;
        color: #e5e9f0 !important;
    }

    section[data-testid="stSidebar"] {
        background: linear-gradient(180deg, #1f2335 0%, #15161e 100%);
    }

    button[role="tab"] {
        border-radius: 999px !important;
        padding: 0.3rem 0.9rem !important;
        margin-right: 0.3rem !important;
        font-weight: 500 !important;
    }
    </style>
    """,
    unsafe_allow_html=True,
)


# -------------------------
# Helper: load local images
# -------------------------
@st.cache_data
def load_image(path: str) -> Optional[Image.Image]:
    try:
        return Image.open(path)
    except Exception:
        return None


# -------------------------
# Header with hero banner
# -------------------------
banner_cols = st.columns([2, 1])

with banner_cols[0]:
    st.title("🧬 Delly Trio Structural Variant Explorer")
    st.markdown(
        """
        <div style="font-size: 0.95rem; color: #c0caf5;">
        An interactive dashboard around for analysing <b>Delly VCF Trios</b>.<br/>
        Upload a multi-sample SV VCF, infer trio roles &amp; sex, and explore inheritance,
        SV type distributions, and putative de novo calls with a visual, talk-ready UI.
        </div>
        """,
        unsafe_allow_html=True,
    )

with banner_cols[1]:
    banner_img = load_image("../assets/trio_banner.png")
    if banner_img is not None:
        st.image(banner_img, use_container_width=True)
    else:
        st.markdown(
            """
            <div style="
                border-radius: 1rem;
                padding: 1rem;
                background: linear-gradient(135deg,#bb9af7,#7dcfff);
                color:#1a1b26;
                font-weight:600;
                text-align:center;">
                Place a banner at<br/><code>../assets/trio_banner.png</code><br/>to see it here.
            </div>
            """,
            unsafe_allow_html=True,
        )

st.markdown("---")

st.markdown(
    """
    <div style="font-size:0.95rem; color:#c0caf5; padding-bottom:0.4rem;">
    <b>📁 Upload a VCF and run the analysis from the sidebar to see results.</b>
    </div>
    """,
    unsafe_allow_html=True,
)

st.markdown("---")

# Pipeline workflow diagram
pipeline_img = load_image("../assets/Pipeline.png")
if pipeline_img is not None:
    st.image(pipeline_img, width=1000)
else:
    st.info("Pipeline.png not found in assets/. Place it at: assets/Pipeline.png")


# -------------------------
# Helpers (analysis logic)
# -------------------------
def save_uploaded_vcf(uploaded_file) -> str:
    """Save uploaded VCF to a temp directory and return the path."""
    if "tmpdir" not in st.session_state:
        st.session_state["tmpdir"] = tempfile.mkdtemp(prefix="trio_analyzer_")
    tmpdir = st.session_state["tmpdir"]

    vcf_path = os.path.join(tmpdir, uploaded_file.name)
    with open(vcf_path, "wb") as f:
        f.write(uploaded_file.getbuffer())
    return vcf_path


def run_first_awk_pass(
    vcf_path: str,
    out_dir: str,
    child: str,
    parent1: str,
    parent2: str,
) -> str:
    """
    Run the first AWK parsing pass (using vca.AWK_PARSE) to generate trio_table.tsv.
    Returns path to trio_table.tsv.
    """
    parsed_tsv = os.path.join(out_dir, "trio_table.tsv")

    with tempfile.NamedTemporaryFile("w", delete=False) as tf:
        tf.write(vca.AWK_PARSE)
        awkfile = tf.name

    with open(parsed_tsv, "w") as o:
        subprocess.run(
            [
                "awk",
                "-v",
                f"sample1={child}",
                "-v",
                f"sample2={parent1}",
                "-v",
                f"sample3={parent2}",
                "-f",
                awkfile,
                vcf_path,
            ],
            stdout=o,
            check=True,
        )

    return parsed_tsv


def compute_sample_to_col(vcf_path: str) -> Dict[str, int]:
    """
    Map sample name -> column index in the VCF (1-based),
    as required by run_awk_denovo.
    """
    header_samples: List[str] = []
    with open(vcf_path) as f:
        for line in f:
            if line.startswith("#CHROM"):
                header_samples = line.strip().split("\t")
                break

    sample_to_col: Dict[str, int] = {}
    for idx, name in enumerate(header_samples):
        if idx >= 9:
            sample_to_col[name] = idx + 1  # 1-based for awk

    return sample_to_col


def run_full_pipeline(
    vcf_path: str,
    child_label: str,
    parent1_label: str,
    parent2_label: str,
) -> Tuple[pd.DataFrame, Optional[pd.DataFrame], Dict, Dict, str]:
    """
    High-level wrapper that:
      - runs AWK first pass
      - infers trio roles + sex
      - annotates SV length, size category, inheritance
      - (optionally) annotates with bedtools
      - runs de novo AWK pass
    """
    if "tmpdir" not in st.session_state:
        st.session_state["tmpdir"] = tempfile.mkdtemp(prefix="trio_analyzer_")
    tmpdir = st.session_state["tmpdir"]

    out_dir = os.path.join(tmpdir, "out")
    os.makedirs(out_dir, exist_ok=True)

    # 1) AWK primary parsing
    parsed_tsv = run_first_awk_pass(
        vcf_path, out_dir, child_label, parent1_label, parent2_label
    )

    # 2) Load trio table
    df = pd.read_csv(parsed_tsv, sep="\t", dtype=str)

    # 3) Extract sample names from *_GT columns
    samples = [c.replace("_GT", "") for c in df.columns if c.endswith("_GT")]

    # 4) Trio role inference
    roles, trio_stats, min_vio = vca.infer_trio(df, samples)
    if roles is None:
        st.error(
            f"""❌ Trio role inference failed.
    
            The minimum Mendelian-violation count across all child/parent assignments
            is **{min_vio}**, which is greater than the allowed threshold (> 4).
            
            This suggests the three selected samples are unlikely to be a valid
            child–parent–parent trio (e.g. wrong samples, unrelated individuals, or bad genotypes).
            """
            )

        # Show the full permutation table so the user can see what was tried
        ts_df = pd.DataFrame(trio_stats)
        st.markdown("#### Tried trio assignments and their Mendelian violation counts")
        st.dataframe(ts_df, use_container_width=True)
    
        # Stop the rest of the analysis for this run
        st.stop()

    # If we get here, roles is valid and we can proceed
    child = roles["child"]
    p1 = roles["parent1"]
    p2 = roles["parent2"]

    # 5) Sex inference from chrX
    sex_calls, x_het_rates = vca.infer_sex_from_chrX(df, samples)

    # 6) SV length + size categorisation
    df = vca.add_sv_length_and_size_category(df)

    # 7) Inheritance annotation (requires trio roles)
    df = vca.annotate_inheritance(df, child=child, p1=p1, p2=p2)

    # 8) bedtools annotations (optional)
    try:
        df = vca.bedtools_annotate_all_variants(df, out_dir)
        # ---- Create clinically_prioritised_SVs.tsv in the SAME out_dir ----
        df_voi = vca.flag_variants_of_interest(df, child)
        voi_path = os.path.join(out_dir, "clinically_prioritised_SVs.tsv")
        df_voi.to_csv(voi_path, sep="\t", index=False)
        
        # Save for UI use
        st.session_state["df_voi"] = df_voi
        st.session_state["voi_path"] = voi_path

        bedtools_ok = True
    except Exception as e:
        bedtools_ok = False
        st.warning(
            f"bedtools annotation skipped (bedtools or BED files/BED paths problem): {e}"
        )
    
    # 9) De novo AWK pass
    df_denovo = None
    try:
        sample_to_col = compute_sample_to_col(vcf_path)
        parent1_col = sample_to_col[p1]
        parent2_col = sample_to_col[p2]
        child_col = sample_to_col[child]

        denovo_txt = os.path.join(out_dir, "de_novo_calls.txt")
        denovo_tsv = os.path.join(out_dir, "de_novo_calls.tsv")

        vca.run_awk_denovo(
            vcf_path,
            denovo_txt,
            denovo_tsv,
            parent1_col,
            parent2_col,
            child_col,
        )
        df_denovo = pd.read_csv(denovo_tsv, sep="\t", dtype=str)
    except Exception as e:
        st.warning(f"De novo AWK pass failed or skipped: {e}")

    # Cache into session_state for later tabs
    st.session_state["analysis_df"] = df
    st.session_state["denovo_df"] = df_denovo
    st.session_state["roles"] = roles
    st.session_state["sex_calls"] = sex_calls
    st.session_state["x_het_rates"] = x_het_rates
    st.session_state["out_dir"] = out_dir
    st.session_state["bedtools_ok"] = bedtools_ok
    st.session_state["trio_stats"] = trio_stats

    return df, df_denovo, roles, sex_calls, out_dir


# -------------------------
# Sidebar – Input controls
# -------------------------
st.sidebar.header("⚙️ Pipeline Controls")

sidebar_img = load_image("assets/sidebar_avatar.png")
if sidebar_img is not None:
    st.sidebar.image(sidebar_img, caption="Trio SV Analyzer", use_container_width=True)

st.sidebar.markdown(
    """
    1. Choose input source  
    2. If uploading: provide multi-sample SV VCF (hg38/GRCh38)  
    3. Pick child & parents  
    4. Click **Run Trio Analysis**
    """.strip()
)

# --- NEW: choose between sample VCF and upload ---
input_mode = st.sidebar.radio(
    "VCF input source",
    ["Use sample VCF", "Upload your own VCF"],
    index=0,
    key="vcf_input_mode",
)

sample_vcf_path = os.path.join(os.path.dirname(__file__), "data", "sample_trio_sv.vcf")
vcf_path = None
samples = []

if input_mode == "Use sample VCF":
    st.sidebar.markdown("Using bundled example: `data/sample_trio_sv.vcf`")
    if os.path.exists(sample_vcf_path):
        vcf_path = sample_vcf_path
        try:
            samples = vca.get_vcf_samples(vcf_path)
        except Exception as e:
            st.sidebar.error(f"Could not read samples from sample VCF: {e}")
    else:
        st.sidebar.error("Sample VCF not found at `data/sample_trio_sv.vcf`. Please create it.")
else:
    uploaded_vcf = st.sidebar.file_uploader(
        "Upload multi-sample SV VCF",
        type=["vcf", "vcf.gz"],
        key="user_vcf_uploader",
    )

    if uploaded_vcf is not None:
        vcf_path = save_uploaded_vcf(uploaded_vcf)
        try:
            samples = vca.get_vcf_samples(vcf_path)
        except Exception as e:
            st.sidebar.error(f"Could not read samples from uploaded VCF: {e}")

# --- Now proceed only if we have a VCF path and samples ---
if vcf_path is not None and samples:
    st.sidebar.markdown("**Detected samples in VCF**")
    st.sidebar.write(samples)

    child_sample = st.sidebar.selectbox("👤 Sample 1", options=samples, index=0, key="child_sample")
    parent1_options = [s for s in samples if s != child_sample]
    parent1_sample = st.sidebar.selectbox("👤 Sample 2", options=parent1_options, key="p1_sample")

    parent2_options = [s for s in samples if s not in (child_sample, parent1_sample)]
    parent2_sample = st.sidebar.selectbox(
        "👤 Sample 3",
        options=parent2_options,
        index=0 if parent2_options else 0,
        key="p2_sample",
    )

    run_btn = st.sidebar.button("🚀 Run Trio Analysis", key="run_trio_btn")

    if run_btn and parent2_options:
        with st.spinner("Running AWK + Python pipeline on selected VCF..."):
            df, df_denovo, roles, sex_calls, out_dir = run_full_pipeline(
                vcf_path,
                child_sample,
                parent1_sample,
                parent2_sample,
            )
    else:
        df = st.session_state.get("analysis_df")
        df_denovo = st.session_state.get("denovo_df")
else:
    df = st.session_state.get("analysis_df")
    df_denovo = st.session_state.get("denovo_df")


# -------------------------
# Main content / tabs
# -------------------------
if df is None:
    empty_cols = st.columns([2, 1])
    with empty_cols[0]:
        st.info("Awaiting VCF upload…")
    with empty_cols[1]:
        overview_img = load_image("../assets/overview_illustration.png")
        if overview_img is not None:
            st.image(overview_img, use_container_width=True)
else:
    roles = st.session_state.get("roles", {})
    sex_calls = st.session_state.get("sex_calls", {})
    x_het_rates = st.session_state.get("x_het_rates", {})
    trio_stats = st.session_state.get("trio_stats", [])
    out_dir = st.session_state.get("out_dir", "")
    bedtools_ok = st.session_state.get("bedtools_ok", False)

    child = roles.get("child")
    p1 = roles.get("parent1")
    p2 = roles.get("parent2")

    # Top summary metrics
    st.subheader("🔍 Trio Summary")
    cols = st.columns(4)
    cols[0].metric("Total SVs", f"{len(df):,}")
    cols[1].metric("Child", f"{child} ({sex_calls.get(child, 'unknown')})")
    cols[2].metric("Parent 1", f"{p1} ({sex_calls.get(p1, 'unknown')})")
    cols[3].metric("Parent 2", f"{p2} ({sex_calls.get(p2, 'unknown')})")

    with st.expander("🧬 Trio role search (Mendelian violations)", expanded=False):
        if trio_stats:
            trio_df = pd.DataFrame(trio_stats)
            st.dataframe(trio_df, use_container_width=True)
        else:
            st.write("No trio statistics found.")

    with st.expander("⚧ Sex inference (chrX heterozygosity)", expanded=False):
        rows = []
        for s in sex_calls:
            rows.append(
                {
                    "Sample": s,
                    "Sex": sex_calls[s],
                    "X_het_rate": x_het_rates.get(s),
                }
            )
        st.dataframe(pd.DataFrame(rows), use_container_width=True)

    # Tabs for deeper exploration
    tab_overview, tab_svtype, tab_inheritance, tab_denovo, tab_raw, tab_annot = st.tabs(
        [
            "📊 Overview",
            "🧩 SV Type Distribution",
            "👨‍👩‍👧 Inheritance Patterns",
            "🌟 De novo Variants",
            "📋 Raw Variant Table",
            "🧪 Annotation & Tiers",
        ]
    )

    # -------------------------
    # Overview tab
    # -------------------------
    with tab_overview:
        col_ov1, col_ov2 = st.columns([2, 1])
        with col_ov1:
            st.markdown("### 📊 Global Overview")

            if "SVTYPE" in df.columns:
                sv_counts = df["SVTYPE"].value_counts().reset_index()
                sv_counts.columns = ["SVTYPE", "COUNT"]

                inner_cols = st.columns(2)
                with inner_cols[0]:
                    st.markdown("**SVTYPE counts**")
                    st.dataframe(sv_counts, use_container_width=True)
                with inner_cols[1]:
                    fig = plot_count_bar(
                        sv_counts,
                        x="SVTYPE",
                        y="COUNT",
                        title="SVTYPE Distribution",
                        palette="set2",
                    )
                    st.plotly_chart(fig, use_container_width=True)

            else:
                st.warning("SVTYPE column not found in trio table.")

            if "SIZE_PRIORITY" in df.columns:
                size_counts = df["SIZE_PRIORITY"].value_counts().reset_index()
                size_counts.columns = ["SIZE_PRIORITY", "COUNT"]
                st.markdown("**Size priority distribution**")
                st.dataframe(size_counts, use_container_width=True)
                fig = plot_count_bar(
                    size_counts,
                    x="SIZE_PRIORITY",
                    y="COUNT",
                    title="Size Priority Distribution",
                    palette="pastel",
                )
                st.plotly_chart(fig, use_container_width=True)

            else:
                st.info("SIZE_PRIORITY not available (SV length annotation may have failed).")

            if bedtools_ok:
                st.markdown(
                    "✅ <b>bedtools</b>-based annotations were applied.",
                    unsafe_allow_html=True,
                )
            else:
                st.markdown(
                    "⚠️ <b>bedtools</b>-based annotations were skipped or failed.",
                    unsafe_allow_html=True,
                )

        with col_ov2:
            overview_img = load_image("../assets/overview_illustration.png")
            if overview_img is not None:
                st.image(
                    overview_img,
                    caption="Conceptual view of trio SV analysis.",
                    use_container_width=True,
                )

    # -------------------------
    # SV Type Distribution tab
    # -------------------------
    with tab_svtype:
        st.markdown("### 🧩 SV Type Distribution")
    
        if "SVTYPE" in df.columns and "CHROM" in df.columns:
            # --- Build ordered chromosome list ---
            raw_chroms = df["CHROM"].dropna().unique().tolist()
            chroms = sorted(raw_chroms, key=chrom_order_key)
    
            sv_types = sorted(df["SVTYPE"].dropna().unique().tolist())
    
            # Filters
            fc, fs = st.columns(2)
            with fc:
                selected_chroms = st.multiselect(
                    "Filter by Chromosome",
                    options=chroms,
                    default=chroms,
                    key="svtab_filter_chrom",
                )
            with fs:
                selected_svtypes = st.multiselect(
                    "Filter by SVTYPE",
                    options=sv_types,
                    default=sv_types,
                    key="svtab_filter_svtype",
                )
    
            # Apply filters
            filtered = df[
                df["CHROM"].isin(selected_chroms)
                & df["SVTYPE"].isin(selected_svtypes)
            ].copy()
    
            st.write(f"Showing **{len(filtered):,}** variants after filters.")
    
            # -------------------------
            # 1) Overall SVTYPE counts (after filters)
            # -------------------------
            sv_counts = (
                filtered["SVTYPE"]
                .value_counts()
                .reset_index()
            )
            sv_counts.columns = ["SVTYPE", "COUNT"]
    
            if not sv_counts.empty:
                st.markdown("#### 📊 SVTYPE distribution (filtered set)")
                fig1 = plot_count_bar(
                    sv_counts,
                    x="SVTYPE",
                    y="COUNT",
                    title="SVTYPE Counts (Filtered Variants)",
                    palette="set2",
                )
                st.plotly_chart(fig1, use_container_width=True)
            else:
                st.info("No variants after filtering to show SVTYPE counts.")
    
            # -------------------------
            # 2) Distribution of SV types across chromosomes
            # -------------------------
            st.markdown("#### 🧬 SVTYPE across chromosomes")
    
            if not filtered.empty:
                sv_by_chrom = (
                    filtered.groupby(["CHROM", "SVTYPE"])
                    .size()
                    .reset_index(name="COUNT")
                )
    
                if not sv_by_chrom.empty:
                    # enforce proper chromosome order for plotting
                    ordered_chroms = sorted(
                        sv_by_chrom["CHROM"].unique().tolist(),
                        key=chrom_order_key,
                    )
                    sv_by_chrom["CHROM"] = pd.Categorical(
                        sv_by_chrom["CHROM"],
                        categories=ordered_chroms,
                        ordered=True,
                    )
    
                    fig2 = px.bar(
                        sv_by_chrom,
                        x="CHROM",
                        y="COUNT",
                        color="SVTYPE",
                        text="COUNT",
                        barmode="group",
                        template="plotly_dark",
                        color_discrete_sequence=PALETTES.get("set3"),
                        title="Distribution of SV Types Across Chromosomes",
                    )
                    fig2.update_traces(
                        textposition="outside",
                        marker=dict(line=dict(width=1, color="rgba(255,255,255,0.3)")),
                        hovertemplate=(
                            "CHROM: %{x}<br>"
                            "SVTYPE: %{legendgroup}<br>"
                            "COUNT: %{y}<extra></extra>"
                        ),
                    )
                    fig2.update_layout(
                        font=dict(size=14),
                        margin=dict(t=60, b=80, l=40, r=20),
                        xaxis_title="Chromosome",
                        yaxis_title="Count",
                        height=500,
                        legend_title_text="SVTYPE",
                    )
                    st.plotly_chart(fig2, use_container_width=True)
    
                    st.markdown("##### Underlying SVTYPE × Chromosome table")
                    st.dataframe(sv_by_chrom, use_container_width=True)
                else:
                    st.info("No SVTYPE × chromosome combinations to plot.")
            else:
                st.info("No variants available after filtering to show chromosome distribution.")
    
            # -------------------------
            # 3) Filtered variant subset table (sorted by CHROM, then POS)
            # -------------------------
            st.markdown("##### Filtered variant subset")
    
            if not filtered.empty:
                # sort chromosomes and positions
                filtered = filtered.copy()
                filtered["CHROM"] = pd.Categorical(
                    filtered["CHROM"],
                    categories=sorted(
                        filtered["CHROM"].unique().tolist(),
                        key=chrom_order_key,
                    ),
                    ordered=True,
                )
    
                if "POS" in filtered.columns:
                    filtered["POS_INT"] = pd.to_numeric(filtered["POS"], errors="coerce")
                    sort_cols = ["CHROM", "POS_INT"]
                else:
                    sort_cols = ["CHROM"]
    
                subset_cols = [
                    c
                    for c in [
                        "CHROM",
                        "POS",
                        "SVTYPE",
                        "SVLEN_BP",
                        "SIZE_PRIORITY",
                        "INHERITANCE",
                    ]
                    if c in filtered.columns
                ]
    
                st.dataframe(
                    filtered.sort_values(sort_cols)[subset_cols],
                    use_container_width=True,
                )
            else:
                st.info("No variants to display in the filtered table.")
        else:
            st.warning("SVTYPE or CHROM column not found in the trio table.")

    # -------------------------
    # Inheritance Patterns tab
    # -------------------------
    with tab_inheritance:
        st.markdown("### 👨‍👩‍👧 Inheritance Patterns")
        if "INHERITANCE" in df.columns:
            inh_counts = df["INHERITANCE"].value_counts(dropna=False).reset_index()
            inh_counts.columns = ["INHERITANCE", "COUNT"]

            cols_inh = st.columns(2)
            with cols_inh[0]:
                st.dataframe(inh_counts, use_container_width=True)
            with cols_inh[1]:
                fig = plot_count_bar(
                    inh_counts,
                    x="INHERITANCE",
                    y="COUNT",
                    title="Inheritance Class Distribution",
                    palette="set3",
                )
                st.plotly_chart(fig, use_container_width=True)

            selected_inh = st.multiselect(
                "Filter variants by inheritance class",
                options=inh_counts["INHERITANCE"].tolist(),
                default=inh_counts["INHERITANCE"].tolist(),
            )
            df_inh = df[df["INHERITANCE"].isin(selected_inh)]

            subset_cols = [
                c
                for c in [
                    "CHROM",
                    "POS",
                    "SVTYPE",
                    "SVLEN_BP",
                    "SIZE_PRIORITY",
                    "INHERITANCE",
                ]
                if c in df_inh.columns
            ]
            st.write(f"Showing {len(df_inh):,} variants.")
            st.dataframe(df_inh[subset_cols], use_container_width=True)
        else:
            st.info("INHERITANCE column not found – inheritance annotation may have failed.")

    # -------------------------
    # De novo tab
    # -------------------------
    with tab_denovo:
        cols_dn = st.columns([2, 1])
        with cols_dn[0]:
            st.markdown("### 🌟 De novo Calls (AWK-based)")
            if df_denovo is not None:
                st.write(f"Total putative de novo calls: {len(df_denovo):,}")
                if "SV_TYPE" in df_denovo.columns:
                    svtypes_dn = sorted(df_denovo["SV_TYPE"].dropna().unique().tolist())
                    sel_dn = st.multiselect(
                        "Filter de novo by SV_TYPE",
                        options=svtypes_dn,
                        default=svtypes_dn,
                    )
                    df_dn = df_denovo[df_denovo["SV_TYPE"].isin(sel_dn)]
                else:
                    df_dn = df_denovo

                st.dataframe(df_dn, use_container_width=True)

                dn_tsv_path = os.path.join(out_dir, "de_novo_calls.tsv")
                if os.path.exists(dn_tsv_path):
                    with open(dn_tsv_path, "rb") as f:
                        st.download_button(
                            "Download de_novo_calls.tsv",
                            data=f,
                            file_name="de_novo_calls.tsv",
                            mime="text/tab-separated-values",
                        )
            else:
                st.info("No de novo dataframe available (AWK pass may have been skipped/failed).")
        with cols_dn[1]:
            dn_img = load_image("../assets/denovo_illustration.png")
            if dn_img is not None:
                st.image(
                    dn_img,
                    caption="Putative de novo SVs in the child.",
                    use_container_width=True,
                )

    # -------------------------
    # Raw variant table tab
    # -------------------------
    with tab_raw:
        st.markdown("### 📋 Full Annotated Variant Table")

        chroms = sorted(df["CHROM"].dropna().unique().tolist()) if "CHROM" in df.columns else []
        svtypes = sorted(df["SVTYPE"].dropna().unique().tolist()) if "SVTYPE" in df.columns else []

        left, right = st.columns(2)
        with left:
            sel_chrom = st.multiselect("Chromosome", options=chroms, default=chroms)
        with right:
            sel_svtype = st.multiselect("SVTYPE", options=svtypes, default=svtypes)

        df_view = df.copy()
        if sel_chrom:
            df_view = df_view[df_view["CHROM"].isin(sel_chrom)]
        if sel_svtype:
            df_view = df_view[df_view["SVTYPE"].isin(sel_svtype)]

        st.write(f"Showing {len(df_view):,} variants after filters.")
        st.dataframe(df_view, use_container_width=True)

        out_tsv = os.path.join(out_dir, "trio_table_annotated_streamlit.tsv")
        df_view.to_csv(out_tsv, sep="\t", index=False)
        with open(out_tsv, "rb") as f:
            st.download_button(
                "Download current view as TSV",
                data=f,
                file_name="trio_table_annotated_streamlit.tsv",
                mime="text/tab-separated-values",
            )

    # -------------------------
    # Annotation & Pathogenicity Tiers tab
    # -------------------------
    with tab_annot:
        st.markdown("### 🧪 Annotation & Pathogenicity Tiers")
    
        df_full = st.session_state.get("analysis_df")
        df_voi = st.session_state.get("df_voi")
        voi_path = st.session_state.get("voi_path")
        out_dir = st.session_state.get("out_dir")
    
        if df_full is None:
            st.info("Run the analysis first.")
            st.stop()
    
        # -------------------------------------------------
        # Filters (Chromosome + SVTYPE)
        # -------------------------------------------------
        st.markdown("#### 🔍 Filters")
        chroms = sorted(df_full["CHROM"].dropna().unique().tolist())
        svtypes = sorted(df_full["SVTYPE"].dropna().unique().tolist())
    
        fc, ft = st.columns(2)
        with fc:
            sel_chrom = st.multiselect(
                "Filter by Chromosome",
                chroms,
                chroms,
                key="annot_filter_chrom"
            )
        with ft:
            sel_sv = st.multiselect(
                "Filter by SVTYPE",
                svtypes,
                svtypes,
                key="annot_filter_svtype"
            )
    
        df_filtered = df_full[
            df_full["CHROM"].isin(sel_chrom) &
            df_full["SVTYPE"].isin(sel_sv)
        ].copy()
    
        # -------------------------------------------------
        # Annotation table
        # -------------------------------------------------
        st.markdown("### 📘 Full Annotation Table")
    
        annotation_cols = [
            "CHROM","POS","END_INT","SVTYPE","SVLEN_BP","SIZE_PRIORITY",
            "GENE_LIST","GENE_OVERLAP_COUNT",
            "EXON_OVERLAP_COUNT",
            "CLINGEN_HI_OVERLAP_COUNT",
            "CLINGEN_TS_OVERLAP_COUNT",
            "CLINGEN_REC_CNV_OVERLAP_COUNT",
            "CLINVAR_OVERLAP_COUNT",
            "CLINVAR_PHENOTYPE",
            "INHERITANCE",
            "PATHOGENICITY_TIER",
        ]
        annotation_cols = [c for c in annotation_cols if c in df_filtered.columns]
    
        st.dataframe(df_filtered[annotation_cols], use_container_width=True)
    
        # -------------------------------------------------
        # Pathogenicity tiers (Tier0/1/2)
        # -------------------------------------------------
        st.markdown("### ⭐ Clinically Prioritised SVs (Tier 0 / 1 / 2)")
    
        if df_voi is None or df_voi.empty:
            st.warning("No Tier0/1/2 variants produced for this VCF.")
        else:
            st.success(f"Found **{len(df_voi)}** clinically prioritised variants.")
    
            tier_counts = df_voi["PATHOGENICITY_TIER"].value_counts().reset_index()
            tier_counts.columns = ["PATHOGENICITY_TIER", "COUNT"]
    
            t1, t2 = st.columns(2)
            with t1:
                st.dataframe(tier_counts, use_container_width=True)
            with t2:
                fig = plot_count_bar(
                    tier_counts,
                    x="PATHOGENICITY_TIER",
                    y="COUNT",
                    title="Pathogenicity Tier Breakdown",
                    palette="bold",
                )
                st.plotly_chart(fig, use_container_width=True)
    
            # Tier filter
            selected_tiers = st.multiselect(
                "Filter by tier",
                tier_counts["PATHOGENICITY_TIER"].tolist(),
                tier_counts["PATHOGENICITY_TIER"].tolist(),
                key="annot_filter_tier"
            )

            df_tier_filtered = df_voi[df_voi["PATHOGENICITY_TIER"].isin(selected_tiers)]
    
            st.markdown(f"Showing **{len(df_tier_filtered)}** prioritised variants")
            st.dataframe(df_tier_filtered, use_container_width=True)
    
            # Download button
            if voi_path and os.path.exists(voi_path):
                with open(voi_path, "rb") as f:
                    st.download_button(
                        "⬇️ Download clinically_prioritised_SVs.tsv",
                        data=f,
                        file_name="clinically_prioritised_SVs.tsv",
                        mime="text/tab-separated-values",
                    )


# -------------------------
# Footer: logo + author credit
# -------------------------
st.markdown("---")

footer_cols = st.columns([1, 3])

with footer_cols[0]:
    footer_logo = load_image("../assets/logo.png")
    if footer_logo is not None:
        st.image(
            footer_logo,
            use_container_width=True,
        )
    else:
        st.markdown(
"""
<div style="
    border-radius: 0.75rem;
    padding: 0.6rem;
    background: linear-gradient(135deg,#7aa2f7,#bb9af7);
    color:#1a1b26;
    font-weight:600;
    text-align:center;
    font-size:0.8rem;">
    Add a logo at<br/><code>../assets/logo.png</code>
</div>
""",
            unsafe_allow_html=True,
        )

with footer_cols[1]:
    st.markdown(
"""
<style>
/* Footer link row layout */
.footer-links {
  display: flex;
  flex-direction: row;
  align-items: center;
  gap: 14px;
  margin-top: 8px;
}

/* Base link styling */
.footer-links a {
  text-decoration: none;
  color: #7aa2f7;
  display: flex;
  align-items: center;
  font-size: 0.85rem;
  transition: all 0.2s ease-out;
}

/* Icon styling */
.footer-links a svg {
  margin-right: 6px;
  fill: #7aa2f7;
  transition: all 0.2s ease-out;
}

/* Separator between links */
.footer-links span.sep {
  color: #565f89;
}

/* Hover animation + glow */
.footer-links a:hover {
  color: #c0caf5;
  text-shadow: 0 0 6px rgba(122, 162, 247, 0.8);
  transform: translateY(-1px);
}

/* Icon glow + lift on hover */
.footer-links a:hover svg {
  filter: drop-shadow(0 0 6px rgba(122, 162, 247, 0.9));
  transform: translateY(-1px);
}
</style>

<div style="font-size: 0.85rem; color: #9aa5ce; padding-top: 0.4rem;">
  <b>VCF Trio SV Analyzer</b><br/>
  Author: <span style="color:#c0caf5;">Ekarsi Lodh</span><br/>
  Built with <span style="color:#7dcfff;">Python</span> &amp;
  <span style="color:#7dcfff;">Streamlit</span> • Designed to be
  <i>presentation-ready</i> for talks, posters, and documentation.
</div>

<div style="
  margin-top: 0.8rem;
  font-size: 0.85rem;
  color:#c0caf5;
  border-top: 1px solid #3b4261;
  padding-top: 0.75rem;
">
  <span style="font-weight:600;">Links &amp; contact:</span>

  <div class="footer-links">

<a href="https://github.com/ekarsilodh" target="_blank">
  <svg xmlns="http://www.w3.org/2000/svg" width="20" height="20"
       viewBox="0 0 16 16">
    <path d="M8 0C3.58 0 0 3.58 0 8c0 3.54 2.29 6.53 5.47 
             7.59.4.07.55-.17.55-.38 0-.19-.01-.82-.01-1.49
             -2.01.37-2.53-.49-2.69-.94-.09-.23-.48-.94-.82-1.13
             -.28-.15-.68-.52-.01-.53.63-.01 1.08.58 1.23.82
             .72 1.21 1.87.87 2.33.66.07-.52.28-.87.51-1.07
             -1.78-.2-3.64-.89-3.64-3.95 0-.87.31-1.59.82-2.15
             -.08-.2-.36-1.02.08-2.12 0 0 .67-.21 2.2.82
             .64-.18 1.32-.27 2-.27s1.36.09 2 .27
             c1.53-1.04 2.2-.82 2.2-.82.44 1.1.16 1.92.08 2.12
             .51.56.82 1.27.82 2.15 0 3.07-1.87 3.75-3.65 3.95
             .29.25.54.73.54 1.48 0 1.07-.01 1.93-.01 2.19
             0 .21.15.46.55.38A8.013 8.013 0 0 0 16 8
             c0-4.42-3.58-8-8-8z"/>
  </svg>
  GitHub
</a>

<span class="sep">|</span>

<a href="https://www.linkedin.com/in/ekarsi24" target="_blank">
  <svg xmlns="http://www.w3.org/2000/svg" width="20" height="20"
       viewBox="0 0 16 16">
    <path d="M1.146 1.146C1.542.75 2.07.5 2.625.5s1.083.25 1.48.646c.396.396.646.924.646 1.48
             0 .555-.25 1.083-.646 1.48-.397.396-.925.646-1.48.646S1.542 4.002 1.146 3.606
             .5 2.57.5 2.125c0-.555.25-1.083.646-1.48zM.5 5.5h4.25v10H.5v-10zm5.25 0h4.077v1.36h.058
             c.568-1.076 1.955-2.21 4.023-2.21C15.71 4.65 16 7.01 16 9.979V15.5h-4.25V10.59
             c0-1.17-.021-2.67-1.63-2.67-1.63 0-1.88 1.27-1.88 2.58V15.5H5.75v-10z"/>
  </svg>
  LinkedIn
</a>

<span class="sep">|</span>

<a href="https://orcid.org/0009-0000-7462-3217" target="_blank">
  <svg xmlns="http://www.w3.org/2000/svg" width="20" height="20"
       viewBox="0 0 256 256">
    <circle cx="128" cy="128" r="120"/>
    <path d="M95 175h-16V81h16v94zm8-94h34c30 0 49 19 49 47 0 29-19 47-49 47h-34V81zm34 79
             c20 0 32-12 32-32 0-19-12-32-32-32h-18v64h18z" fill="#1a1b26"/>
  </svg>
  ORCID
</a>

<span class="sep">|</span>

<a href="https://mail.google.com/mail/?view=cm&fs=1&to=exl533@student.bham.ac.uk"
   target="_blank">
  <svg xmlns="http://www.w3.org/2000/svg" width="20" height="20"
       viewBox="0 0 16 16">
    <path d="M0 4a2 2 0 0 1 2-2h12a2 2 0 0 1 2 2v.5L8 9 0 4.5V4z"/>
    <path d="M0 5.697V12a2 2 0 0 0 2 2h12a2 2 0 0 0 2-2V5.697L8 10.25 0 5.697z"/>
  </svg>
  Email
</a>

</div>
</div>
""",
        unsafe_allow_html=True,
    )
