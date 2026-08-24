import streamlit as st
import base64
import pandas as pd
from io import StringIO, BytesIO

from utils.api import deseq2_api


def render_deseq2_ui():
    # --- Layout Definition: Left Panel (Inputs) vs Right Panel (Results) ---
    col_sidebar, col_main = st.columns([1, 2])

    with col_sidebar:
        st.header("⚙️ Pipeline Configuration")
        
        # 1. File Uploads
        st.subheader("📁 Dataset Input Uploads")
        count_file = st.file_uploader(
            "1. Raw Counts Matrix File (TSV)",
            type=["tsv"]
        )

        meta_file = st.file_uploader(
            "2. Sample Metadata Table (CSV)",
            type=["csv"]
        )
        
        st.divider()
        
        # 2. Threshold Adjustments
        st.subheader("🎯 Statistical Cutoffs")

        p_adj = st.number_input(
            "Adjusted p-value (FDR Threshold)",
            min_value=0.000,
            max_value=1.000,
            value=0.05,
            step=0.01,
            format="%.3f"
        )

        log2fc = st.number_input(
            "Log2 Fold Change Threshold",
            min_value=0.0,
            max_value=10.0,
            value=1.0,
            step=0.1
        )
        
        st.divider()
        
        # 3. Microservice Matrix Mappings
        st.subheader("🔬 Metadata & Design Layout")

        design = st.text_input(
            "Experimental Design Formula",
            value="~condition"
        )

        min_gene_count = st.number_input(
            "Minimum Pre-filtering Gene Count",
            min_value=0,
            value=10,
            step=1
        )

        min_samples = st.number_input(
            "Minimum Required Samples",
            min_value=1,
            value=3,
            step=1
        )
        
        meta_index = st.text_input(
            "Metadata Sample ID Column Name",
            value="sample_id"
        )

        count_index = st.text_input(
            "Counts Matrix Gene ID Column Name",
            value="gene_id"
        )

        label_col = st.text_input(
            "Target Experimental Group Column",
            value="condition"
        )

        control_label = st.text_input(
            "Baseline/Control Group Reference Label",
            value="control"
        )
        
        st.divider()
        
        # 4. Action Button Execution Hook
        run_pipeline = st.button(
            "🚀 Execute DESeq2 Pipeline",
            type="primary",
            width='stretch'
        )


    with col_main:
        st.header("📊 Interactive Analysis Terminal")
        
        if run_pipeline:

            # Client-side validation
            if not (count_file and meta_file):

                st.error(
                    "❌ Action Blocked: You must upload all 2 tracking files "
                    "(Counts, Metadata) before running."
                )

            else:

                with st.spinner(
                    "Processing RNA-Seq matrices server-side... Please wait."
                ):

                    try:

                        # =====================================================
                        # GRADIO CLIENT API CALL
                        # =====================================================

                        data = deseq2_api(
                            count_file,
                            meta_file,
                            p_adj,
                            log2fc,
                            design,
                            min_gene_count,
                            min_samples,
                            meta_index,
                            count_index,
                            label_col,
                            control_label
                        )

                        st.success(
                            "✅ DESeq2 Pipeline completed successfully!"
                        )

                        # =====================================================
                        # RESULTS
                        # =====================================================

                        # --- Tab Element Setup for Results Segmentation ---
                        tab_summary, tab_plots, tab_tables = st.tabs(
                            [
                                "📋 Run Summary",
                                "🖼️ Graphics & Plots",
                                "📑 Tabular Data Exports"
                            ]
                        )

                        # =====================================================
                        # TAB 1 — SUMMARY
                        # =====================================================

                        with tab_summary:

                            if data.get("qc_summary"):

                                st.subheader(
                                    "Quality Control Metrics Overview"
                                )

                                st.json(
                                    data["qc_summary"]
                                )

                            if data.get("manifest"):

                                st.subheader(
                                    "Run Configuration History Manifest"
                                )

                                st.json(
                                    data["manifest"]
                                )

                        # =====================================================
                        # TAB 2 — PLOTS
                        # =====================================================

                        with tab_plots:

                            plots_dict = data.get("plots", {})

                            if plots_dict:

                                for plot_name, base64_str in plots_dict.items():

                                    st.subheader(
                                        f"📊 Plot: {plot_name}"
                                    )

                                    try:

                                        image_bytes = base64.b64decode(
                                            base64_str
                                        )

                                        st.image(
                                            image_bytes,
                                            width='stretch'
                                        )

                                        st.download_button(
                                            label=f"💾 Download {plot_name}",
                                            data=image_bytes,
                                            file_name=f"{plot_name}.png",
                                            mime="image/png",
                                            key=f"dl_{plot_name}"
                                        )

                                    except Exception as img_err:

                                        st.error(
                                            f"Could not render graphic "
                                            f"'{plot_name}': {str(img_err)}"
                                        )

                            else:

                                st.info(
                                    "No plot assets were returned from "
                                    "the microservice payload execution."
                                )

                        # =====================================================
                        # TAB 3 — TABLES
                        # =====================================================

                        with tab_tables:

                            # 3a. Differentially Expressed Genes
                            deg_dict = data.get(
                                "deg_csvs",
                                {}
                            )

                            if deg_dict:

                                st.subheader(
                                    "Differential Expression Statistics Profiles"
                                )

                                for fname, csv_content in deg_dict.items():

                                    df = pd.read_csv(
                                        StringIO(csv_content)
                                    )

                                    st.markdown(
                                        f"**Data Preview:** `{fname}` "
                                        f"({len(df)} records)"
                                    )

                                    st.dataframe(
                                        df.head(10),
                                        width='stretch'
                                    )

                                    st.download_button(
                                        label=f"📥 Download {fname}",
                                        data=csv_content,
                                        file_name=fname,
                                        mime="text/csv",
                                        key=f"dl_{fname}"
                                    )

                            # 3b. VST Counts
                            vst_content = data.get(
                                "vst_counts_csv"
                            )

                            if vst_content:

                                st.divider()

                                st.subheader(
                                    "Variance Stabilized Transformation "
                                    "(VST) Normalized Counts"
                                )

                                df_vst = pd.read_csv(
                                    StringIO(vst_content)
                                )

                                st.dataframe(
                                    df_vst.head(10),
                                    width='stretch'
                                )

                                st.download_button(
                                    label="📥 Download VST Normalized Counts Matrix",
                                    data=vst_content,
                                    file_name="vst_normalized_counts.csv",
                                    mime="text/csv"
                                )

                    except Exception as connection_error:

                        st.error(
                            "🚨 DESeq2 API call failed. "
                            "Please verify that the Hugging Face Space "
                            "is running."
                        )

                        st.exception(
                            connection_error
                        )

        else:

            st.info(
                "💡 Dashboard Idle: Tweak your parameters in the left "
                "config tray and hit 'Execute DESeq2 Pipeline' to begin "
                "processing."
            )