import io
import pandas as pd
import numpy as np
from itertools import combinations
from itertools import combinations

from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats
from pydeseq2.default_inference import DefaultInference

from modules.validators import safe_filename
from modules.gene_mapper import map_genes
from modules.plots import plot_volcano, plot_pca

async def run_deseq2_microservice(
    counts_df,
    meta_df,
    p_adj,
    log2fc,
    mapper_file,
    labels,
    processing_stats,
    design="~disease_state"
):
    # This dictionary will hold ALL the outputs in-memory to send back to the UI
    output_payload = {
        "deg_csvs": {},        # Will store filename -> CSV string content
        "mapped_csvs": {},     # Will store filename -> CSV string content
        "plots": {},           # Will store plotname -> base64 image strings
        "manifest": None,      # Will store manifest JSON dictionary
        "qc_summary": None,    # Will store the summary row dictionary
        "vst_counts_csv": None # Will store VST data as CSV text
    }

    inference = DefaultInference()
    dds = DeseqDataSet(
        counts=counts_df.T,
        metadata=meta_df,
        design=design,
        refit_cooks=True,
        inference=inference
    )

    print("Running DESeq2 normalization...")
    dds.deseq2()

    out = list(combinations(labels, 2))
    contrasts = []
    for c in out:
        if 'Control' in c:
            control_idx = c.index('Control')
            contrast = ['disease_state', c[1-control_idx], c[control_idx]]
            contrasts.append(contrast)

    if len(contrasts) == 0:
        raise ValueError("No Control group contrasts found.")

    total_deg_count = 0
    total_up = 0
    total_down = 0
    manifest_data = None

    for contrast in contrasts:
        print(f"Running contrast: {contrast[1]} vs {contrast[2]}")
        de_stats = DeseqStats(dds, contrast=contrast, inference=inference)
        de_stats.summary()
        results = de_stats.results_df.copy()

        if results.empty:
            print(f"WARNING: Empty results for {contrast}")
            continue

        sig = results[
            (results['padj'] < p_adj) &
            (abs(results['log2FoldChange']) > log2fc)
        ].copy()

        sig["direction"] = np.where(sig["log2FoldChange"] > 0, "UP", "DOWN")
        upregulated = (sig["direction"] == "UP").sum()
        downregulated = (sig["direction"] == "DOWN").sum()

        total_deg_count += sig.shape[0]
        total_up += upregulated
        total_down += downregulated

        # =================================================
        # IN-MEMORY REPLACEMENT FOR: sig.to_csv(sig_file)
        # =================================================
        csv_buffer = io.StringIO()
        sig.to_csv(csv_buffer)
        deg_filename = f"{safe_filename(contrast[1])}_vs_{safe_filename(contrast[2])}_DEGs.csv"
        # Store raw text content inside our payload map
        output_payload["deg_csvs"][deg_filename] = csv_buffer.getvalue()

        print(f"Significant genes: {sig.shape[0]}")

        mapped_results = await map_genes(results, mapper_file)

        # =================================================
        # MODIFIED UTILITIES
        # Modify your custom `save_results` tool to output in-memory data
        # =================================================
        res_buffer = io.StringIO()
        mapped_results.to_csv(res_buffer)
        mapped_filename = f"{safe_filename(contrast[1])}_vs_{safe_filename(contrast[2])}_MAPPED.csv"
        output_payload["mapped_csvs"][mapped_filename] = res_buffer.getvalue()

        # =================================================
        # PLOT REPLACEMENT
        # Modify your `plot_volcano` function to return base64 string
        # using buffer.getvalue() instead of plt.savefig(path)
        # =================================================
        plot_key = f"{safe_filename(contrast[1])}_vs_{safe_filename(contrast[2])}_volcano"
        output_payload["plots"][plot_key] = plot_volcano(mapped_results, contrast)

        # =================================================
        # MANIFEST IN-MEMORY
        # =================================================
        manifest_data = {
            "disease": contrast[1],
            "n_cases": int((meta_df["disease_state"] == contrast[1]).sum()),
            "n_controls": int((meta_df["disease_state"] == "Control").sum()),
            "significant_genes": int(sig.shape[0]),
            "upregulated": int(upregulated),
            "downregulated": int(downregulated),
            "padj_threshold": p_adj,
            "logfc_threshold": log2fc,
            "mean_counts": float(counts_df.values.mean()),
            "median_counts": float(np.median(counts_df.values)),
            "mean_library_size": float(counts_df.sum(axis=0).mean())
        }
    
    output_payload["manifest"] = manifest_data

    # =====================================================
    # VST IN-MEMORY
    # =====================================================
    print("Generating VST counts...")
    dds.vst(use_design=False)
    vst_counts = dds.layers["vst_counts"]

    vst_df = pd.DataFrame(vst_counts, index=meta_df.index, columns=counts_df.index)
    
    vst_buffer = io.StringIO()
    vst_df.to_csv(vst_buffer)
    output_payload["vst_counts_csv"] = vst_buffer.getvalue()

    # =====================================================
    # PCA IN-MEMORY
    # Ensure your custom `plot_pca` tool returns (metrics, plot_base64)
    # =====================================================
    pca_metrics, pca_plot_base64 = plot_pca(vst_counts, meta_df)
    output_payload["plots"]["pca_plot"] = pca_plot_base64

    # =====================================================
    # QC SUMMARY IN-MEMORY
    # =====================================================
    qc_summary = {
    # 1. Cast incoming stats from the loader file
    "samples_removed": int(processing_stats["samples_removed"]),
    "samples_retained": int(processing_stats["samples_retained"]),
    "genes_filtered": int(processing_stats["genes_filtered"]),
    
    # 2. Cast loop accumulation counters (total_deg_count etc.)
    "total_DEGs": int(total_deg_count),
    "upregulated": int(total_up),
    "downregulated": int(total_down),
    
    # 3. Explicitly read standard floats and strings from PCA block
    "variance_explained_PC1": float(pca_metrics["pc1_variance"]),
    "variance_explained_PC2": float(pca_metrics["pc2_variance"]),
    "cluster_separation": float(pca_metrics["cluster_separation"]) if pca_metrics["cluster_separation"] is not None else None,
    "batch_effect": str(pca_metrics["batch_effect"]),
    "significant_DEGs": int(total_deg_count)
    }
    
    output_payload["qc_summary"] = qc_summary

    print("DESeq2 pipeline completed successfully.")
    return output_payload # Hand the entire package directly over to FastAPI
