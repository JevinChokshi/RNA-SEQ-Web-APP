import io
import os
import pandas as pd
import numpy as np
import base64
import matplotlib.pyplot as plt
import seaborn as sns
from adjustText import adjust_text
from sklearn.decomposition import PCA
from sklearn.metrics import silhouette_score
from matplotlib.patches import Ellipse
from scipy.spatial.distance import pdist

import matplotlib
matplotlib.use('Agg')

def plot_volcano(
    results_mapped: pd.DataFrame,
    contrast: list,
    padj_threshold: float = 0.05,
    log2fc_threshold: float = 1.5,
    top_n_labels: int = 8
) -> None:

    """
    Publication-grade volcano plot.
    """

    # ======================================
    # CLEAN DATA
    # ======================================

    df = results_mapped.copy()

    df = df.replace([np.inf, -np.inf], np.nan)

    df = df.dropna(
        subset=[
            "padj",
            "log2FoldChange"
        ]
    )

    df["padj"] = df["padj"].clip(lower=1e-300)

    df["neglog10_padj"] = -np.log10(df["padj"])

    # ======================================
    # SIGNIFICANCE GROUPS
    # ======================================

    df["category"] = "Not Significant"

    up_mask = (
        (df["padj"] < padj_threshold)
        &
        (df["log2FoldChange"] >= log2fc_threshold)
    )

    down_mask = (
        (df["padj"] < padj_threshold)
        &
        (df["log2FoldChange"] <= -log2fc_threshold)
    )

    df.loc[up_mask, "category"] = "Upregulated"

    df.loc[down_mask, "category"] = "Downregulated"

    # ======================================
    # COLORS
    # ======================================

    colors = {

        "Not Significant": "#BDBDBD",

        "Upregulated": "#D62728",

        "Downregulated": "#1F77B4"
    }

    # ======================================
    # FIGURE
    # ======================================

    plt.figure(figsize=(10, 8))

    ax = plt.gca()

    # ======================================
    # PLOT POINTS
    # ======================================

    for category in [

        "Not Significant",

        "Upregulated",

        "Downregulated"
    ]:

        subset = df[
            df["category"] == category
        ]

        ax.scatter(

            subset["log2FoldChange"],

            subset["neglog10_padj"],

            c=colors[category],

            s=28,

            alpha=0.75,

            edgecolors="none",

            label=category,

            rasterized=True
        )

    # ======================================
    # THRESHOLD LINES
    # ======================================

    ax.axvline(
        x=log2fc_threshold,
        linestyle="--",
        linewidth=1.2,
        color="black",
        alpha=0.7
    )

    ax.axvline(
        x=-log2fc_threshold,
        linestyle="--",
        linewidth=1.2,
        color="black",
        alpha=0.7
    )

    ax.axhline(
        y=-np.log10(padj_threshold),
        linestyle="--",
        linewidth=1.2,
        color="black",
        alpha=0.7
    )

    # ======================================
    # LABEL TOP GENES ONLY
    # ======================================

    label_df = df[
        df["category"] != "Not Significant"
    ].copy()

    label_df = label_df.sort_values(
        by="neglog10_padj",
        ascending=False
    )

    label_df = label_df.head(top_n_labels)

    texts = []

    for _, row in label_df.iterrows():

        gene = row.get("Symbol")

        if pd.isna(gene):

            gene = row.get("GeneID")

        txt = ax.text(

            row["log2FoldChange"],

            row["neglog10_padj"],

            str(gene),

            fontsize=11,

            fontweight="bold"
        )

        texts.append(txt)

    adjust_text(

        texts,

        arrowprops=dict(

            arrowstyle="-",

            color="black",

            lw=0.8
        )
    )

    # ======================================
    # STYLING
    # ======================================

    ax.set_xlabel(

        r"$\log_2$ Fold Change",

        fontsize=15,

        fontweight="bold"
    )

    ax.set_ylabel(

        r"$-\log_{10}$(Adjusted P-value)",

        fontsize=15,

        fontweight="bold"
    )

    ax.set_title(

        f"{contrast[1]} vs {contrast[2]}",

        fontsize=17,

        fontweight="bold"
    )

    ax.tick_params(

        axis="both",

        labelsize=12
    )

    # ======================================
    # REMOVE TOP/RIGHT SPINES
    # ======================================

    ax.spines["top"].set_visible(False)

    ax.spines["right"].set_visible(False)


    # ======================================
    # LEGEND
    # ======================================

    legend = ax.legend(

        frameon=False,

        fontsize=11,

        loc="center left",

        bbox_to_anchor=(1.02, 0.5),

        borderaxespad=0
    )

    # ======================================
    # ADD RIGHT MARGIN FOR LEGEND
    # ======================================

    plt.subplots_adjust(right=0.82)

    # ======================================
    # LAYOUT
    # ======================================

    plt.tight_layout()

    # ======================================
    # SAVE
    # ======================================

    img_buffer = io.BytesIO()
    plt.savefig(img_buffer, format='png', bbox_inches='tight')
    plt.close()
    img_buffer.seek(0)

    img_str = base64.b64encode(img_buffer.read()).decode('utf-8')
    return img_str



def detect_pca_outliers(df_pc, z_threshold=3):

    outliers = []

    for axis in ["PC1", "PC2"]:

        z_scores = (
            (df_pc[axis] - df_pc[axis].mean())
            / df_pc[axis].std()
        )

        axis_outliers = df_pc[
            np.abs(z_scores) > z_threshold
        ]["sample"].tolist()

        outliers.extend(axis_outliers)

    return list(set(outliers))


def estimate_batch_effect(df_pc):

    distances = pdist(
        df_pc[["PC1", "PC2"]].values
    )

    mean_distance = np.mean(distances)

    if mean_distance < 2:
        return "HIGH"

    elif mean_distance < 5:
        return "MODERATE"

    return "LOW"


def confidence_ellipse(
    x,
    y,
    ax,
    n_std=2.0,
    facecolor='none',
    **kwargs
):

    if x.size != y.size:

        raise ValueError(
            "x and y must be same size"
        )

    cov = np.cov(x, y)

    vals, vecs = np.linalg.eigh(cov)

    order = vals.argsort()[::-1]

    vals = vals[order]

    vecs = vecs[:, order]

    theta = np.degrees(

        np.arctan2(
            *vecs[:, 0][::-1]
        )
    )

    width, height = 2 * n_std * np.sqrt(vals)

    ellipse = Ellipse(

        xy=(np.mean(x), np.mean(y)),

        width=width,

        height=height,

        angle=theta,

        facecolor=facecolor,

        **kwargs
    )

    ax.add_patch(ellipse)

    return ellipse


def plot_pca(
    vst_counts,
    meta_df,
):

    pca = PCA(n_components=2)

    pc = pca.fit_transform(vst_counts)

    explained = (
        pca.explained_variance_ratio_ * 100
    )

    df_pc = pd.DataFrame({

        "sample": meta_df.index,

        "PC1": pc[:, 0],

        "PC2": pc[:, 1],

        "condition":
            meta_df[
                "disease_state"
            ].values
    })


    # =====================================================
    # METRICS
    # =====================================================

    outliers = detect_pca_outliers(df_pc)

    cluster_sep = None

    if len(df_pc["condition"].unique()) > 1:

        try:

            cluster_sep = silhouette_score(
                df_pc[["PC1", "PC2"]],
                df_pc["condition"]
            )

        except:
            cluster_sep = None

    batch_effect = estimate_batch_effect(df_pc)

    pca_metrics = {
        "pc1_variance": round(float(explained[0]), 3),
        "pc2_variance": round(float(explained[1]), 3),
        "cluster_separation": round(float(cluster_sep), 3),
        "outlier_samples": outliers,
        "batch_effect": batch_effect
    }

    # =====================================================
    # COLORS
    # =====================================================

    conditions = sorted(
        df_pc["condition"].unique()
    )

    palette = sns.color_palette(
        "Set2",
        n_colors=len(conditions)
    )

    color_map = dict(
        zip(conditions, palette)
    )

    plt.figure(figsize=(9, 8))

    ax = plt.gca()

    for cond in conditions:

        subset = df_pc[
            df_pc["condition"] == cond
        ]

        ax.scatter(

            subset["PC1"],

            subset["PC2"],

            s=120,

            alpha=0.9,

            color=color_map[cond],

            edgecolor="black",

            linewidth=0.8,

            label=cond
        )

        if subset.shape[0] >= 3:

            confidence_ellipse(

                subset["PC1"].values,

                subset["PC2"].values,

                ax,

                edgecolor=color_map[cond],

                linestyle="--",

                linewidth=1.5,

                alpha=0.7
            )

    ax.set_xlabel(
        f"PC1 ({explained[0]:.1f}%)",
        fontsize=15,
        fontweight="bold"
    )

    ax.set_ylabel(
        f"PC2 ({explained[1]:.1f}%)",
        fontsize=15,
        fontweight="bold"
    )

    ax.set_title(
        "Principal Component Analysis",
        fontsize=18,
        fontweight="bold"
    )

    ax.grid(
        alpha=0.2,
        linestyle="--"
    )

    ax.spines["top"].set_visible(False)

    ax.spines["right"].set_visible(False)

    ax.tick_params(
        axis="both",
        labelsize=12
    )

    ax.legend(
        frameon=False,
        fontsize=11,
        loc="best"
    )

    
    buf = io.BytesIO()
    plt.savefig(buf, format='png', bbox_inches='tight')
    plt.close()
    buf.seek(0)
    
    # Convert binary image data into a clean text string for JSON transfer
    img_str = base64.b64encode(buf.read()).decode('utf-8')

    print(
        "Publication-grade PCA plot saved"
    )

    return pca_metrics, img_str

