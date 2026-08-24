# DESeq2 RNA-Seq Analysis Platform

A web-based RNA-seq differential gene expression analysis platform built on **PyDESeq2**. The application provides an end-to-end workflow for running DESeq2 analysis directly from uploaded raw count data and sample metadata.

The platform automatically performs:

* Differential gene expression analysis using **PyDESeq2**
* Gene filtering and preprocessing
* **Ensembl Gene ID → Gene Symbol** mapping
* DEG identification
* **Volcano plot** generation
* **PCA plot** generation
* VST count generation
* Quality-control metric calculation
* Downloadable analysis results and visualizations

## Live Application

### Frontend

**Streamlit Application:**
https://deseq2-by-jevin.streamlit.app/

### API

**Hugging Face Space API:**
https://huggingface.co/spaces/Jevin23/DeSEQ2_API

---

# Features

* 🧬 PyDESeq2-based differential expression analysis
* 📊 Automatic DEG identification
* 🗺️ Automatic Ensembl Gene ID to Gene Symbol mapping
* 🌋 Interactive and downloadable Volcano plots
* 📈 PCA visualization and sample clustering
* 🧪 Variance Stabilizing Transformation (VST)
* 📥 Downloadable DEG tables
* 📥 Downloadable VST counts
* 📥 Downloadable PCA and Volcano plots
* ⚙️ Configurable statistical and preprocessing parameters
* 📋 Automatic quality-control summary
* 📜 Run configuration history manifest
* 🚀 Web API deployment through Hugging Face Spaces
* 🖥️ User-friendly Streamlit frontend

---

# How to Use the Web Application

## Step 1: Open the Application

Go to:

https://deseq2-by-jevin.streamlit.app/

---

## Step 2: Upload Input Files

Upload:

1. **Raw Counts File**
2. **Run Table / Sample Metadata File**

The count matrix should contain:

* Genes in rows
* Samples in columns
* Raw integer counts

The run table should contain:

* Sample identifiers
* Experimental condition or disease status
* Any metadata required for the DESeq2 design

The sample identifiers in the run table must match the sample column names in the count matrix.

---

## Step 3: Configure Analysis Parameters

Fill in the following parameters:

| Parameter                        | Description                                                | Example             |
| -------------------------------- | ---------------------------------------------------------- | ------------------- |
| Adjusted p-value Threshold       | Significance threshold                                     | `0.05`              |
| Log2 Fold Change Threshold       | Minimum absolute log2 fold change                          | `1.5`               |
| Design Formula                   | DESeq2 experimental design                                 | `~disease_state`    |
| Minimum Pre-filtering Gene Count | Minimum count threshold for gene filtering                 | `10`                |
| Minimum Required Samples         | Number of samples required to pass the filtering threshold | `2`                 |
| Run Table Sample ID Column       | Sample identifier column in metadata                       | `sample_id`         |
| Count Matrix Gene ID Column      | Gene identifier column                                     | `gene_id`           |
| Target Column                    | Experimental condition column                              | `disease_state`     |
| Control Value                    | Reference/control group                                    | `Healthy Volunteer` |

### Example Configuration

```text
Adjusted p-value threshold: 0.05
Log2FC threshold: 1.5
Design formula: ~disease_state
Minimum pre-filtering gene count: 10
Minimum required samples: 2
Run table sample ID column: sample
Count matrix Gene ID column: gene_id
Target column: disease_state
Control value: HIV non-infected
```

> **Important:** The Count Matrix Gene ID column should not be left blank.

---

## Step 4: Run the Pipeline

Click:

```text
Execute DESeq2 Pipeline
```

Wait a few minutes for the analysis to complete.

The exact runtime depends on:

* Number of genes
* Number of samples
* Size of the uploaded files
* Available server resources

---

# Results

The platform generates multiple result types.

## 1. Run Summary

The application generates a summary of the completed analysis, including quality-control metrics and run configuration.

The summary can be copied with a single click.

### Quality Control Metrics Overview

Example:

```json
{
    "samples_removed": 1,
    "samples_retained": 14,
    "genes_filtered": 22886,
    "total_DEGs": 5454,
    "upregulated": 2369,
    "downregulated": 3085,
    "variance_explained_PC1": 55.598,
    "variance_explained_PC2": 15.499,
    "cluster_separation": 0.23,
    "batch_effect": "LOW",
    "significant_DEGs": 5454
}
```

### Metric Descriptions

* **samples_removed** — Samples removed during preprocessing
* **samples_retained** — Samples retained for analysis
* **genes_filtered** — Genes remaining after filtering
* **total_DEGs** — Total significant differentially expressed genes
* **upregulated** — Significantly upregulated genes
* **downregulated** — Significantly downregulated genes
* **variance_explained_PC1** — Variance explained by the first principal component
* **variance_explained_PC2** — Variance explained by the second principal component
* **cluster_separation** — Measure of separation between experimental groups
* **batch_effect** — Estimated batch-effect assessment
* **significant_DEGs** — Total genes passing the configured statistical thresholds

---

## 2. Run Configuration History Manifest

Example:

```json
{
    "disease": "Type 2 Diabetes",
    "n_cases": 5,
    "n_controls": 5,
    "significant_genes": 2729,
    "upregulated": 1145,
    "downregulated": 1584,
    "padj_threshold": 0.05,
    "logfc_threshold": 1.5,
    "mean_counts": 1260.549813739929,
    "median_counts": 75,
    "mean_library_size": 20786466.42857143
}
```

This manifest summarizes:

* Disease or experimental condition
* Number of cases
* Number of controls
* Significant genes
* Upregulated genes
* Downregulated genes
* Statistical thresholds used
* Mean and median counts
* Mean sequencing library size

---

## 3. Volcano Plot

The platform generates a Volcano plot showing:

* Upregulated genes
* Downregulated genes
* Statistical significance
* Log2 fold changes

The plot is downloadable from the web application.

---

## 4. PCA Plot

Principal Component Analysis is performed to visualize:

* Sample clustering
* Group separation
* Experimental variation
* Major sources of variance

The PCA plot is downloadable.

---

## 5. Differentially Expressed Genes

The DEG results can be downloaded as CSV files.

Results include:

* Gene identifiers
* Gene symbols
* Log2 fold changes
* Statistical significance
* Adjusted p-values
* Expression statistics

The application automatically maps **Ensembl Gene IDs to gene symbols** where possible.

---

## 6. VST Counts

Variance Stabilizing Transformation counts are generated for downstream analysis and visualization.

The VST count matrix is downloadable as CSV.

---

# Running the Backend Locally

## 1. Clone the Repository

```bash
git clone https://github.com/YOUR_USERNAME/YOUR_REPOSITORY.git
```

Move into the project directory:

```bash
cd YOUR_REPOSITORY
```

---

## 2. Navigate to the Server

```bash
cd server
```

---

## 3. Install Dependencies

```bash
pip install -r requirements.txt
```

---

## 4. Start the FastAPI Server

```bash
uvicorn main:app
```

The backend should be available at:

```text
http://127.0.0.1:8000
```

---

# Running the API

Send a POST request to:

```text
http://127.0.0.1:8000/run_deseq2/
```

The API returns an output payload structured as:

```python
output_payload = {
    "deg_csvs": {},        # Filename -> CSV string content
    "mapped_csvs": {},     # Filename -> mapped CSV string content
    "plots": {},           # Plot name -> base64 image string
    "manifest": None,      # Run manifest JSON dictionary
    "qc_summary": None,    # QC summary dictionary
    "vst_counts_csv": None # VST data as CSV text
}
```

---

# Using the Hosted Hugging Face API Directly

You can use the hosted API without cloning the repository.

## 1. Install the Gradio Client

```bash
pip install gradio_client
```

---

## 2. Connect to the API

```python
from gradio_client import Client, handle_file

client = Client("Jevin23/DeSEQ2_API")
```

---

## 3. Run DESeq2 Analysis

```python
from gradio_client import Client, handle_file

client = Client("Jevin23/DeSEQ2_API")

result = client.predict(
    count_file=handle_file("path/to/raw_counts.csv"),
    meta_file=handle_file("path/to/run_table.csv"),
    p_adj=0.05,
    log2fc=1.5,
    design="~disease_state",
    min_gene_count=10,
    min_samples=2,
    meta_index="sample",
    count_index="gene",
    label_col="condition",
    control_label="Control",
    api_name="/run_deseq2",
)

print(result)
```

Replace:

```text
path/to/raw_counts.csv
```

with your raw count matrix.

Replace:

```text
path/to/run_table.csv
```

with your sample metadata/run table.

---

# API Parameters

| Parameter        | Description                           |
| ---------------- | ------------------------------------- |
| `count_file`     | Raw gene count matrix                 |
| `meta_file`      | Sample metadata/run table             |
| `p_adj`          | Adjusted p-value threshold            |
| `log2fc`         | Log2 fold-change threshold            |
| `design`         | DESeq2 design formula                 |
| `min_gene_count` | Minimum gene count for filtering      |
| `min_samples`    | Minimum required samples              |
| `meta_index`     | Sample ID column in the metadata file |
| `count_index`    | Gene ID column in the count matrix    |
| `label_col`      | Target/condition column               |
| `control_label`  | Control/reference group               |
| `api_name`       | Gradio API endpoint                   |

---

# Technical Architecture

## Backend

The backend is built using:

* **PyDESeq2** for differential expression analysis
* **FastAPI** for backend API functionality
* **Gradio** for the hosted API interface
* **Hugging Face Spaces** for backend deployment

The Hugging Face Space is deployed using the **Gradio SDK**.

## Frontend

The frontend is built using:

* **Streamlit**
* Deployed on **Streamlit Community Cloud**

The Streamlit application communicates with the backend to execute the DESeq2 pipeline and display results.

---

# Deployment Architecture

```text
                    ┌───────────────────────┐
                    │ Streamlit Community   │
                    │        Cloud          │
                    └───────────┬───────────┘
                                │
                                │ API Request
                                ▼
                    ┌───────────────────────┐
                    │  Hugging Face Spaces  │
                    │                       │
                    │ Gradio / FastAPI API  │
                    └───────────┬───────────┘
                                │
                                ▼
                    ┌───────────────────────┐
                    │       PyDESeq2        │
                    │                       │
                    │  DEG Analysis Engine  │
                    └───────────┬───────────┘
                                │
              ┌─────────────────┼─────────────────┐
              │                 │                 │
              ▼                 ▼                 ▼
          DEG Results       PCA Plot        Volcano Plot
              │                 │                 │
              └─────────────────┼─────────────────┘
                                │
                                ▼
                       Streamlit Results
```

---

# Hugging Face Space Configuration

The backend Hugging Face Space is configured using:

```yaml
---
title: Rna Seq Api
emoji: 🧬
colorFrom: blue
colorTo: indigo
sdk: gradio
sdk_version: 6.25.0
app_file: app.py
python_version: "3.10"
pinned: false
---
```

---

# Technology Stack

| Component               | Technology                      |
| ----------------------- | ------------------------------- |
| Differential Expression | PyDESeq2                        |
| Backend API             | FastAPI                         |
| API Interface           | Gradio                          |
| Backend Hosting         | Hugging Face Spaces             |
| Frontend                | Streamlit                       |
| Frontend Hosting        | Streamlit Community Cloud       |
| Data Processing         | Pandas / NumPy                  |
| Visualization           | Matplotlib                      |
| Gene Mapping            | Ensembl Gene IDs → Gene Symbols |

---

# Notes

* The input count matrix should contain **raw counts**, not normalized expression values.
* The design formula should use a column available in the run table.
* The sample identifiers in the metadata file must match the sample column names in the count matrix.
* The Gene ID column should not be left blank.
* The control value should exactly match the value present in the target column.
* Ensembl Gene IDs are automatically mapped to gene symbols where possible.
* Processing time depends on dataset size and server availability.

---

# Author

**Jevin Chokshi**

* GitHub: https://github.com/Jevin23
* Streamlit App: https://deseq2-by-jevin.streamlit.app/
* Hugging Face API: https://huggingface.co/spaces/Jevin23/DeSEQ2_API

---

# License

This project is intended for research and educational purposes. Please ensure that your use of the software and any associated datasets complies with applicable licenses and data usage policies.
