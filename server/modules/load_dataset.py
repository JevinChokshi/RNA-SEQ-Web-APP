import io
import pandas as pd
from dotenv import load_dotenv
from server.modules.validators import validate_dataframe_not_empty, validate_required_columns, validate_numeric_counts, validate_sample_overlap
from fastapi import UploadFile

load_dotenv()

async def load_and_process_data(
    meta_file: UploadFile,
    counts_file: UploadFile,
    min_gene_count: int,
    min_samples: int,
    meta_index: str,
    count_index: str,
    label_col: str,
    control_label: str
):
    # 1. Asynchronously read raw bytes from the FastAPI UploadFile objects
    if isinstance(meta_file, str):
        with open(meta_file, "rb") as f:
            meta_bytes = f.read()
    else:
        meta_bytes = await meta_file.read()


    if isinstance(counts_file, str):
        with open(counts_file, "rb") as f:
            counts_bytes = f.read()
    else:
        counts_bytes = await counts_file.read()

    # 2. Wrap bytes in io.BytesIO so pandas can parse them natively
    meta_df = pd.read_csv(
        io.BytesIO(meta_bytes),
        index_col=meta_index
    )

    counts_df = pd.read_csv(
        io.BytesIO(counts_bytes),
        sep='\t',
        index_col=count_index
    )

    # --- Rest of your logic stays exactly the same ---
    validate_dataframe_not_empty(meta_df, "Metadata")
    validate_dataframe_not_empty(counts_df, "Counts")

    validate_required_columns(
        meta_df,
        [label_col],
        "Metadata"
    )

    original_samples = meta_df.shape[0]
    original_genes = counts_df.shape[0]

    meta_df.rename(
        columns={label_col: 'disease_state'},
        inplace=True
    )

    meta_df.loc[
        meta_df['disease_state'] == control_label,
        'disease_state'
    ] = 'Control'

    meta_df = meta_df.dropna(
        subset=['disease_state']
    )

    labels = meta_df['disease_state'].unique().tolist()

    meta_df['disease_state'] = pd.Categorical(
        meta_df['disease_state'],
        categories=labels,
        ordered=True
    )

    meta_df = meta_df[
        ~meta_df.index.duplicated(keep="first")
    ]

    validate_sample_overlap(
        counts_df,
        meta_df
    )

    common_samples = counts_df.columns.intersection(
        meta_df.index
    )

    counts_df = counts_df[common_samples]
    meta_df = meta_df.loc[common_samples]

    counts_df = counts_df.loc[
        (counts_df.sum(axis=1) > 0),
        :
    ]

    counts_df = counts_df[
        (counts_df >= min_gene_count).sum(axis=1)
        >= min_samples
    ]

    validate_numeric_counts(counts_df)

    filtered_genes = (
        original_genes - counts_df.shape[0]
    )

    removed_samples = (
        original_samples - meta_df.shape[0]
    )

    processing_stats = {
    "original_samples": int(original_samples),
    "samples_retained": int(meta_df.shape[0]),
    "samples_removed": int(removed_samples),
    "original_genes": int(original_genes),
    "genes_retained": int(counts_df.shape[0]),
    "genes_filtered": int(filtered_genes)
}

    return (
        counts_df,
        meta_df,
        labels,
        processing_stats
    )