import io
import pandas as pd
import numpy as np

async def map_genes(results: pd.DataFrame, mapper_file):

    # Gradio sends mapper_file as a file path / NamedString
    # with open(str(mapper_file), "rb") as f:
    #     mapper_bytes = f.read()

    mapper = pd.read_csv(
        mapper_file,
        sep="\t",
        low_memory=False
    )

    mapper = mapper[
        [
            "GeneID",
            "Symbol",
            "EnsemblGeneID",
            "Description",
            "GeneType"
        ]
    ].drop_duplicates()

    mapper["GeneID"] = mapper["EnsemblGeneID"].astype(str)

    results = results.reset_index()

    results["GeneID"] = results["GeneID"].astype(str)

    results_mapped = results.merge(
        mapper,
        on="GeneID",
        how="left"
    )

    results_mapped["direction"] = np.where(
        results_mapped["log2FoldChange"] > 0,
        "UP",
        "DOWN"
    )

    return results_mapped