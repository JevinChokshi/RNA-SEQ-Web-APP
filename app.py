import gradio as gr
import spaces
from pathlib import Path
from server.modules.load_dataset import load_and_process_data
from server.modules.deseq2_runner import run_deseq2_microservice


# ---------------------------------------------------------
# Dummy ZeroGPU function
# This exists only to satisfy Hugging Face ZeroGPU startup
# ---------------------------------------------------------
@spaces.GPU
def ping_gpu():
    return "ZeroGPU is working!"


# ---------------------------------------------------------
# DESeq2 API function
# This remains CPU-based
# ---------------------------------------------------------
async def run_deseq2_api(
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
):
    BASE_DIR = Path(__file__).resolve().parent

    mapper_file = (
        BASE_DIR
        / "server"
        / "data"
        / "Human.GRCh38.p13.annot.tsv"
    )
    

    count_df, meta_df, labels, processing_stats = \
        await load_and_process_data(
            meta_file,
            count_file,
            min_gene_count=min_gene_count,
            min_samples=min_samples,
            meta_index=meta_index,
            count_index=count_index,
            label_col=label_col,
            control_label=control_label
        )

    result = await run_deseq2_microservice(
        counts_df=count_df,
        meta_df=meta_df,
        labels=labels,
        processing_stats=processing_stats,
        mapper_file=str(mapper_file),
        p_adj=p_adj,
        log2fc=log2fc,
        design=design
    )

    return result


# ---------------------------------------------------------
# Gradio API interface
# ---------------------------------------------------------
demo = gr.Interface(
    fn=run_deseq2_api,

    inputs=[
        gr.File(type="filepath", label="Count Matrix"),
        gr.File(type="filepath", label="Metadata"),

        gr.Number(
            value=0.05,
            label="Adjusted p-value threshold"
        ),

        gr.Number(
            value=1.5,
            label="Log2 Fold Change threshold"
        ),

        gr.Textbox(
            value="~disease_state",
            label="Design"
        ),

        gr.Number(
            value=10,
            label="Minimum gene count"
        ),

        gr.Number(
            value=2,
            label="Minimum samples"
        ),

        gr.Textbox(
            value="sample",
            label="Metadata index"
        ),

        gr.Textbox(
            value="gene",
            label="Count index"
        ),

        gr.Textbox(
            value="condition",
            label="Label column"
        ),

        gr.Textbox(
            value="Control",
            label="Control label"
        ),
    ],

    outputs=gr.JSON(),

    api_name="run_deseq2",

    title="PyDESeq2 Analysis API",

    description=(
        "Run DESeq2 analysis and generate DEG results, "
        "VST counts, PCA and volcano plot outputs."
    )
)


# ---------------------------------------------------------
# Queue
# ---------------------------------------------------------
demo.queue(default_concurrency_limit=1)

if __name__ == "__main__":
    demo.launch()