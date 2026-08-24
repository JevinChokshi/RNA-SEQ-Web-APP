from modules.deseq2_runner import run_deseq2_microservice
from modules.load_dataset import load_and_process_data
# 1. Import Form from fastapi
from fastapi import APIRouter, UploadFile, Form 
from fastapi.responses import JSONResponse

router = APIRouter()

@router.post('/run_deseq2/')
async def run_deseq2(
    count_file: UploadFile, 
    meta_file: UploadFile, 
    # mapper_file: UploadFile, 
    # 2. Wrap the text/number fields in Form() so FastAPI parses them correctly
    p_adj: float = Form(...),  # p_adj is usually a float (e.g., 0.05)
    log2fc: float = Form(...), # log2fc is usually a float/int
    design: str = Form(...), 
    min_gene_count: int = Form(...), 
    min_samples: int = Form(...), 
    meta_index: str = Form(...), 
    count_index: str = Form(...), 
    label_col: str = Form(...), 
    control_label: str = Form(...)
):
    try:
        mapper_file = 'data/Human.GRCH38.p13.annot.tsv'
        count_df, meta_df, labels, processing_stats = await load_and_process_data(
            meta_file, count_file, min_gene_count=min_gene_count, 
            min_samples=min_samples, meta_index=meta_index, 
            count_index=count_index, label_col=label_col, control_label=control_label
        )
        output_payload = await run_deseq2_microservice(
            counts_df=count_df, meta_df=meta_df, labels=labels, 
            processing_stats=processing_stats, mapper_file=mapper_file,
            p_adj=p_adj, log2fc=log2fc, design=design
        )

        return JSONResponse(status_code=200, content=output_payload)
    except Exception as e:
        return JSONResponse(status_code=500, content={'error': str(e)})
