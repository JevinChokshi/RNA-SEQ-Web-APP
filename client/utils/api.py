import os
import tempfile
from config import SPACE_ID
from gradio_client import Client, handle_file




client = Client(SPACE_ID)


def deseq2_api(
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
    temp_files = []

    try:
        # Create temporary local files from Streamlit UploadedFile objects
        for uploaded_file in [count_file, meta_file]:

            suffix = os.path.splitext(uploaded_file.name)[1]

            temp = tempfile.NamedTemporaryFile(
                delete=False,
                suffix=suffix
            )

            temp.write(uploaded_file.getvalue())
            temp.close()

            temp_files.append(temp.name)

        # Send temporary file paths to Gradio
        result = client.predict(
            count_file=handle_file(temp_files[0]),
            meta_file=handle_file(temp_files[1]),
            p_adj=p_adj,
            log2fc=log2fc,
            design=design,
            min_gene_count=min_gene_count,
            min_samples=min_samples,
            meta_index=meta_index,
            count_index=count_index,
            label_col=label_col,
            control_label=control_label,
            api_name="/run_deseq2",
        )

        return result

    finally:
        # Delete temporary files after upload to Gradio
        for file_path in temp_files:
            try:
                os.remove(file_path)
            except OSError:
                pass