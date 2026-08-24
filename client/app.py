import streamlit as st
from components.deseq_ui import render_deseq2_ui


st.set_page_config(layout="wide", page_title="RNA-Seq DESeq2 Pipeline", page_icon="🧬")

st.title("🧬 RNA-Seq Differential Expression Web App")
st.markdown("Run your DESeq2 microservice, visualize quality control metrics, and download differential expression profiles.")
render_deseq2_ui()