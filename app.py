import streamlit as st
import pandas as pd
import matplotlib.pyplot as plt
from Bio import SeqIO
from io import StringIO
from fpdf import FPDF


st.sidebar.title("🧬 FASTA Analyzer")
st.sidebar.markdown("### Navigation")

section = st.sidebar.radio(
    "Go to",
    ["Upload & Overview", "Detailed Analysis", "Export Report"]
)

# --- Helper Functions ---
def read_fasta(uploaded_file):
    stringio = StringIO(uploaded_file.getvalue().decode("utf-8"))
    records = list(SeqIO.parse(stringio, "fasta"))
    return records

    


def gc_content(seq):
    return float(seq.count("G") + seq.count("C")) / len(seq) * 100

def create_report(data, filename="report.pdf"):
    pdf = FPDF()
    pdf.add_page()
    pdf.set_font("Arial", size=12)
    pdf.cell(200, 10, txt="FASTA Analysis Report", ln=True, align="C")
    pdf.ln(10)
    for idx, row in data.iterrows():
        pdf.cell(200, 10, txt=f"{row['ID']} - Length: {row['Length']} - GC%: {row['GC_Content']:.2f}", ln=True)
    pdf.output(filename)



# --- logo ---



st.markdown("""
<style>
body {
    background-color: #0b0f19;
}

.block-container {
    padding: 2rem 3rem;
}

h1, h2, h3 {
    font-weight: 600;
}

[data-testid="stMetric"] {
    background-color: #111827;
    padding: 15px;
    border-radius: 12px;
    text-align: center;
}

.stButton>button {
    border-radius: 8px;
    background: linear-gradient(90deg, #2563eb, #1e40af);
    color: white;
    font-weight: 600;
}

.stDownloadButton>button {
    border-radius: 8px;
    background: #1e293b;
    color: white;
}

hr {
    border: 1px solid #1f2937;
}
</style>
""", unsafe_allow_html=True)




# --- Streamlit UI ---
if section == "Upload & Overview":

    st.title("🧬FASTA Analyzer")
    st.caption("📤Upload and analyze genomic sequences instantly")

    uploaded_file = st.file_uploader("Upload FASTA file", type=["fasta"])

    if uploaded_file:
        with st.spinner("Processing..."):
            records = read_fasta(uploaded_file)

            data = pd.DataFrame({
                "ID": [rec.id for rec in records],
                "Length": [len(rec.seq) for rec in records],
                "GC_Content": [gc_content(rec.seq) for rec in records]
            })

        st.session_state["data"] = data

        st.success("Analysis completed!")

        # --- METRICS ---
        col1, col2, col3 = st.columns(3)
        col1.metric("Sequences", len(data))
        col2.metric("Avg Length", int(data["Length"].mean()))
        col3.metric("Avg GC%", round(data["GC_Content"].mean(), 2))


elif section == "Detailed Analysis":

    st.title("📊Detailed Analysis")

    if "data" not in st.session_state:
        st.warning("Please upload a file first.")
    else:
        data = st.session_state["data"]

        st.dataframe(data, use_container_width=True)

        col1, col2 = st.columns(2)

        with col1:
            st.subheader("📊GC Distribution")
            fig1, ax = plt.subplots()
            ax.hist(data["GC_Content"], bins=10)

            st.pyplot(fig1)

            fig1.savefig("gc_dist.png")
            with open("gc_dist.png", "rb") as f:
                st.download_button("📥Download Chart", f, "gc_distribution.png")

        with col2:
            st.subheader("📊GC Comparison")
            fig2, ax = plt.subplots(figsize = (10, 6.45))
            ax.bar(data["ID"], data["GC_Content"])
            plt.xticks(rotation=60)
            st.pyplot(fig2)

            fig2.savefig("gc_compare.png")
            with open("gc_compare.png", "rb") as f:
                st.download_button("📥Download Chart", f, "gc_comparison.png")

        st.subheader("📈GC Trend")
        fig3, ax = plt.subplots()
        ax.plot(data["GC_Content"], marker='o')
        st.pyplot(fig3)

        fig3.savefig("gc_trend.png")
        with open("gc_trend.png", "rb") as f:
            st.download_button("📥Download Chart", f, "gc_trend.png")


elif section == "Export Report":

    st.title("📤Export Report")

    if "data" not in st.session_state:
        st.warning("No data available. Upload file first.")
    else:
        data = st.session_state["data"]

        if st.button("🔄Generate PDF"):
            create_report(data)
            st.success("✅Report generated!")

            with open("report.pdf", "rb") as f:
                st.download_button("🧾Download PDF", f, "FASTA_Report.pdf")






# ------------------- FOOTER -------------------
st.markdown("---")
st.caption("FASTA Analyzer • Built for genomic data insights")
st.markdown(
    "<hr><center style='color:white !important;'>🧬Powered by BioCode Innovators | Developed by Abeera Iftikhar</center>",
    unsafe_allow_html=True
)