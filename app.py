import streamlit as st
import pandas as pd
import matplotlib.pyplot as plt
from Bio import SeqIO
from io import StringIO
from fpdf import FPDF


st.sidebar.title("🧬 FASTA Analyzer")
st.sidebar.markdown("### Dashboard Menu")

section = st.sidebar.radio(
    "Go to",
    ["📤 Upload & Overview", "📊 Detailed Analysis", "🧾Export Report", "ℹ️About Us"]
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



# --- CSS ---


# --- Hide Streamlit default header & footer ---
hide_streamlit_style = """
    <style>
    #MainMenu {visibility: hidden;}     /* Hamburger menu */
    header {visibility: hidden;}        /* Top header */
    footer {visibility: hidden;}        /* Footer */
    </style>
"""
st.markdown(hide_streamlit_style, unsafe_allow_html=True)

st.markdown("""
<style>
@import url('https://fonts.googleapis.com/css2?family=Inter:wght@300;400;600&family=Space+Grotesk:wght@500;700&display=swap');
""", unsafe_allow_html=True)


st.markdown("""
<style>

/* Import fonts */
@import url('https://fonts.googleapis.com/css2?family=Inter:wght@300;400;600&family=Space+Grotesk:wght@500;700&display=swap');

/* Global font */
html, body, [class*="css"] {
    font-family: 'Inter', sans-serif;
}

/* MAIN TITLE (FASTA Analyzer) */
h1 {
    font-family: 'Space Grotesk', sans-serif;
    font-weight: 700;
    letter-spacing: 0.5px;
}

/* Subheaders */
h2, h3 {
    font-family: 'Inter', sans-serif;
    font-weight: 600;
}

/* Sidebar text */
[data-testid="stSidebar"] {
    font-family: 'Inter', sans-serif;
}

</style>
""", unsafe_allow_html=True)



st.markdown("""
<style>

/* Main App Background */
[data-testid="stAppViewContainer"] {
    background-color: #0b0f19;
}

/* Sidebar */
[data-testid="stSidebar"] {
    background-color: #111827;
}

/* Titles */
h1, h2, h3 {
    color: #ffffff;
    font-weight: 600;
}

/* Buttons */
button {
    border-radius: 8px !important;
    background: linear-gradient(90deg, #2563eb, #1e40af) !important;
    color: white !important;
    font-weight: 600 !important;
}

/* Metrics Cards */
[data-testid="stMetric"] {
    background-color: #111827;
    padding: 15px;
    border-radius: 12px;
    text-align: center;
}

/* Dataframe */
[data-testid="stDataFrame"] {
    border-radius: 10px;
}

/* Section spacing */
.block-container {
    padding-top: 2rem;
    padding-left: 3rem;
    padding-right: 3rem;
}

/* Divider */
hr {
    border: 1px solid #1f2937;
}

/* Footer */
footer {
    visibility: hidden;
}

</style>
""", unsafe_allow_html=True)






# --- Streamlit UI ---
if section == "📤 Upload & Overview":

    st.markdown("""
    <h1 style='font-family: Space Grotesk;'>
    🧬 FASTA Analyzer
    </h1>
    """, unsafe_allow_html=True)
    st.caption("📤Upload and analyze genomic sequences instantly")

    uploaded_file = st.file_uploader("Upload FASTA file", type=["fasta"])

    if uploaded_file:
        with st.spinner("Processing..."):
            records = read_fasta(uploaded_file)

            data = pd.DataFrame({
                "Protein_Name": [rec.description for rec in records],
                "ID": [rec.id for rec in records],
                "Length": [len(rec.seq) for rec in records],
                "GC_Content": [gc_content(rec.seq) for rec in records]
            })

        st.session_state["data"] = data
        st.dataframe(data)

        st.success("Analysis completed!")


        # --- CSV Download Button ---
        csv = data.to_csv(index=False).encode("utf-8")
        st.download_button(
            label="📥 Download CSV",
            data=csv,
            file_name="FASTA_Analysis.csv",
            mime="text/csv"
        )

        # --- METRICS ---
        col1, col2, col3 = st.columns(3)
        col1.metric("Sequences", len(data))
        col2.metric("Avg Length", int(data["Length"].mean()))
        col3.metric("Avg GC%", round(data["GC_Content"].mean(), 2))


elif section == "📊 Detailed Analysis":

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
            ax.set_title("GC Content Distribution")
            ax.set_xlabel("GC%")
            ax.set_ylabel("Frequency")

            st.pyplot(fig1)

            fig1.savefig("gc_dist.png")
            with open("gc_dist.png", "rb") as f:
                st.download_button("📥Download Chart", f, "gc_distribution.png")

        with col2:
            st.subheader("📊GC Comparison")
            fig2, ax = plt.subplots(figsize = (10, 5.8))
            ax.set_title("GC Content Comparison Across Multiple FASTA Files")
            ax.set_xlabel("Gene")
            ax.set_ylabel("GC Content (%)")
            ax.bar(data["ID"], data["GC_Content"])
            plt.xticks(rotation=60)
            st.pyplot(fig2)

            fig2.savefig("gc_compare.png")
            with open("gc_compare.png", "rb") as f:
                st.download_button("📥Download Chart", f, "gc_comparison.png")

        st.subheader("📈GC Trend")
        fig3, ax = plt.subplots()
        ax.plot(data["GC_Content"], marker='o')
        ax.set_title("GC Content Trend Across Sequences")
        ax.set_xlabel("Sequence Index")
        ax.set_ylabel("GC %")
        st.pyplot(fig3)

        fig3.savefig("gc_trend.png")
        with open("gc_trend.png", "rb") as f:
            st.download_button("📥Download Chart", f, "gc_trend.png")


        # --- New Section: Additional Plots ---
        st.subheader("📊Additional GC Content Analysis")

        fig4, axs = plt.subplots(2, 2, figsize=(14,10))

        # Top-left: GC Content (%) per gene
        axs[0,0].bar(data["ID"], data["GC_Content"])
        axs[0,0].set_title("GC Content (%)")
        axs[0,0].set_ylabel("GC %")
        axs[0,0].tick_params(axis='x', rotation=60)

        # Top-right: Gene Length
        axs[0,1].bar(data["ID"], data["Length"])
        axs[0,1].set_title("Gene Length")
        axs[0,1].set_ylabel("Length")
        axs[0,1].tick_params(axis='x', rotation=60)

        # Bottom-left: GC Content Distribution
        axs[1,0].hist(data["GC_Content"], bins=10)
        axs[1,0].set_title("GC Content Distribution")
        axs[1,0].set_xlabel("GC %")
        axs[1,0].set_ylabel("Frequency")

        # Bottom-right: GC Content Spread (Boxplot)
        axs[1,1].boxplot(data["GC_Content"], patch_artist=True, boxprops=dict())
        axs[1,1].set_title("GC Content Spread")
        axs[1,1].set_ylabel("GC %")

        plt.tight_layout()
        st.pyplot(fig4)

        fig4.savefig("additional_plots.png")
        with open("additional_plots.png", "rb") as f:
            st.download_button("📥Download Chart", f, "additional_plots.png")


elif section == "🧾Export Report":

    st.title("🧾Export Report")

    if "data" not in st.session_state:
        st.warning("No data available. Upload file first.")
    else:
        data = st.session_state["data"]

        if st.button("Generate PDF"):
            create_report(data)
            st.success("Report generated!")

            with open("report.pdf", "rb") as f:
                st.download_button("Download PDF", f, "FASTA_Report.pdf")









elif section == "ℹ️About Us":

    # --- New Section: About Biocode Innovators ---
    st.markdown("---")  # separator line

    st.header("About Biocode Innovators")

    st.markdown("""
    **Biocode Innovators** bridges the gap between biology and computational technology.  
    We provide hands-on training, workshops, and resources for students, researchers, and professionals in bioinformatics, genomics, and data-driven life sciences.  

    **Founded by:** Abeera Iftikhar  
    **Mission:** Make biological data accessible, understandable, and actionable.  
    **Slogan:**  *Let's decode biology together!*  

    **Areas of Focus:**
    - Computational Biology  
    - Data Analysis  
    - Biological Data Visualization  
    - In-silico Analysis  
    - Python for Biology  
    - R for Biology  

    **Industry:** E-Learning Providers  
    **Company Size:** 2–10 employees (4 associated members)  
    **Founded:** 2025  
    **Phone:** 03159633608  
    """)







# ------------------- FOOTER -------------------
st.markdown("---")
st.caption("FASTA Analyzer • Built for genomic data insights")
st.markdown(
    "<hr><center style='color:white !important;'>🧬Powered by BioCode Innovators | Developed by Abeera Iftikhar</center>",
    unsafe_allow_html=True
)