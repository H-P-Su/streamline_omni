import streamlit as st
from bioinformatics import translate

pages = {
    "Utilities":[
        st.Page("Utilities/converter.py", title="Quick Converters")
    ],
    "Bioinformatics":[
        st.Page("bioinformatics/translate.py", title="DNA Translate"), 
        st.Page("bioinformatics/codon_usage.py", title="Calculate Codon Usage"), 
        st.Page("bioinformatics/mol_wt.py", title="Calculate Protein Molecular Weight"),
        st.Page("bioinformatics/fragment_wts.py", title="Id fragments with molecular weight"),            
    ],
}
pg = st.navigation(pages)
# Run the app
pg.run()

