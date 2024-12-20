import streamlit as st
from bioinformatics import translate

pages = {
    "Bioinformatics":[
        st.Page("bioinformatics/translate.py", title="DNA Translate"), 
        st.Page("bioinformatics/codon_usage.py", title="Calculate Codon Usage"), 
        st.Page("bioinformatics/mol_wt.py", title="Calculate Protein Molecular Weight"),
        st.Page("bioinformatics/fragment_wts.py", title="Id fragments with molecular weight"),            
    ],
    "Utilities":[
        st.Page("Utilities/converter.py")
    ]
}
pg = st.navigation(pages)
# Run the app
pg.run()

