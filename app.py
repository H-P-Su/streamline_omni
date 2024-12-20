import streamlit as st
import home, page1, page2
from bioinformatics import translate

pg = st.navigation([st.Page("home.py"),
                    st.Page("page1.py"), 
                    st.Page("page2.py"), 
                    st.Page("bioinformatics/translate.py"), 
                    st.Page("bioinformatics/mol_wt.py"), 
                    st.Page("bioinformatics/fragment_wts.py"), 
                    ])
# Run the app
pg.run()

