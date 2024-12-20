import re, sys
import streamlit as st

sys.path.append('./bioinformatics')
from reference import amino_acid_weights, weight_of_water

st.session_state.info = {}
aa_sequence=st.text_area("Protein sequence", value=st.session_state.info.get("sequence", "")) 
submit = st.button("Submit")
if submit:
    st.session_state.info["sequence"] = aa_sequence



def get_title(sequence):
    sequence_lines = sequence.split("\n")
    if re.match(">", sequence_lines[0]): 
        title = sequence_lines[0]
        sequence = "\n".join(sequence_lines[1:])
    else:
        title = None
    
    return title, sequence 

def calculate_mw(sequence):
    sequence = "".join([char for char in sequence if char.upper() in "ACDEFGHIKLMNPQRSTVWY"])
    weight = weight_of_water
    while len(sequence) > 0: 
        aa = sequence[0].upper()
        sequence = sequence[1:]
        if aa in amino_acid_weights.keys():
            weight += amino_acid_weights[aa]
        else:
            print(f"No weight for {aa}")
    return weight

title, aa_sequence = get_title(aa_sequence)
molecular_weight = calculate_mw(aa_sequence)

with st.expander("Input"):
    st.write(f"Title: {title}")
    st.code(aa_sequence)

st.metric("Molecular Weight", value = "{:.2f}".format(molecular_weight))
