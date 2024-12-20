import re, sys
import streamlit as st
import pandas as pd
from Bio import SeqIO

sys.path.append('./bioinformatics')
from reference import genetic_code

st.header("Translate DNA to Protein")

st.session_state.info = {}
dna_sequence=st.text_area("DNA sequence", value=st.session_state.info.get("sequence", ""))  
submit = st.button("Submit")
if submit:
    st.session_state.info["sequence"] = dna_sequence

st.divider()

def translate(sequence):
    sequence = sequence.upper().replace("T", "U")
    aa_sequence = ""
    while len(sequence) > 2:
        codon = sequence[0:3]
        sequence = sequence[3:]
        if codon in genetic_code.keys():
            aa_sequence += genetic_code[codon]
        else:   
            print(codon)
    return aa_sequence


def get_title(dna_sequence):
    dna_lines = dna_sequence.split("\n")
    if re.match(">", dna_lines[0]): 
        title = dna_lines[0]
        dna_sequence = "\n".join(dna_lines[1:])
    else:
        title = None
    dna_sequence = "".join([char for char in dna_sequence if char in "acgtuACGTU"])
    
    return title, dna_sequence 


if dna_sequence:
    with st.expander("Input"):
        st.code(dna_sequence)

    title, dna_sequence = get_title(dna_sequence)
    st.subheader("Translation")
    st.text(title)
    st.text(translate(dna_sequence))


