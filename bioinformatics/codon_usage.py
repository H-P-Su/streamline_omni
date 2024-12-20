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


def count_codons(sequence):
    sequence = sequence.upper().replace("T", "U")
    codons = dict()
    for k in genetic_code.keys():
        codons[k] = 0
        
    while len(sequence) > 2:
        codon = sequence[0:3]
        sequence = sequence[3:]
        if codon in codons.keys():
            codons[codon] += 1  
    return codons

def count_codons_percentages(codon_count):
    codon_percentages = dict()
    
    total = 0
    for k,v in codon_count.items():
        total += v

    for k,v in codon_count.items():
        codon_percentages[k] = v/total

    return codon_percentages


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

    with st.expander("Output"):
        st.text(title)
        st.write(translate(dna_sequence))

    codon_count = count_codons(dna_sequence)
    codon_percentages = count_codons_percentages(codon_count)

    codon_stats = []
    codons_by_aa = [[value, key] for key, value in genetic_code.items() if value != "*"]
    codon_stats = sorted(codons_by_aa)
    for i in codon_stats:
        i.append(codon_count[i[1]])
        i.append(100*codon_percentages[i[1]])

    df = pd.DataFrame(codon_stats, columns=["AA", "Codon", "Count", "Percentage"])
    # df['Percentage'] = df["Percentage"].round(1)
    st.table(df)
