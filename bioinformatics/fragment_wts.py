import re, sys
import streamlit as st

sys.path.append('./bioinformatics')
from reference import amino_acid_weights, weight_of_water

#set initial session state default values
if "sequence" not in st.session_state.keys():
    st.session_state.sequence = ""  
if "target_wt" not in st.session_state.keys():
    st.session_state.target_wt = 0
if "adjust_wt" not in st.session_state.keys():
    st.session_state.adjust_wt = 0
if "delta" not in st.session_state.keys():
    st.session_state.delta = 20

# funciton for reset button
def reset_values():
    st.session_state.sequence = ""  
    st.session_state.target_wt = 0
    st.session_state.adjust_wt = 0
    st.session_state.delta = 20

# function for submit button
def set_values(aa_sequence, target_wt, adjust_wt, delta_wt):
    st.session_state.sequence = aa_sequence
    st.session_state.target_wt = target_wt
    st.session_state.adjust_wt = adjust_wt
    st.session_state.delta = delta_wt

# st.write(st.session_state)

aa_sequence=st.text_area("Protein sequence", key="sequence", value=st.session_state.sequence) 

col1, col2, col3=st.columns(3)
with col1:
    target_wt= st.text_input("Target Weight (Da)", value=st.session_state.target_wt)
with col2:
    adjust_wt = st.text_input("Adjust Weight (Da)", value=st.session_state.adjust_wt)
with col3:
    delta_wt = st.text_input("Error in Weight (Da)", value=st.session_state.delta)

col1, col2, col3, col4, col5=st.columns(5)
with col1:
    # col21, col22=st.columns(2)
    # with col21:
    reset = st.button("Reset", on_click=reset_values)
with col5:
    # with col22:
        submit = st.button("Submit", on_click=set_values, args=(aa_sequence, target_wt, adjust_wt, delta_wt))

st.divider()


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

def filter_fragments(aa_sequence, state_dictionary):
    try:
        sequence = aa_sequence
        target = float(state_dictionary.target_wt) - float(state_dictionary.adjust_wt)
        delta = float(state_dictionary.delta)

        good_fragments= [["start", "end", "size"]]
       

        for start in range(0,len(sequence)):
            for end in range(len(sequence)-1, start-1,-1):
                
                # print(f"start: {start} - end: {end}")
                mass = calculate_mw(sequence[start:end])
                if mass < target-delta:
                    break
                if mass > target+delta:
                    continue
                good_fragments.append([str(start+1), str(end+1), f"{mass: .1f}"])

        return good_fragments
        return [target, delta]
    
    except:
        return None

def format_text(sequence, block_length=10, line_length=50):
    formatted_text = ""
    counter = 0
    while sequence:
        counter += 1
        formatted_text += sequence[0]
        sequence = sequence[1:]
        if counter % 100 == 0:
            formatted_text += "\n\n"
        elif counter % 50 == 0:
            formatted_text += "\n"
        elif counter % 10 == 0:
            formatted_text += " "
    return formatted_text

title, aa_sequence = get_title(aa_sequence)
molecular_weight = calculate_mw(aa_sequence)
just_sequence = "".join([char for char in aa_sequence if char.upper() in "ACDEFGHIKLMNPQRSTVWY"])

with st.expander("Input"):
    st.write(f"{title}")
    st.code(aa_sequence)
    st.write(f"The target molecular weight is {target_wt}")
    st.write(f"{adjust_wt} Da will be added to the calculated molecular weight.")
    st.write(f"A sequence withing {delta_wt} Da of the target molecular weight will be accepted.")

with st.expander("Formatted"):
    st.code(format_text(just_sequence))
    st.write(f"Full length molecular weight: {molecular_weight:.2f}")
    st.write(f"Full length: {len(just_sequence)} residues")

st.subheader("Potential fragments")
st.table(filter_fragments(just_sequence, st.session_state))