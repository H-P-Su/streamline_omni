import streamlit as st

st.header("Converters")

st.subheader("Temperature")
col1, col2, col3 = st.columns([2,2,2])

if "degree_F" not in st.session_state.keys():
    st.session_state.degree_F = ""
if "degree_C" not in st.session_state.keys():
    st.session_state.degree_C = ""
if "degree_K" not in st.session_state.keys():
    st.session_state.degree_K = "" 
if "miles" not in st.session_state.keys():
    st.session_state.miles = "" 
if "km" not in st.session_state.keys():
    st.session_state.km = "" 
if "feet" not in st.session_state.keys():
    st.session_state.feet = "" 
if "inch" not in st.session_state.keys():
    st.session_state.inch = "" 
if "cm" not in st.session_state.keys():
    st.session_state.cm = "" 

def reset_page():
    st.session_state.degree_F = ""
    st.session_state.degree_C = ""
    st.session_state.degree_K = "" 
    st.session_state.miles = "" 
    st.session_state.km = "" 
    st.session_state.feet = "" 
    st.session_state.inch = "" 
    st.session_state.cm = "" 
    return None


def update_temp_F():
    if len(st.session_state.degree_F) == 0:
        f = 0 
    else:
        f = float(st.session_state.degree_F)
    c = (f-32) *  (5/9)
    k = c + 273.15
    st.session_state.degree_C = f"{c: .2f}"
    st.session_state.degree_K = f"{k: .2f}"
    return None

def update_temp_C():
    if len(st.session_state.degree_C) == 0:
        c = 0 
    else:
        c = float(st.session_state.degree_C)

    f = (c * 9/5) + 32 
    k = c + 273.15
    st.session_state.degree_F = f"{f: .2f}"
    st.session_state.degree_K = f"{k: .2f}"
    return None

def update_temp_K():
    if len(st.session_state.degree_K) == 0:
        k = 0 
    else:
        k = float(st.session_state.degree_K)

    c = k - 273.15
    f = (c * 9/5) +32 

    st.session_state.degree_C = f"{c: .2f}"
    st.session_state.degree_F = f"{f: .2f}"
    return None

col1.text_input("Fahrenheit", key="degree_F", value = st.session_state.degree_F, on_change=update_temp_F)
col2.text_input("Celsius", key="degree_C", value = st.session_state.degree_C, on_change=update_temp_C)
col3.text_input("Kelvin", key="degree_K", value = st.session_state.degree_K, on_change=update_temp_K)

st.divider()
st.subheader("Distances")

col1, col2 = st.columns([2,2])

def update_temp_miles():
    m = float(st.session_state.miles)
    st.session_state.km = f"{1.60934 * m: .2f}"
    return None

def update_temp_km():
    km = float(st.session_state.km)
    st.session_state.miles = f"{km/1.60934: .2f}"
    return None

col1.text_input("Miles", key="miles", value = st.session_state.miles, on_change=update_temp_miles)
col2.text_input("Kilometers", key="km", value = st.session_state.km, on_change=update_temp_km)

col1, col2, col3 = st.columns([1, 1, 2])


def update_temp_inch():
    if len(st.session_state.feet) == 0:
        f = 0
    else:
        f = float(st.session_state.feet)

    if len(st.session_state.inch) == 0:
        i = 0
    else:
        i = float(st.session_state.inch)

    st.session_state.cm = f"{2.54 * (f * 12 + i): .2f}"
    return None

def update_temp_cm():
    cm = float(st.session_state.cm)
    feet = cm/2.54 // 12
    inch = cm/2.54 % 12
    st.session_state.feet = f"{int(feet)}"
    st.session_state.inch = f"{inch: .2f}"
    return None

col1.text_input("Feet", key="feet", value = st.session_state.feet, on_change=update_temp_inch)
col2.text_input("Inches", key="inch", value = st.session_state.inch, on_change=update_temp_inch)
col3.text_input("Centimeters", key="cm", value = st.session_state.cm, on_change=update_temp_cm)

st.button("Reset", on_click=reset_page)