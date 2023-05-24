import streamlit as st

# ===== States =====
import app_utils
if "settings_initialized" not in st.session_state:
    app_utils.initialize()
    st.session_state["settings_initialized"] = True
st.session_state.state["reload_batch_evaluation"] = True

# ===== Main =====

st.markdown("# FAQ and Clarifications")
st.markdown("""
##### 1. Could you further clarify how to rate the molecules?
If you do choose to rate the monomers, for every monomer and highlighted functional group, you are given the opportunity to:
1. Rate the functional group specifically for the identified polymerization type. I.e., fairly local considerations of the monomer.
2. Rate the monomer overall

You can always `skip` or mark a given (functional group, reaction) pairing as `incorrect`.

The idea is that a functional group might be pretty decent in isolation (in terms of sterics and bulkiness), but the monomer overall might be bad (e.g. too many functional groups, or too bulky). We want to capture these two levels of information.

Additionally, you are also given the opportunity to define other polymerization reactions that a given functional group might be good for, or leave general comments.

Roughly think of the rating scale as follows:  
`5`. Definitely works  
`4`.   
`3`. Potentially works (probably involves a good amount of modification to the monomer)  
`2`. 
`1`. Impossible (don't bother trying to polymerize it)



**IMPORTANT:** If the polymerization scheme you are thinking of involves significantly modifying or changing the monomer, there is an option in the general comments field that you should select.


##### 2. How are the ring-opening polymerizations defined?
The ring-opening polymerization motifs were taken from the [polymer database](https://polymerdatabase.com/polymer%20chemistry/ROMP_table.html).  

##### 3. I noticed something buggy or confusing, where should I report it?
Go to the `Settings` page and submit a general comment.
Other options:
- submit an issue on [github](https://github.com/EqualAPriori/polymaps/issues) (preferred).
- e-mail <kevinshen@ucsb.edu>




""")

st.markdown("---")
cols = st.columns(2)
with cols[0]:
    st.markdown("# Updates")
    st.markdonw("""### 2023.05.23
- full (automatic) loading for all comments
- updated comment fields
""")
    st.markdown("""### 2023.05.17
- Streamlined UX for clarity
- Simplified single molecule page into view and comment
""")
    st.markdown("""### 2023.05.10  
- Updated batch evaluation interface
- Refined functional group evaluation for cleaner data set
- Refactored to run faster
""")
            
with cols[1]:
    st.markdown("# Wishlist")
    st.markdown("""### 2023.05.10  
- Divide up dataset into even blocks (w/ appropriate settings) for people to evaluate
- Tutorial/guide (e.g. case studies to illustrate how to rate)
- Monomer -> polymer and property prediction
- Overview/statistics of aggregated ratings
""")