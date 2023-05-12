import streamlit as st

st.markdown("# FAQ and Clarifications")
st.markdown("""
##### 1. Could you further clarify how to rate the molecules?
For every monomer and highlighted functional group, you are given the opportunity to:
1. Rate the functional group specifically for the identified polymerization type. I.e., fairly local considerations of the monomer.
2. Rate the monomer overall

The idea is that a functional group might be pretty decent in isolation (in terms of sterics and bulkiness), but the monomer overall might be bad (e.g. too many functional groups, or too bulky). We want to capture these two levels of information.

Additionally, you are also given the opportunity to define other polymerization reactions that a given functional group might be good for, or leave general comments.


Roughly think of the rating scale as follows:
1. Don't bother trying to polymerize it
2. This monomer would need *a lot* of work to polymerize
3. This monomer would need some work to polymerize
4. This monomer needs some work to polymerize, but is interesting enough to pursue
5. This monomer probably works

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