import streamlit as st
st.set_page_config(layout="wide")
import pandas as pd
import os, time
from rdkit import Chem
from rdkit.Chem import Descriptors
import st_utils, app_utils
import gspread_pdlite as gspdl # database connection

if "settings_initialized" not in st.session_state:
    app_utils.initialize()
    st.session_state["settings_initialized"] = True
st_utils.set_font(widget=16)


# ===== Setup =====
def set_update_data_flag(flag):
    st.session_state["b_update_data"] = flag
    st.session_state["b_update_data_batch"] = flag


# ===== ===== BEGIN LAYOUT ===== =====
st.markdown("# Welcome!")
cols = st.columns(2)
with cols[0]:
    st.markdown("- On this page, you can choose some basic filters for the molecules to view, as well as leave general comments about the process, monomers, etc.")
    st.markdown("- Ratings can be given for for the polymerization potential of both specific functional groups as well as the monomer overall.")
    st.markdown("- The scale is `1 (impossible)` - `3 (potentially workable)` - `5 (probably works)`. There is also an option to skip answering, or mark an identified (functional group, reaction) pairing as N/A.")
with cols[1]:
    st.markdown("- On the `Batch Evaluation` page, you can rate multiple functional groups and molecules at a time.")
    st.markdown("- On the `Evaluation` page, you will rate one molecule at a time, and have the opportunity to give more thorough assessments of the molecule.")
    st.markdown("- On the `Results` tab, you can view your top-rated molecules, and optionally download a `.csv` of all molecules that you've rated.")

#st.markdown("- On the `Evaluation` tab (in the sidebar), you will be shown monomers, with one highlighted functional group at a time, and the opportunity to rate the suitability of the functional groups for polymerization. \n  - *The first load may take a few minutes!* \n")


cols = st.columns(2)

with cols[0]:
    user_info = st_utils.persist_widget( st.text_input, "email",
                                key = "userinfo", val0="",
                                on_change = lambda: None)

    if "userinfo" in st.session_state \
        and st.session_state["userinfo"] not in ["",None] \
        and "@" in st.session_state["userinfo"]:
        load_previous = st.button("(Optional): load all previous session results",
                                on_click=app_utils.load_user_session)

    # ===== Navigation =====
    #st.markdown("## Navigation & Settings")
    #st.markdown("### On next molecule, show...")
    rxn_selection = st_utils.persist_widget( st.selectbox, "polymerization motif",
                                key = "rxn_selection", val0 = "choose for me!",
                                options = ("choose for me!",*st.session_state["rxn_types"]),
                                on_change = lambda: set_update_data_flag(True))

    iteration_selection = st_utils.persist_widget( st.selectbox,
                                        "molecule presentation mode: (sequential or random)",
                                        options = ("random","sequential"),
                                        key="iteration_selection",
                                        val0 = "random",
                                        on_change = lambda: set_update_data_flag(True))

with cols[1]:
    # MW
    slider_MW = st_utils.persist_widget(st.slider, "MW range",
                                min_value = 0., max_value = st.session_state.max_MW,
                                key = "slider_MW",
                                val0 = (0.,250.),
                                on_change = lambda: set_update_data_flag(True) )

    # Simplicity/Complexity
    # (# of polymerizations identified, # functional groups, # subsitutents)
    slider_num_ftn_specific = st_utils.persist_widget(st.slider,"Number of selected potentially polymerizable functional groups. If no polymerization type selected, same as next filter.",
                                min_value = 1, max_value = st.session_state.max_numftn,
                                on_change = lambda: set_update_data_flag(True),
                                key="slider_num_ftn_specific",
                                val0 = (1,5) )

    slider_num_ftn = st_utils.persist_widget(st.slider,"Number of (any) potentially polymerizable functional groups",
                                min_value = 1, max_value = st.session_state.max_numftn,
                                on_change = lambda: set_update_data_flag(True),
                                key="slider_num_ftn",
                                val0 = (1,5) )

# Bulkiness

# ===== Comments =====
# Comments/bugs
st.markdown("---")

def log_general_comment():
    """_summary_
    Notes:
        Requires the following: 
    """
    comment_dict = {}
    molid, molidx = st.session_state["data_index"]
    comment_dict["molid"],comment_dict["molidx"] = molid, molidx
    comment_dict["smiles"] = st.session_state["data_rxn"].loc[molid].smiles
    comment_dict["userinfo"] = st.session_state["userinfo"]
    comment_dict["timestamp"] = pd.Timestamp(time.time(),unit="s")
    comment_dict["comments_general"] = st.session_state["comments_general"]

    st.session_state["eval_general"] = pd.concat([ st.session_state["eval_general"], 
                                                    pd.DataFrame([comment_dict]) ], ignore_index=True)
    
    # Save
    st.session_state["eval_general"].to_csv("eval_general.csv",index=False)

    ws    = st.session_state["backend_sheet"].worksheet( st.secrets.data.name_general )
    gspdl.worksheet_append( ws, pd.Series(comment_dict) )

with st.sidebar:
    if "userinfo" not in st.session_state \
        or st.session_state["userinfo"] in ["",None] \
        or "@" not in st.session_state["userinfo"]:
        st.markdown("##### Please enter a valid e-mail to the right in order to submit general comments")
    else:
        with st.form("general comments/bugs",clear_on_submit=True):
            comment_area = st.text_area("General comments?","",key="comments_general")
            submitted = st.form_submit_button("submit",on_click=log_general_comment)

# Preload data
app_utils.first_load()
app_utils.update_filters()
st.markdown(f"#### Molecules left after filtering: `{st.session_state.prev_data.index.unique('molid').size}`")
st.markdown("**Distribution of potential polymerization motifs:**")

# Visualize
counts = []
multi_index, data_rxn = st.session_state["prev_data"], st.session_state["data_rxn"]
molids = multi_index.index.unique("molid")
filtered_data_rxn = data_rxn.loc[molids]
for rxn_name in st.session_state["rxn_types"]:
    valid_molids = app_utils.get_valid_reactions( filtered_data_rxn, (1,50), rxn_name )
    counts.append(valid_molids.size)

import numpy as np
ks = np.array( st.session_state["rxn_types"] )
vs = np.array(counts)
sort_index = np.argsort( np.array(counts) )
import matplotlib.pyplot as plt
cols = st.columns(2)

with cols[0]:
    # Top 10
    fig, ax = plt.subplots(figsize=(2, 2))
    y_pos = np.arange(len(counts))
    rs = st.session_state["rxn_types"]
    ktmp,vtmp = ks[sort_index[-1:-11:-1]], vs[sort_index[-1:-11:-1]] 
    ytmp = np.arange(len(ktmp))
    ax.barh(ytmp,vtmp)
    ax.set_yticks(ytmp, ktmp,fontsize=6)
    ax.invert_yaxis()
    ax.set_xlabel("count")
    st.pyplot(fig)

with cols[1]:
    # Bottom
    fig, ax = plt.subplots(figsize=(2, 2))
    y_pos = np.arange(len(counts))
    rs = st.session_state["rxn_types"]
    ktmp,vtmp = ks[sort_index[-11::-1]], vs[sort_index[-11::-1]] 
    ytmp = np.arange(len(ktmp))
    ax.barh(ytmp,vtmp)
    ax.set_yticks(ytmp, ktmp,fontsize=6)
    ax.invert_yaxis()
    ax.set_xlabel("count")
    st.pyplot(fig)
