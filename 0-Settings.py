import streamlit as st
import pandas as pd
import os, time
from rdkit import Chem
from rdkit.Chem import Descriptors
import st_utils, app_utils
import gspread_pdlite as gspdl # database connection

st.session_state["settings_initialized"] = True
st_utils.set_font(widget=16)


# ===== Setup =====
if "b_update_data" not in st.session_state:
    st.session_state["b_update_data"] = False
if "b_update_data_batch" not in st.session_state:
    st.session_state["b_update_data_batch"] = False
def set_update_data_flag(flag):
    st.session_state["b_update_data"] = flag
    st.session_state["b_update_data_batch"] = flag


# ===== ===== BEGIN LAYOUT ===== =====
st.markdown("# Welcome!")
st.markdown("- On this page, you can choose some basic filters for the molecules to view, as well as leave general comments about the process, monomers, etc.  \n")
st.markdown("- On the `Evaluation` tab (in the sidebar), you will be shown monomers, with one highlighted functional group at a time, and the opportunity to rate the suitability of the functional groups for polymerization. \n  - Rating is on a scale of `1 (bad)` - `3 (interesting)` - `5 (good)`. Default is `0 to skip` judging. \n  - *The first load may take a few minutes!* \n")
st.markdown("- On the `Results` tab, you can view your top-rated molecules, and optionally download a `.csv` of all molecules that you've rated.")

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
                                     val0 = "sequential",
                                     on_change = lambda: set_update_data_flag(True))

# MW
slider_MW = st_utils.persist_widget(st.slider, "MW range",
                            min_value = 0., max_value = st.session_state.max_MW,
                            key = "slider_MW",
                            val0 = (0.,250.),
                            on_change = lambda: set_update_data_flag(True) )

# Simplicity/Complexity
# (# of polymerizations identified, # functional groups, # subsitutents)
slider_num_ftn_specific = st_utils.persist_widget(st.slider,"number of selected (above) potentially polymerizable functional groups. If none selected, defaults to any identified functinoal group.",
                            min_value = 1, max_value = st.session_state.max_numftn,
                            on_change = lambda: set_update_data_flag(True),
                            key="slider_num_ftn_specific",
                            val0 = (1,3) )

slider_num_ftn = st_utils.persist_widget(st.slider,"number of (any) potentially polymerizable functional groups",
                            min_value = 1, max_value = st.session_state.max_numftn,
                            on_change = lambda: set_update_data_flag(True),
                            key="slider_num_ftn",
                            val0 = (1,4) )

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
