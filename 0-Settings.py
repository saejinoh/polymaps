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
    st.session_state["b_update_data_single"] = flag

st.markdown("""
<style>
.big-font {
    font-size:20px !important;
}
</style>
""", unsafe_allow_html=True)


from PIL import Image
image = Image.open("bio-pacific-announcement-banner.jpg")
st.image(image, caption='BioPACIFIC Banner')


# ===== ===== BEGIN LAYOUT ===== =====
st.markdown("# Welcome to MAPS!")
intro_text = """##### The goal of MAPS (see [slides](https://docs.google.com/presentation/d/1wbf9w7fOz857bkP1-JDzYartVtSpfVswCPmVemt8rW4/edit?usp=sharing)) is to crowd source the ***first*** chemical dataset with an explicit focus on polymerization feasibility. 
##### The molecules here are pulled from the [natural product database](https://www.npatlas.org/). 
"""
###### to prospectively find natural products that might be polymerized into interesting new materials. 
##### This project will
###### , and 2) can be used by scientists to conduct initial explorations of natural products, with explicit focus on polymerization.
###### However, to date there is ***no*** dataset that documents polymerization reactions, rendering it difficult to do reliable computational screens of the large natural product space. 
###### On this page there are some settings to filter that data set by polymerization functional motif, MW, etc.  You can submit your feedback on the evaluation pages.
#st.markdown('<p class="big-font">'+intro_text+'</p>',unsafe_allow_html=True)
st.markdown(intro_text)

cols = st.columns(2)
with cols[0]:
    #expander_label = "Details -- what do I do with this webapp?"
    #with st.expander(expander_label, expanded = False):
        st.markdown("- On this page, you can choose some basic filters for the molecules to view, as well as leave general comments about the process, etc.")
        st.markdown("- Ratings can be given for the polymerization potential of both specific functional groups as well as the monomer overall.")
        st.markdown("- You can always `skip` answering, or mark an identified pairing of (functional group, reaction) as `incorrectly ID'd`.")
        st.markdown("- Scale: `1 (impossible)` - `3 (potentially works)` - `5 (definitely works)`.")
    #st_utils.change_widget_fontsize(expander_label,"18px")
with cols[1]:
    #expander_label = "What are all the pages in the side bar?"
    #with st.expander(expander_label, expanded = False):
        st.markdown("- On the `Batch Evaluation` page, you can rate multiple functional groups and molecules at a time.")
        st.markdown("- On the `Single Molecule Viewer` page, you can view individual molecules one at a time (e.g. at random or you can input a molecule ID).")
        st.markdown("- On the `Results` page, you can view your top-rated molecules, and optionally download a `.csv` of all molecules that you've rated.")
        st.markdown("- On the `FAQ` page you'll find additional information.")
    #st_utils.change_widget_fontsize(expander_label,"18px")
#st.markdown("- On the `Evaluation` tab (in the sidebar), you will be shown monomers, with one highlighted functional group at a time, and the opportunity to rate the suitability of the functional groups for polymerization. \n  - *The first load may take a few minutes!* \n")


#expander_label = "Username and Filters"
#with st.expander(expander_label, expanded=False):
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
    slider_num_ftn_specific = st_utils.persist_widget(st.slider,f"Number of `{st.session_state.rxn_selection}` polymerizable functional groups. (If no polymerization type selected, same as next filter.)",
                                min_value = 1, max_value = st.session_state.max_numftn,
                                on_change = lambda: set_update_data_flag(True),
                                key="slider_num_ftn_specific",
                                val0 = (1,5) )

    slider_num_ftn = st_utils.persist_widget(st.slider,"Total number of (any) potentially polymerizable functional groups",
                                min_value = 1, max_value = st.session_state.max_numftn,
                                on_change = lambda: set_update_data_flag(True),
                                key="slider_num_ftn",
                                val0 = (1,5) )
#st_utils.change_widget_fontsize(expander_label,"18px")
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
    try:
        molid, molidx = st.session_state["data_index"]
    except:
         molid, molidx = 0,0
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

# Visualize distribution
counts = []
multi_index, data_rxn = st.session_state["prev_data"], st.session_state["data_rxn"]
if multi_index.size == 0:
     st.markdown("##### No molecules left to analyze. Adjust filters.")
else:
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

# ===== States =====
st.session_state.state["reload_batch_evaluation"] = True