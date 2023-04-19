import streamlit as st
import pandas as pd
import os, time
homedir = os.path.dirname(__file__)

st.set_page_config(layout="wide")
st.markdown(
    """
    <style>
    [data-testid="stSidebar"][aria-expanded="true"]{
        max-width: 650px;
        min-width: 550px;
    }
    """,
    unsafe_allow_html=True,
)   

def persist_widget(widget,text,key,val0,action,*args,**kwargs):
    tmp_key = key + "_"
    # set up keys
    if key not in st.session_state:
        if tmp_key not in st.session_state:
            st.session_state[tmp_key] = val0
            st.session_state[key] = val0
        else:
            st.session_state[key] = st.session_state[tmp_key]
    else:
        st.session_state[tmp_key] = st.session_state[key]

    # set up callback() and return
    if widget in [st.button, st.download_button, st.form_submit_button]:
        def callback():
            st.session_state[key] = st.session_state[tmp_key]
            action()
        return widget(text,*args,key=tmp_key,on_click=callback,**kwargs)
    else:
        def callback():
            st.session_state[key] = st.session_state[tmp_key]
            action()
        return widget(text,*args,key=tmp_key,on_change=callback,**kwargs)
    
# ===== Setup =====
filename = homedir + "/../data/" + "rxntypes_2023-04-10.csv"
data_rxn = pd.read_csv(filename,index_col=False)
rxn_types = data_rxn.columns
# TMP: eliminate step reactions
rxn_types = [x for x in rxn_types if not x.startswith("step") and x != "smiles"]

if "b_update_data" not in st.session_state:
    st.session_state["b_update_data"] = False
def set_update_data_flag(flag):
    st.session_state["b_update_data"] = flag

if "max_MW" not in st.session_state:
    st.session_state["max_MW"] = data_rxn.MW.max()
if "max_numftn" not in st.session_state:
    st.session_state["max_numftn"] = int(data_rxn.num_ftn.max())

if "data_index" not in st.session_state:
    st.session_state["data_index"] = (0,"(0,)")

if "eval_details" not in st.session_state:
    #store mol-level under rxn_name = "general"
    st.session_state["eval_details"] = pd.DataFrame(columns=["molid","molidx","rxn_name","smiles","userinfo","timestamp","rating"])
if "eval_mol" not in st.session_state:
    st.session_state["eval_mol"] = pd.DataFrame(columns=["molid","molidx","smiles","userinfo","timestamp","comments_ftn","comments_mol","rating_mol"])
if "eval_general" not in st.session_state:
    st.session_state["eval_general"] = pd.DataFrame(columns=["molid","molidx","smiles","userinfo","timestamp","comments_general"])


# ===== ===== BEGIN LAYOUT ===== =====
st.markdown("# Welcome!")
user_info = persist_widget( st.text_input, "email",
                               key = "userinfo", val0="",
                               action = lambda: None)

# ===== Navigation =====
#st.markdown("## Navigation & Settings")
#st.markdown("### On next molecule, show...")

rxn_selection = persist_widget( st.selectbox, "reaction type",
                               key = "rxn_selection", val0 = "choose for me!",
                               options = ("choose for me!",*rxn_types),
                               action = lambda: set_update_data_flag(True))

iteration_selection = persist_widget( st.selectbox,
                                     "molecule iteration mode:",
                                     options = ("random","sequential"),
                                     key="iteration_selection",
                                     val0 = "sequential",
                                     action = lambda: set_update_data_flag(True))

# MW
slider_MW = persist_widget(st.slider, "MW range",
                            min_value = 0., max_value = st.session_state.max_MW,
                            key = "slider_MW",
                            val0 = (10.,st.session_state.max_MW),
                            action= lambda: set_update_data_flag(True) )

# Simplicity/Complexity
# (# of polymerizations identified, # functional groups, # subsitutents)
slider_num_ftn = persist_widget(st.slider,"number functional groups",
                            min_value = 1, max_value = st.session_state.max_numftn,
                            action = lambda: set_update_data_flag(True),
                            key="slider_num_ftn",
                            val0 = (1,st.session_state.max_numftn) )

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
    comment_dict["smiles"] = data_rxn.loc[molid].smiles
    comment_dict["userinfo"] = st.session_state["userinfo"]
    comment_dict["timestamp"] = pd.Timestamp(time.time(),unit="s")
    comment_dict["comments_general"] = st.session_state["comments_general"]

    st.session_state["eval_general"] = pd.concat([ st.session_state["eval_general"], 
                                                    pd.DataFrame([comment_dict]) ], ignore_index=True)

with st.form("general comments/bugs",clear_on_submit=True):
    comment_area = st.text_area("General comments?","",key="comments_general")

    submitted = st.form_submit_button("submit",on_click=log_general_comment)
    
    
