import streamlit as st
import pandas as pd
import os
homedir = os.path.dirname(__file__)

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


st.markdown("# Welcome!")
if "b_update_data" not in st.session_state:
    st.session_state["b_update_data"] = False
def set_update_data_flag(flag):
    st.session_state["b_update_data"] = flag

if "max_MW" not in st.session_state:
    st.session_state["max_MW"] = data_rxn.MW.max()
if "max_numftn" not in st.session_state:
    st.session_state["max_numftn"] = int(data_rxn.num_ftn.max())

if "data_index" not in st.session_state:
    st.session_state["data_index"] = ("0","simple")

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
with st.form("general comments/bugs"):
    comment_area = st.text_area("General comments?","")
    submitted = st.form_submit_button("submit")
    