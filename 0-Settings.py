import streamlit as st
import pandas as pd
import os, time
homedir = os.path.dirname(__file__)

st.set_page_config(layout="wide")
st.session_state["settings_initialized"] = True
#st.secrets["gcp_service_account"]

# database connection
import gspread_pdlite as gspdl
@st.cache_resource
def load_sheet():
    sheet = gspdl.open_sheet(st.secrets.data.google_key, st.secrets["gcp_service_account"])
    return sheet
sheet = load_sheet()

# utilities
def persist_widget(widget,*args,key=None,val0=None,action=lambda: None,**kwargs):
    """Persist a widget's state by using a dummy key and a callback to ensure syncing of values.

    Args:
        widget (_type_): _description_
        key (_type_): _description_
        val0 (_type_): _description_
        action (_type_): _description_

    Returns:
        _type_: _description_
    
    Notes:
        syntax is minimally invasive:
        ```
        persist_widget(widget, *(args as would be usually used), **(kwargs) as would be usually used )
        ```
        in addition, must provide the mandatory `key` and `val0` arguments explicitly by name instead of positionally.
    """
    if key is None:
        raise ValueError("must provide a key to persist the widget")
    if val0 is None:
        raise ValueError("must provide initial value of widget, via val0 argument")

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
        # technically shouldn't ever be persisting a button...
        raise ValueError("Does it really make sense to persist a button value?")
        #def callback():
        #    st.session_state[key] = st.session_state[tmp_key]
        #    action()
        #return widget(*args,key=tmp_key,on_click=callback,**kwargs)
    else: #should be `on_change``
        action = kwargs["on_change"]
        def callback():
            st.session_state[key] = st.session_state[tmp_key]
            action()
        kwargs["on_change"] = callback
        return widget(*args,key=tmp_key,**kwargs)
    
# ===== Setup =====
#filename = homedir + "/../data/" + "rxntypes_2023-04-10.csv"
#data_rxn = pd.read_csv(filename,index_col=False)
url = gspdl.urlfy(st.secrets.data.data_rxn_key)
data_rxn = gspdl.url_to_df(url)

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
                               on_change = lambda: None)

# ===== Navigation =====
#st.markdown("## Navigation & Settings")
#st.markdown("### On next molecule, show...")

rxn_selection = persist_widget( st.selectbox, "reaction type",
                               key = "rxn_selection", val0 = "choose for me!",
                               options = ("choose for me!",*rxn_types),
                               on_change = lambda: set_update_data_flag(True))

iteration_selection = persist_widget( st.selectbox,
                                     "molecule iteration mode:",
                                     options = ("random","sequential"),
                                     key="iteration_selection",
                                     val0 = "sequential",
                                     on_change = lambda: set_update_data_flag(True))

# MW
slider_MW = persist_widget(st.slider, "MW range",
                            min_value = 0., max_value = st.session_state.max_MW,
                            key = "slider_MW",
                            val0 = (10.,st.session_state.max_MW),
                            on_change = lambda: set_update_data_flag(True) )

# Simplicity/Complexity
# (# of polymerizations identified, # functional groups, # subsitutents)
slider_num_ftn = persist_widget(st.slider,"number functional groups",
                            min_value = 1, max_value = st.session_state.max_numftn,
                            on_change = lambda: set_update_data_flag(True),
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
    
    # Save
    st.session_state["eval_general"].to_csv("eval_general.csv",index=False)
    st.write(pd.Series(comment_dict))
    st.write(pd.DataFrame([comment_dict]))
    ws    = sheet.worksheet( st.secrets.data.name_general )
    gspdl.worksheet_append( ws, pd.Series(comment_dict) )

if "userinfo" not in st.session_state \
    or st.session_state["userinfo"] in ["",None] \
    or "@" not in st.session_state["userinfo"]:
    st.markdown("##### Please enter a valid e-mail above in order to submit general comments")
else:
    with st.form("general comments/bugs",clear_on_submit=True):
        comment_area = st.text_area("General comments?","",key="comments_general")
        submitted = st.form_submit_button("submit",on_click=log_general_comment)
    
    
