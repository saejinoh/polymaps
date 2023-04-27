import streamlit as st
import pandas as pd
import os, time
from rdkit import Chem
from rdkit.Chem import Descriptors

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


rxn_types = data_rxn.columns.values
# TMP: eliminate step reactions
rxn_types = [x for x in rxn_types if not x.startswith("step") and x != "smiles"]
rxn_types = rxn_types[:-3] #last two elements should be MW and rxn
rxn_name_alias = {"simple":"alkene-linear"}
rxn_name_alias_reverse = {"alkene-linear":"simple"}
names_to_alias = []
for ix,x in enumerate(rxn_types):
    if x in rxn_name_alias: #e.g. simple -> alkene-linear
        names_to_alias.append( (ix,x) )
for ix,x in names_to_alias:
    rxn_types[ix] = rxn_name_alias[x]

if "b_update_data" not in st.session_state:
    st.session_state["b_update_data"] = False
def set_update_data_flag(flag):
    st.session_state["b_update_data"] = flag

if "max_MW" not in st.session_state:
    st.session_state["max_MW"] = float(data_rxn.MW.max())
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
st.markdown("- On this page, you can choose some basic filters for the molecules to view, as well as leave general comments about the process, monomers, etc.  \n")
st.markdown("- On the `Evaluation` tab (in the sidebar), you will be shown monomers, with one highlighted functional group at a time, and the opportunity to rate the suitability of the functional groups for polymerization. \n  - Rating is on a scale of `1 (bad)` - `3 (interesting)` - `5 (good)`. Default is 0 to skip judging. \n  - The first load may take a few minutes!* \n")
st.markdown("- On the `Results` tab, you can view your top-rated molecules, and optionally download a `.csv` of all molecules that you've rated.")

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


# Preload data
@st.cache_data(ttl=300)
def load_data():
    if "base_data" in st.session_state:
        return st.session_state["base_data"][0], st.session_state["base_data"][0]

    # functional group data
    #filename = homedir + "/../../data/" + "fg_analysis_2023-04-21.csv"
    #data_df = pd.read_csv(filename,index_col=False)
    url = gspdl.urlfy(st.secrets.data.fgroups_key)
    data_df = gspdl.url_to_df(url) #still slow (~30 seconds)
    #sheet = gspdl.open_sheet(st.secrets.data.fgroups_key,st.secrets["gcp_service_account"])
    #ws = sheet.get_worksheet(0)
    #data_df = gspdl.worksheet_to_df(ws)

    #def toset(mystr):
    #   return frozenset( ast.literal_eval(mystr))
    #data_df["matchidx"] = data_df["matchidx"].map(toset)
    multi = data_df.set_index(["molid","matchidx","rxn_name"])

    # functional group count data; makes initial filtering easier
    #filename = homedir + "/../../data/" + "rxntypes_2023-04-21.csv"
    #data_rxn = pd.read_csv(filename,index_col=False)
    url = gspdl.urlfy(st.secrets.data.data_rxn_key)
    data_rxn = gspdl.url_to_df(url)

    # evaluations
    # evaluate #ftnl groups identified
    if "num_ftn" not in data_rxn.columns:
        data_rxn["num_ftn"] = 0.

        num_ftn = multi.reset_index().set_index("molid").groupby(["molid"]).agg("nunique").matchidx
        data_rxn.loc[multi.index.unique("molid"),"num_ftn"] = num_ftn

    # evaluate MW
    if "MW" not in data_rxn.columns:
        data_rxn["MW"] = 0.
        def get_MW(row):
            mol = Chem.MolFromSmiles(row.smiles)
            row.MW = Descriptors.MolWt(mol)
            return row
        data_rxn = data_rxn.apply(get_MW,axis=1)
    
    return multi,data_rxn

if "base_data" not in st.session_state:
    multi,data_rxn = load_data()
    st.session_state["base_data"] = (multi,data_rxn)