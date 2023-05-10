#Main imports
import streamlit as st
import time
import os
import ast
import copy
homedir = os.path.dirname(__file__)

st.set_page_config(layout="wide")

# database connection
import gspread_pdlite as gspdl
@st.cache_resource
def load_sheet():
    sheet = gspdl.open_sheet(st.secrets.data.google_key, st.secrets["gcp_service_account"])
    return sheet
sheet = load_sheet()

# Analysis imports
import random
import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.Chem import rdChemReactions
from rdkit.Chem import Draw
import mychem
iterate_by_matchidx = False

# ===== Loading
@st.cache_data(ttl=600)
def load_data():
    if "base_data" in st.session_state:
        return st.session_state["base_data"][0], st.session_state["base_data"][0]

    # functional group data
    #filename = homedir + "/../../data/" + "fg_analysis_2023-04-21.csv"
    url = gspdl.urlfy(st.secrets.data.fgroups_key)
    data_df = gspdl.url_to_df(url)

    if iterate_by_matchidx:
        multi   = data_df.set_index(["molid","rxn_name","ftn_id","matchid"])
    else:
        multi = data_df.set_index(["molid","matchidx","rxn_name"])

    url = gspdl.urlfy(st.secrets.data.data_rxn_key)
    data_rxn = gspdl.url_to_df(url)

    # evaluations
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
else:
    multi,data_rxn = st.session_state["base_data"]

# ===== Chemical Utilities
rxn_name_alias = {"simple":"alkene-linear"}
rxn_name_alias_reverse = {"alkene-linear":"simple"}


def filter_rxn(data, data_rxn, rxn_name = None):
    """return molids that are consistent"""
    # TMP, remove step reactions
    rxn_names = data.index.get_level_values("rxn_name")
    step_rxn_names = [x for x in rxn_names if x.startswith("step")]

    # filter by reaction name. TODO: step reactions needs adjustment
    if rxn_name is None or rxn_name == "choose for me!":
        # TMP, remove step reactions
        sub_data = data.iloc[ ~data.index.get_level_values("rxn_name").isin(step_rxn_names) ]
        sub_data = data
        inds = sub_data.index.get_level_values("molid").unique().values
    else:
        # filter by rxn_name
        if rxn_name in rxn_name_alias:
            rxn_name = rxn_name_alias_reverse[rxn_name]

        inds = np.argwhere( (data_rxn[rxn_name]>0).values ).flatten() #alternative way to filter
        #sub_data = data.query("molid in @inds")
        rxns_of_interest = [rxn_name]
        sub_data = data.query("rxn_name in @rxns_of_interest")
        #sub_data = sub_data[ sub_data.index.get_level_values(rxn_name).isin([rxn_name])] 

    if "prev_data" in st.session_state and "slider_MW" in st.session_state:
        # filter by MW
        inds_where_MW_range = data_rxn.loc[ 
            (data_rxn.MW>= st.session_state.slider_MW[0]) & 
            (data_rxn.MW <= st.session_state.slider_MW[1]) ].index.values

        sub_data = sub_data.query( "molid in @inds_where_MW_range")

        # filter by number of functional groups
        inds_where_num_ftn_range = data_rxn.loc[ 
            (data_rxn.num_ftn>= st.session_state.slider_num_ftn[0]) & 
            (data_rxn.num_ftn <= st.session_state.slider_num_ftn[1]) ].index.values
        sub_data = sub_data.query( "molid in @inds_where_num_ftn_range")

    # return
    return sub_data

# ===== Load Data
if "settings_initialized" not in st.session_state:
    st.markdown("# Upon browser refresh, please revisit Settings page first.")


if "prev_data" not in st.session_state: #first time through
    if "rxn_selection" in st.session_state:
        multi_filtered = filter_rxn(multi,data_rxn,st.session_state.rxn_selection)
    else:
        multi_filtered = filter_rxn(multi,data_rxn,None)
    if multi_filtered.size == 0:
        multi_filtered = multi
    st.session_state["prev_data"] = multi_filtered

    #if iterate_by_matchidx:
    #    st.session_state["data_index"] = generate_index_by_matchid(multi_filtered)
    #else:
    #    st.session_state["data_index"] = generate_index(multi_filtered)
else:
    multi_filtered = st.session_state["prev_data"]
molids = multi_filtered.index.get_level_values("molid").unique()
molnum = molids.values.size

multi_filtered0 = filter_rxn(multi,data_rxn,None) #default data set if filters fail


if "b_update_data_batch" not in st.session_state:
    st.session_state["b_update_data_batch"] = False
def set_update_data_flag(flag):
    st.session_state["b_update_data_batch"] = flag


def update_filters():
    update_index = False
    if st.session_state["b_update_data_batch"]:
        tmp_multi_filtered = filter_rxn(multi,data_rxn,st.session_state.rxn_selection)
        if multi.size == 0:
            st.write("##### ERROR: Filter too strict, returning 0 molecules. Returning to previous data set.")
            if st.session_state["prev_data"].size == 0:
                multi_filtered = multi_filtered0
                update_index = True
            else:
                multi_filtered = st.session_state["prev_data"]
                update_index = False
        else:
            multi_filtered = tmp_multi_filtered
            if multi_filtered.size == 0:
                st.write("##### ERROR: Filter too strict, returning 0 molecules. Returning to default data set.")
                multi_filtered = multi_filtered0
            update_index = True
        st.session_state["b_update_data_batch"] = False
        st.session_state["prev_data"] = multi_filtered
    else:
        multi_filtered = st.session_state["prev_data"]
        update_index = False
        if multi_filtered.size == 0:
            st.write("##### ERROR: Filter too strict, returning 0 molecules. Returning to default data set.")
            multi_filtered = multi_filtered0
            update_index = True
    set_update_data_flag(False)
    return update_index

# prune multi_filtered (keep only one entry per unique multi index)
multi_filtered = multi_filtered[ ~multi_filtered.index.duplicated(keep='first') ]


# update proxy index
update_index = update_filters()

if update_index == True \
    or "proxy_indices" not in st.session_state \
    or "prev_iteration_mode" not in st.session_state \
    or st.session_state.prev_iteration_mode != st.session_state.iteration_selection \
    or len(st.session_state.proxy_indices) != multi_filtered.shape[0]:
    #the above if statements test for:
    # 0) data has been updated for filters
    # 1) proxy_indices have been initialized
    # 2) iteration_mode has been set
    # 3) iteration_mode has been toggled
    # 4) light check that data has been resized/re-filtered by single-mol evaluation
    # mostly, should rely on `b_update_data_batch` to know if data has been updated or not
    # i.e. don't have to know about the other page!

    proxy_indices = np.arange(multi_filtered.shape[0])
    if st.session_state.iteration_selection == "random":
        np.random.shuffle(proxy_indices)
    st.session_state["proxy_indices"] = proxy_indices
else:
    proxy_indices = st.session_state["proxy_indices"]


if "batch_page" not in st.session_state:
    st.session_state["batch_page"] = 1 #change to 0

if "prev_iteration_mode" not in st.session_state \
    or st.session_state.prev_iteration_mode != st.session_state.iteration_selection:
    st.session_state["batch_page"] = 1 #change to 0

st.session_state["prev_iteration_mode"] = st.session_state.iteration_selection


# ===== Page number
if "mols_per_page" not in st.session_state:
    st.session_state["mols_per_page"] = 24
mols_per_page = st.session_state.mols_per_page #TODO: need to persist this widget...!
last_page_index = int(np.floor( len(proxy_indices)/mols_per_page ))

# retrieve block indices
def get_page_indices(page):
    """Get begin & end proxy indices for a page

    Args:
        page (int): page number

    Note:
        acts on `proxy_indices` array and last_page_index defined globally
        idea is to use .iloc(ind0,ind1) to retrieve entries
    """
    ind0 = min( (page-1) * mols_per_page, last_page_index*mols_per_page )
    ind1 = min( ind0 + mols_per_page, len(proxy_indices) )

    return ind0, ind1

ind0,ind1 = get_page_indices( st.session_state.batch_page )


# generate images of molecules
data_slice = multi_filtered.iloc[proxy_indices[ind0:ind1]]
data_slice
data_rxn.iloc[data_slice.index.unique("molid")]

mols = []
svgs = []
molids = []
rxn_names = []
for row in data_slice.iterrows():
    mols.append( Chem.MolFromSmiles(row[1].smiles) )
    molids.append( row[0][0] )
    rxn_names.append( row[0][2] )
    ftn_group_ids = ast.literal_eval(row[0][1]) #depends on ordering of index; should get ids of atoms in ftnl group.
    
    highlightcolors,highlightbonds = mychem.color_ftn(mols[-1],ftn_group_ids)
    im = mychem.highlight_draw(mols[-1],highlightcolors,highlightbonds,format="svg",
                               wd=250,ht=250)
    svgs.append(im)

# ===== Display
with st.sidebar:
    def get_closest_page():
        st.session_state.mols_per_page = st.session_state.mols_per_page_widget
        st.session_state["batch_page"] = 1 + int( np.floor( ind0/st.session_state.mols_per_page ) )

    st.select_slider("**Results per page**",(12,24,48,96),value = st.session_state.mols_per_page,
                     key="mols_per_page_widget",on_change=get_closest_page)
    
    def store_page():
        st.session_state.batch_page = int(st.session_state.batch_page_text)
        st.session_state.batch_page_text = ""
    st.text_input("**Jump to page**",value = "", key="batch_page_text",
                  on_change=store_page)
    
    def page_forward():
        st.session_state.batch_page = min(last_page_index+1, st.session_state.batch_page + 1)
    def page_backward():
        st.session_state.batch_page = max(1, st.session_state.batch_page - 1)
    st.button("next page",on_click = page_forward)
    st.button("previous page",on_click = page_backward)

# Main Text
st.markdown(f"## Page `{st.session_state.batch_page}`/`{last_page_index+1}`")
#st.markdown(f"- **results {ind0+1}~{ind1}**")
#st.write( proxy_indices[ind0:ind1] )

# Entry Area
rating_scale = ("skip","1: bad","2","3: interesting","4","5: good")
n_rows = int(mols_per_page/3)

entries = []
entries_general = []
for ia in range(n_rows):
    with st.container():
        cols = st.columns(3)
        for ic,col in enumerate(cols):
            index_abs = ia*3 + ic
            with col:
                if index_abs < len(svgs):
                    st.image( svgs[index_abs] )
                    st.write(f"case `{index_abs + ind0}`, **mol ID: `{molids[index_abs]}`**")
                    entries.append( st.selectbox(f"**quality for `{rxn_names[index_abs]}` polymerization**",
                                                rating_scale,
                                                key=f"entry_{index_abs}"
                                                ))
                    entries_general.append( st.selectbox("**overall monomer quality**",
                                                rating_scale,
                                                key=f"entry_general_{index_abs}"
                                                ))
                #data_slice.iloc[ia*3 + ic]
