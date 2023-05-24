import streamlit as st
import numpy as np
import pandas as pd
import gspread_pdlite as gspdl
import os

# ===== Notes
# Global variables in session_state:
# < load/initialize once >
# - backend_sheet           output data worksheet
# - data_rxn
# - rxn_types
# - multi
# - base_data               (multi_index bare data, data_rxn)
# - max_MW                  from data_rxn
# - max_numftn              from data_rxn


# < storage >
# - eval_general
# - eval_mol                (can load from backend)
# - eval_details            (can load from backend)
# - prev_data               filtered data

# - proxy_indices           for batch evaluation
# - entry_comments_{ii}     molecule comments
# - entry_general_{ii}      overall molecule rating
# - entry_{ii}              rating for assigned reaction
# - entry_other_{ii}        rating for manually entered reaction
# - entry_other_name_{ii}   manually entered reaction name


# < state trackers >
# - state.reload_batch_evaluation

# - settings_initialized
# - b_update_data
# - data_index              (molid, matchidx)

# - b_update_data_batch         for batch evaluation
# - entry_other_reset       
# - entry_{ii}_load             whether to load
# - entry_general_{ii}_load     whether to try to load old data    
# - batch_page


# < settings (e.g. from widgets) >
# - userinfo
# - rxn_selection           (from a slider)
# - iteration_selection 
# - slider_MW
# - slider_num_ftn_specific
# - slider_num_ftn

# - mols_per_page           for batch evaluation
# - mols_per_page_widget    manual persistence, since need to cast to an integer and store and clear

# < dynamic >

# ===== Constants
homedir = os.path.dirname(__file__)

rxn_name_alias = {"simple":"acyclic alkene","rop-olefin":"ROMP"}
rxn_name_alias_reverse = {"acyclic alkene":"simple","ROMP":"rop-olefin"}

rating_scale = ("skip","1: bad","2","3: interesting","4","5: good")
rating_scale = ("skip (don't answer)","0: N/A","1: impossible","2: bad potential","3: workable potential","4: promising potential","5: probably works")
rating_scale = ("skip (don't answer)","incorrectly ID'd","1","2","3","4","5")
rating_scale_index = {entry:ix for ix,entry in enumerate(rating_scale)}

remove_step = False
remove_rxns = ["rop-thioether","rop-oxazoline","rop-phosphonite","rop-siloxane"]
read_local  = False

# ===== System state
st.session_state["state"] = {}
if "reload_batch_evaluation" not in st.session_state.state:
    st.session_state.state["reload_batch_evaluation"] = True


# ===== Functions

# ----- I/O
# For read-only, use google sheets with export:
# https://stackoverflow.com/questions/36096194/where-can-i-host-a-csv-so-i-can-directly-read-it-into-a-neo4j-database

@st.cache_resource(ttl=600)
def load_sheet():
    sheet = gspdl.open_sheet(st.secrets.data.google_key, st.secrets["gcp_service_account"])
    return sheet
    
@st.cache_data(ttl=600)
def load_data():
    if "base_data" in st.session_state:
        return st.session_state["base_data"][0], st.session_state["base_data"][0]

    # functional group data
    if read_local:
        filename = homedir + "/../data/" + "fg_analysis_2023-05-12b.csv"
        data_df = pd.read_csv(filename,index_col=False)
    else:
        url = gspdl.urlfy(st.secrets.data.fgroups_key)
        data_df = gspdl.url_to_df(url) #still slow (~30 seconds)

    #sheet = gspdl.open_sheet(st.secrets.data.fgroups_key,st.secrets["gcp_service_account"])
    #ws = sheet.get_worksheet(0)
    #data_df = gspdl.worksheet_to_df(ws)

    #def toset(mystr):
    #   return frozenset( ast.literal_eval(mystr))
    #data_df["matchidx"] = data_df["matchidx"].map(toset)
    multi = data_df.set_index(["molid","matchidx","rxn_name"])
    #if iterate_by_matchidx:
    #    multi   = data_df.set_index(["molid","rxn_name","ftn_id","matchid"])
    #else:
    #    multi = data_df.set_index(["molid","matchidx","rxn_name"])


    # functional group count data; makes initial filtering easier
    if read_local:
        filename = homedir + "/../data/" + "rxntypes_2023-05-12b.csv"
        data_rxn = pd.read_csv(filename,index_col=False)
    else:
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

def load_user_session():
    #sheet = gspdl.open_sheet(st.secrets.data.google_key, st.secrets["gcp_service_account"])
    sheet = st.session_state["backend_sheet"]

    ws = sheet.worksheet(st.secrets.data.name_mol)
    tmp_df = gspdl.worksheet_to_df(ws)
    st.session_state["eval_mol"] = tmp_df[ tmp_df.userinfo == st.session_state.userinfo ]
    
    ws = sheet.worksheet(st.secrets.data.name_details)
    tmp_df = gspdl.worksheet_to_df(ws)
    st.session_state["eval_details"] = tmp_df[ tmp_df.userinfo == st.session_state.userinfo ]


# ----- Analysis
def get_valid_reactions(data_rxn, bounds, rxn_name = None):
    """Returns molids of molecules that satisfy rxn_name
    Arguments:
        bounds (list, tuple): range to sit within
    """
    if rxn_name in rxn_name_alias_reverse:        
        rxn_name = rxn_name_alias_reverse[rxn_name]
        

    if rxn_name is None or "choose" in rxn_name:
        # In this case, count unique ftnl groups via multi_index data
        #multi_onlymolidmatchidx = multi.reset_index()[["molid","matchidx"]]
        #res = multi_onlymolidmatchidx.reset_index().set_index("molid").groupby(["molid"]).agg("nunique").matchidx
        #valid_molids = res.index.values
        valid_molids = data_rxn.loc[ 
            (data_rxn.num_ftn >= bounds[0]) & 
            (data_rxn.num_ftn <= bounds[1]) ].index.values
    elif "step" not in rxn_name:
        #valid_molids = data_rxn[ data_rxn[rxn_name] > 0 ].index.values
        valid_molids = data_rxn.loc[ 
            (data_rxn[rxn_name] >= bounds[0]) & 
            (data_rxn[rxn_name] <= bounds[1]) ].index.values
    else:
        # step reactions, need to count tuples that are not (0, 0) or '(0, 0)'
        tmp = data_rxn[rxn_name]
        import ast
        def transform(x):
            if isinstance(x,str):
                return ast.literal_eval(x)
            else:
                return x
        tmp = tmp.map(transform)

        def validate(x):
            # a bit ambiguous what it means for a step ftn, which has multiple ftnl groups, to satisfy the condition
            # here, we interpret as either ftnl group satisfying the criterion
            # e.g. for diels alder would want an & condition
            if x[0] >= bounds[0] and x[0] <= bounds[1] \
                or x[1] >= bounds[0] and x[1] <= bounds[1]:
                return True
            else:
                return False
        valid_molids = tmp[tmp.map(validate)].index.values
        #valid_molids = data_rxn[ (tmp != (0, 0)) & (tmp != "(0, 0)") ].index.values
        
    return valid_molids

def filter_rxn(data, data_rxn, rxn_name = None):
    """return molids that are consistent"""
    # TMP, remove step reactions
    rxn_names = data.index.get_level_values("rxn_name")
    step_rxn_names = [x for x in rxn_names if x.startswith("step")]

    # filter by reaction name. TODO: step reactions needs adjustment
    if rxn_name is None or rxn_name == "choose for me!":
        # TMP, remove step reactions
        subdata = data
        if len(remove_rxns) > 0:
            sub_data = subdata.iloc[ ~subdata.index.get_level_values("rxn_name").isin(remove_rxns) ]
        if remove_step:
            sub_data = subdata.iloc[ ~subdata.index.get_level_values("rxn_name").isin(step_rxn_names) ]
        sub_data = data
        inds = sub_data.index.get_level_values("molid").unique().values
    else:
        # filter by rxn_name
        if rxn_name in rxn_name_alias_reverse:
            rxn_name = rxn_name_alias_reverse[rxn_name]
        
        #inds = np.argwhere( (data_rxn[rxn_name]>0).values ).flatten() #alternative way to filter
        #inds = get_valid_reactions(data_rxn,(0,50),rxn_name)

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

        # filter by number of specified functional groups
        inds_where_num_ftn_specific = get_valid_reactions(
            data_rxn,
            (st.session_state.slider_num_ftn_specific[0],
            st.session_state.slider_num_ftn_specific[1]),
            rxn_name)
        sub_data = sub_data.query( "molid in @inds_where_num_ftn_specific")
    
    # return
    return sub_data


# ----- Processing
def update_filters():
    data_rxn = st.session_state["data_rxn"]
    multi = st.session_state["multi"]
    multi_filtered0 = filter_rxn(multi,data_rxn,None)
    if multi_filtered0.size == 0:
        multi_filtered0 = multi
    
    if st.session_state["b_update_data"]:
        tmp_multi_filtered = filter_rxn(multi,data_rxn,st.session_state.rxn_selection)

        if tmp_multi_filtered.size == 0:
            st.write("##### ERROR: Filter too strict, returning 0 molecules. Returning to previous data set.")
            #if st.session_state["prev_data"].size == 0:
            #    multi_filtered = multi_filtered0
            #else:
            #    multi_filtered = st.session_state["prev_data"]
            multi_filtered = tmp_multi_filtered
        else:
            multi_filtered = tmp_multi_filtered
            if multi_filtered.size == 0:
                st.write("##### ERROR: Filter too strict, returning 0 molecules. Returning to default data set.")
                #multi_filtered = filter_rxn(multi,data_rxn,None)
        st.session_state["prev_data"] = multi_filtered
    else:
        multi_filtered = st.session_state["prev_data"]
        if multi_filtered.size == 0:
            st.write("##### ERROR: Filter too strict, returning 0 molecules. Returning to default data set.")
            #multi_filtered = multi_filtered0
    st.session_state["b_update_data"] = False


def update_filters_legacy():
    """Old update filter used for batch_updates

    Returns:
        bool: whether or not to update proxy index

    Notes:
        Decided against fine-grained case-control, since generating proxy indices is pretty cheap anyway
    """
    update_index = False
    if st.session_state["b_update_data"]:
        tmp_multi_filtered = app_utils.filter_rxn(multi,data_rxn,st.session_state.rxn_selection)
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


# ===== Preprocessing
def initialize():
    if "b_update_data" not in st.session_state:
        st.session_state["b_update_data"] = True
    if "b_update_data_batch" not in st.session_state:
        st.session_state["b_update_data_batch"] = True
    if "b_update_data_single" not in st.session_state:
        st.session_state["b_update_data_single"] = True

    # backend data connection
    if "backend_sheet" not in st.session_state:
        st.session_state["backend_sheet"] = load_sheet()


    # user log that will be stored to backend
    if "eval_details" not in st.session_state:
        #store mol-level under rxn_name = "general"
        st.session_state["eval_details"] = pd.DataFrame(columns=["molid","molidx","rxn_name","smiles","userinfo","timestamp","rating"])
    if "eval_mol" not in st.session_state:
        st.session_state["eval_mol"] = pd.DataFrame(columns=["molid","molidx","smiles","userinfo","timestamp","comments_ftn","comments_mol","rating_mol"])
    if "eval_general" not in st.session_state:
        st.session_state["eval_general"] = pd.DataFrame(columns=["molid","molidx","smiles","userinfo","timestamp","comments_general"])


    # data_rxn
    if "data_rxn" not in st.session_state:
        if read_local:
            filename = homedir + "/../data/" + "rxntypes_2023-05-12b.csv"
            data_rxn = pd.read_csv(filename,index_col=False)
            st.session_state["data_rxn"] = data_rxn
        else:
            url = gspdl.urlfy(st.secrets.data.data_rxn_key)
            st.session_state["data_rxn"] = gspdl.url_to_df(url)
    data_rxn = st.session_state["data_rxn"]
    data_rxn

    # rxn_types
    if "rxn_types" not in st.session_state:
        rxn_types = data_rxn.columns.values
        # TMP: eliminate step reactions
        if remove_step:
            rxn_types = [x for x in rxn_types if not x.startswith("step")]
            #rxn_types = rxn_types[:-2]
        if len(remove_rxns) > 0:
            rxn_types = [x for x in rxn_types if x not in remove_rxns]
        rxn_types = rxn_types[:-3] #last two elements should be MW and rxn and smiles
        names_to_alias = []
        for ix,x in enumerate(rxn_types):
            if x in rxn_name_alias: #e.g. simple -> alkene-linear
                names_to_alias.append( (ix,x) )
        for ix,x in names_to_alias:
            rxn_types[ix] = rxn_name_alias[x]
        st.session_state["rxn_types"] = rxn_types
    rxn_types = st.session_state["rxn_types"]

    # some data_rxn metrics
    if "max_MW" not in st.session_state:
        st.session_state["max_MW"] = float(data_rxn.MW.max())
    if "max_numftn" not in st.session_state:
        st.session_state["max_numftn"] = int(data_rxn.num_ftn.max())

# molecule data
def first_load():
    if "base_data" not in st.session_state:
        multi,data_rxn = load_data()
        st.session_state["base_data"] = (multi,data_rxn)
        st.session_state["multi"] = multi
        st.session_state["data_rxn"] = data_rxn 
        #note that this data_rxn will have `num_ftn` and `MW` imputed.
        #in principle, should *already* be in data_rxn 

        # for batch evaluation, don't care about substituents. 
        # prune first (keep only one entry per unique multi index).
        # st.session_state["multi_rxn_unique"] = multi[ ~multi.index.duplicated(keep='first') ] 

    if "prev_data" not in st.session_state: #first filter
        if "rxn_selection" in st.session_state:
            multi_filtered = filter_rxn(multi,data_rxn,st.session_state.rxn_selection)
        else:
            multi_filtered = filter_rxn(multi,data_rxn,None)
        if multi_filtered.size == 0:
            multi_filtered = multi
        st.session_state["prev_data"] = multi_filtered


