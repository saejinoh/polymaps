# Main imports
import streamlit as st
import time
import os
import ast
import copy
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

# database connection
import gspread_pdlite as gspdl
@st.cache_resource
def load_sheet():
    sheet = gspdl.open_sheet(st.secrets.data.google_key, st.secrets["gcp_service_account"])
    return sheet
sheet = load_sheet()

# For read-only, use google sheets with export:
# https://stackoverflow.com/questions/36096194/where-can-i-host-a-csv-so-i-can-directly-read-it-into-a-neo4j-database


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


@st.cache_data(ttl=600)
def load_data():
    if "base_data" in st.session_state:
        return st.session_state["base_data"][0], st.session_state["base_data"][0]

    # functional group data
    #filename = homedir + "/../../data/" + "fg_analysis_2023-04-21.csv"
    #data_df = pd.read_csv(filename,index_col=False)
    url = gspdl.urlfy(st.secrets.data.fgroups_key)
    #st.write(url)
    data_df = gspdl.url_to_df(url)
    #sheet = gspdl.open_sheet(st.secrets.data.fgroups_key,st.secrets["gcp_service_account"])
    #ws = sheet.get_worksheet(0)
    #st.write(ws)
    #data_df = gspdl.worksheet_to_df(ws)

    #def toset(mystr):
    #   return frozenset( ast.literal_eval(mystr))
    #data_df["matchidx"] = data_df["matchidx"].map(toset)
    if iterate_by_matchidx:
        multi   = data_df.set_index(["molid","rxn_name","ftn_id","matchid"])
    else:
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
else:
    multi,data_rxn = st.session_state["base_data"]

# Key variable initializations
if "max_MW" not in st.session_state:
    st.session_state["max_MW"] = data_rxn.MW.max()
if "max_numftn" not in st.session_state:
    st.session_state["max_numftn"] = int(data_rxn.num_ftn.max())

quality_dict = {0:"poor",1:"decent",2:"promising"}
eprop_dict = {-1:"withdrawing", 1:"donating", 0:"undetermined"}
index = {"molid":0, "rxn_name":1, "ftn_id":2, "matchid":3}
#rxn_types = data_rxn.columns
# TMP: eliminate step reactions
#rxn_types = [x for x in rxn_types if not x.startswith("step") and x != "smiles"]

if "b_update_data" not in st.session_state:
    st.session_state["b_update_data"] = False
def set_update_data_flag(flag):
    st.session_state["b_update_data"] = flag


# ===== Chemical Utilities
rxn_name_alias = {"simple":"alkene-linear"}
rxn_name_alias_reverse = {"alkene-linear":"simple"}
def random_molecule(data,rxn_name = None):
    """return a random molecule from the set"""
    molids = data.index.levels[0].unique().values
    num_mols = len(molids)

    # choose random mol
    molid = np.random.randint(0,num_mols)

    subdata = multi.loc[(molid,slice(None),slice(None),slice(None))]

    return molid, num_mols, subdata

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

def generate_index(multi_filtered,trial_molid=None):
    """Generate next molecule's index (and data) if needed

    Args:
        multi (_type_): _description_
        data_rxn (_type_): _description_
        rxn_selection (_type_): _description_

    Returns:
        tuple: (index tuple), filtered_data

    Notes:
        To simplify interface, makes use of st.session_state:
            - ftn_tracking
            - prev_data
        And also saves to st.session_state, both the above, and:
            - data_index
            - data_index_description
        (may refactor to clean up in the future)

        Global variables used:
            user inputs
    """
    
    if "initialized" not in st.session_state: #Only the first run
        molids = multi_filtered.index.get_level_values("molid").unique()
        molnum = molids.values.size
        mol_subid = np.random.randint(0,molnum) #default to random molecule
        molid = molids[mol_subid]

        mol_specific_data = multi_filtered.query("molid == @molid")
        num_ftnl_groups = mol_specific_data.reset_index().agg("nunique").matchidx
        ftn_subid = 0 #default init
        ftn_group_ids = mol_specific_data.index.get_level_values("matchidx").unique()[ftn_subid]
        
        description_base = ""
        st.session_state["initialized"] = True
    elif trial_molid is not None:
        #multi_filtered = st.session_state.prev_data

        if trial_molid not in multi_filtered.index.get_level_values("molid").unique().values:
            st.markdown(f"##### Input molid `{st.session_state.molid_input}` not in current filters. Ignoring.")
            molid = trial_molid
            ftn_subid = st.session_state["ftn_tracking"][0]
            description_base = ""
        else: 
            molid = trial_molid
            molids = multi_filtered.index.get_level_values("molid").unique()
            molnum = molids.values.size
            mol_subid = np.argwhere( molids == molid )[0][0]

            ftn_subid = 0 #default init
            
            description_base = f"##### Manually selected molecule `{trial_molid}`  \n"
    else:
        #prev_mol_subid = st.session_state["mol_subid"]
        prev_molid = st.session_state["data_index"][0]
        prev_ftn_subid, prev_num_ftnl_groups = st.session_state["ftn_tracking"]


        if prev_molid not in multi_filtered.index.get_level_values("molid").unique():
            b_mol_out_of_set = True
        else:
            b_mol_out_of_set = False

        if (prev_ftn_subid >= prev_num_ftnl_groups - 1):
            b_finished_mol = True
        else:
            b_finished_mol = False

        if b_mol_out_of_set:
            molids = multi_filtered.index.get_level_values("molid").unique()
            molnum = molids.values.size

            description_base = "##### previous molecule out of filtered set, returning random molecule  \n"
            mol_subid = np.random.randint(0,molnum) 
            molid = molids[mol_subid]
            ftn_subid = 0 #always default to first ftnl group
        elif b_finished_mol:
            molids = multi_filtered.index.get_level_values("molid").unique()
            molnum = molids.values.size

            description_base = "##### new molecule  \n"
            current_mol_subid = np.argwhere( molids == prev_molid )[0][0]
            molid = prev_molid

            reached_last_mol = True if (current_mol_subid >= molnum - 1) else False

            if st.session_state["iteration_selection"] == "random":
                mol_subid = np.random.randint(0,molnum) #IMPORTANT
                molid = molids[mol_subid]
            elif st.session_state["iteration_selection"] == "sequential":
                if reached_last_mol:
                    mol_subid = 0
                    molid = molids[mol_subid]
                    description_base += "previous molecule was first in set, returning to first molecule  \n"
                else:
                    mol_subid = current_mol_subid + 1
                    molid = molids[mol_subid]
                
            ftn_subid = 0 #always default to first ftnl group
        else: #only generate new ftn index
            multi_filtered = st.session_state["prev_data"]
            molids = multi_filtered.index.get_level_values("molid").unique()
            molnum = molids.values.size
            mol_subid = np.argwhere( molids == prev_molid )[0][0] #prev_mol_subid
            molid = molids[mol_subid]
            description_base = "##### continuing with molecule  \n"

            ftn_subid = prev_ftn_subid + 1
            
    # Common actions
    mol_specific_data = multi_filtered.query("molid == @molid")
    num_ftnl_groups = mol_specific_data.reset_index().agg("nunique").matchidx
    ftn_group_ids = mol_specific_data.index.get_level_values("matchidx").unique()[ftn_subid]
    
    # Save and return
    numdigits=len(str(molnum))
    numdigits_str = "0"+str(numdigits)

    if "rxn_selection" not in st.session_state: #to avoid circularity during first call before navigation bar is loaded
        description = description_base + f"##### **Molecule ID:**\t\t `{molid:{numdigits_str}}` ({mol_subid+1}/{molnum} monomers identified for the chosen reaction type)  \n##### **Showing:**\t\t potential functional group `{ftn_subid+1}`/`{num_ftnl_groups}` for rxn type `any`"
    else:
        description = description_base + f"##### **Molecule ID:**\t\t `{molid:{numdigits_str}}`   \n({mol_subid+1}/{molnum} monomers identified for the chosen reaction type)  \n##### **Potential functional group:**\t\t `{ftn_subid+1}`/`{num_ftnl_groups}`  \nfor rxn type `{st.session_state.rxn_selection}`"
    st.session_state["prev_data"] = multi_filtered
    st.session_state["data_index"] = (molid, ftn_group_ids)
    st.session_state["ftn_tracking"] = (ftn_subid,num_ftnl_groups)
    #st.session_state["mol_subid"] = mol_subid
    st.session_state["data_index_description"] = description

    return molid,ftn_group_ids
    
def generate_index_by_matchid(multi_filtered):
    """Right now just chooses random reaction groups. Incomplete. Here for legacy reasons."""

    # utility
    molids = multi.index.get_level_values("molid").unique()
    molnum = molids.values.size
    numdigits=len(str(molnum))
    numdigits_str = "0"+str(numdigits)

    # generate
    if "data_index" in st.session_state:
        prev_molid,prev_rxn_name,prev_ftnid,prev_matchid = st.session_state["data_index"]
        first_time = False
    else: #only the first time through
        prev_molid = -1
        first_time = True

    subid = np.random.randint(0,molnum) #IMPORTANT
    molid = molids[subid]

    mol_specific_data = multi.query(f"molid == {molid}")

    # if "choose for me!", choose a random reaction
    if st.session_state.rxn_selection == "choose for me!":
        rxn_types = mol_specific_data.index.get_level_values("rxn_name").unique().values
        rxn_name = random.choice(rxn_types)
    else:
        rxn_name = mol_specific_data.index[0][1] # in my scheme, should all be the same reaction name???

    #if not rxn_name.startswith("step"): # TMP filter out step
    ftn_id = 0
    match_totals = len(mol_specific_data.index.get_level_values("matchid").unique().values)
    match_id = np.random.randint(0,match_totals) #chooses a random match

    # store and return
    st.session_state["data_index_description"] = f"**Molecule:**\t\t `{molid:{numdigits_str}}` ({subid}/{molnum} monomers identified for the chosen reaction type)  \n**Rxn type:**\t\t `{rxn_name}`  \n**Showing:**\t\t `{match_id+1}`/`{match_totals}` reactive sites identified"

    return molid,rxn_name,ftn_id,match_id

def characterize_substituents(match_specific_data):
    """_summary_

    Args:
        match_specific_data (_type_): _description_

    Returns:
        _type_: _description_

    Notes:
        match_specific_data *should* only have 1 reaction name

        Which can be done with:
        match_specific_data = match_specific_data.loc[ match_specific_data.index.values[0] ]

        Have to be careful to turn Series back into multiindex
    """
    ftn_group_ids = ast.literal_eval(match_specific_data.index.get_level_values("matchidx").values[0])

    substituents = []
    for row in match_specific_data.iterrows():
        tmp_index, tmp_data = row
        substituents.append( (tmp_data.subid, eprop_dict[tmp_data.electronics], tmp_data.bulkiness))

    n_substituents = len(substituents)
    evaluation = f"Functional group (colored in `green`) w/ atoms `{ftn_group_ids}`\n"
    evaluation +=f"    with `{n_substituents}` substituents.\n\n"
    description = ""
    for entry in substituents:
        description += f"- root atom `{entry[0]}`\n  - electronic properties: `{entry[1]}`\n  - bulkiness: `{entry[2]:0.2}`\n"
    return evaluation, description

def characterize_substituents_by_matchid(match_specific_data):
    substituents = []
    for row in match_specific_data.iterrows():
        tmp_index, tmp_data = row
        substituents.append( (tmp_data.subid, eprop_dict[tmp_data.electronics], tmp_data.bulkiness))

    ftn_group_ids = ast.literal_eval(match_specific_data.iloc[0].matchidx)

    n_substituents = len(substituents)
    evaluation = f"ID'd functional group (colored in `green`) w/ atoms `{ftn_group_ids}`\n"
    evaluation +=f"    with `{n_substituents}` substituents.\n\n"
    description = ""
    for entry in substituents:
        description += f"- root atom `{entry[0]}`\n  - electronic properties: `{entry[1]}`\n  - bulkiness: `{entry[2]:0.2}`\n"
    return evaluation,description


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

    if iterate_by_matchidx:
        st.session_state["data_index"] = generate_index_by_matchid(multi_filtered)
        
    else:
        st.session_state["data_index"] = generate_index(multi_filtered)
else:
    multi_filtered = st.session_state["prev_data"]
molids = multi_filtered.index.get_level_values("molid").unique()
molnum = molids.values.size

multi_filtered0 = filter_rxn(multi,data_rxn,None)


# ===== Sidebar and submission logic
def update_filters():
    if st.session_state["b_update_data"]:
        tmp_multi_filtered = filter_rxn(multi,data_rxn,st.session_state.rxn_selection)
        if multi.size == 0:
            st.write("##### ERROR: Filter too strict, returning 0 molecules. Returning to previous data set.")
            if st.session_state["prev_data"].size == 0:
                multi_filtered = multi_filtered0
            else:
                multi_filtered = st.session_state["prev_data"]
        else:
            multi_filtered = tmp_multi_filtered
            if multi_filtered.size == 0:
                st.write("##### ERROR: Filter too strict, returning 0 molecules. Returning to previous data set.")
                multi_filtered = multi_filtered0
        st.session_state["b_update_data"] = False
        st.session_state["prev_data"] = multi_filtered
    else:
        multi_filtered = st.session_state["prev_data"]
        if multi_filtered.size == 0:
            st.write("##### ERROR: Filter too strict, returning 0 molecules. Returning to previous data set.")
            multi_filtered = multi_filtered0
    set_update_data_flag(False)


def log_evaluation():
    """Log user evaluation of molecule

    Notes:
        Requires certain session_state variables to be defined already:
        - userinfo
        - 
    """
    #"molid","molidx","smiles","userinfo","timestamp","comments_ftn","comments_mol","rating_mol"
    #"molid","molidx","smiles","userinfo","timestamp","rxn_name","rating"

    root_dict = {}
    molid, molidx = st.session_state["data_index"]
    root_dict["molid"],root_dict["molidx"] = molid, molidx
    root_dict["smiles"] = data_rxn.loc[molid].smiles
    root_dict["userinfo"] = st.session_state["userinfo"]
    root_dict["timestamp"] = pd.Timestamp(time.time(),unit="s")

    # comments
    comment_dict = copy.copy(root_dict)
    comment_dict["comments_ftn"] = st.session_state["comments_ftn"]
    comment_dict["comments_mol"] = st.session_state["comments_mol"]
    comment_dict["rating_mol"] = st.session_state["rating_mol"]
    comment_dict_df = pd.DataFrame([comment_dict])
    st.session_state["eval_mol"] = pd.concat([ st.session_state["eval_mol"], 
                                                 comment_dict_df],
                                                ignore_index=True)

    # rxn-specific data
    rxn_ratings = []
    for rxn_name_entry_in_session_state in st.session_state["rxns_for_this_ftn_group"]:
        rxn_rating_dict = copy.copy(root_dict)
        rxn_name = rxn_name_entry_in_session_state.removeprefix("rating_")
        rxn_rating_dict["rxn_name"] = rxn_name
        rxn_rating_dict["rating"] = st.session_state[rxn_name_entry_in_session_state]

        rxn_ratings.append(rxn_rating_dict)
    rxn_ratings_dict_df = pd.DataFrame(rxn_ratings)
    st.session_state["eval_details"] = pd.concat([ st.session_state["eval_details"], 
                                                rxn_ratings_dict_df ],
                                                ignore_index=True)

    # saving
    st.session_state["eval_mol"].to_csv("eval_mol.csv",index=False)
    st.session_state["eval_details"].to_csv("eval_details.csv",index=False)

    ws = sheet.worksheet(st.secrets.data.name_mol)
    gspdl.worksheet_append( ws, comment_dict_df )
    ws = sheet.worksheet(st.secrets.data.name_details)
    gspdl.worksheet_append( ws, rxn_ratings_dict_df )
    

def submit_update(molid=None,log=True):
    # log answers
    if log:
        log_evaluation()

    # only update data after submission
    update_filters()
    multi_filtered = st.session_state["prev_data"]

    # update index
    if iterate_by_matchidx:
        st.session_state["data_index"]= generate_index_by_matchid(multi_filtered)
    else:
        st.session_state["data_index"] = generate_index(multi_filtered,molid)

def clear_input():
    if "settings_initialized" in st.session_state:
        submit_update(log=False)
        try:
            trial_molid = int(st.session_state.molid_input)
            update_filters()
            generate_index(st.session_state.prev_data,trial_molid)
        except:
            st.markdown(f"Input molid `{st.session_state.molid_input}` invaild, ignoring.")
        st.session_state["molid_input"] = "" #callbacks are executed before everything, so the session_state can be set *before* the input field is declared/defined.


if st.session_state["b_update_data"]: #in multipage form, make sure we always update when we come back from the settings page, IF THE SETTINGS WERE CHANGED
    submit_update(log=False)


with st.sidebar:
    molid_input = st.text_input("specify molecule ID (optional)",key="molid_input",on_change=clear_input)

    # Submission
    if "userinfo" not in st.session_state \
        or st.session_state["userinfo"] in ["",None] \
        or "@" not in st.session_state["userinfo"]:
        st.markdown("**Enter a valid e-mail on Settings page to submit evaluations.**")
    else:
        with st.form("evaluation",clear_on_submit = True):
            st.markdown("**Functional group quality for...**")
            st.markdown("(scoring scale: `0` to skip, `1`: bad, `3`: interesting, `5`: good)")
            radio_quality_list = []

            multi_filtered = st.session_state["prev_data"]
            match_specific_data = multi_filtered.loc[ st.session_state["data_index"] ]

            rxn_types = match_specific_data.index.unique("rxn_name")

            st.session_state["rxns_for_this_ftn_group"] = []
            for rxn_name in rxn_types:  
                if rxn_name in rxn_name_alias:
                    rxn_name = rxn_name_alias[rxn_name]
                keyname = "rating_" + rxn_name
                st.session_state["rxns_for_this_ftn_group"].append(keyname)
                #radio_quality_list.append( st.radio(rxn_name + " polymerization",
                #                                    ("skip","bad","interesting","good"),
                #                                    horizontal=True,
                #                                    key = keyname)
                #)
                radio_quality_list.append( st.select_slider("**" + rxn_name + " polymerization**",
                                                            ("skip","1: bad","2","3: interesting","4","5: good"),
                                                            key=keyname))

                #radio_quality = st.radio("**Ftnl group quality**",("skip","bad","interesting","good"),horizontal=True)
            #radio_quality_list.append( st.radio("other polymerization",
            #                                    ("skip","bad","interesting","good"),
            #                                    horizontal=True,
            #                                    key="rating_other") 
            #)
            with st.expander("manually enter other polymerization for this functional group",expanded=False):
                radio_quality_list.append( st.select_slider("**other polymerization**",
                                                        ("skip", "1: bad", "2", 
                                                        "3: interesting",
                                                        "4","5: good"),
                                                        key="rating_other"))
                st.session_state["rxns_for_this_ftn_group"].append("rating_other")

                text_form = st.text_area("comments on the highlighted functional group: (use atom indices if needed)","",key="comments_ftn")

            #radio_quality = st.radio("**Overall monomer quality**",("no comment","bad","interesting","good"),horizontal=True,key="rating_mol")
            radio_quality = st.select_slider("**Overall monomer quality**",
                                                        ("skip","1: bad","2",
                                                        "3: interesting",
                                                        "4","5: good"),
                                                        key="rating_mol")
            text_form = st.text_area("comments on the monomer: (use atom indices if needed)","",key="comments_mol")

            submitted = st.form_submit_button("submit",on_click=submit_update)
    #end if userinfo clause

# ===== Actual main content
# retrieve index, molecule, match
if iterate_by_matchidx:
    molid,rxn_name,ftn_id,match_id = st.session_state["data_index"]
    match_specific_data = multi_filtered.loc[(molid,rxn_name,ftn_id,match_id)]
    if isinstance(match_specific_data,pd.Series):
        match_specific_data = match_specific_data.to_frame().transpose()
        match_specific_data.index.names = ["molid","rxn_name","ftn_id","matchid"]

    evaluation,description = characterize_substituents_by_matchid(match_specific_data)

    # characterize and draw
    compound_smiles = data_rxn.loc[molid].smiles
    mol = Chem.MolFromSmiles(compound_smiles)

    ftn_group_ids = ast.literal_eval(match_specific_data.iloc[0].matchidx)
    #core_bool = mychem.get_max_scaffold(mol,ftn_group_ids)
    #highlightcolors,highlightbonds = mychem.color_scaffold(mol,core_bool)
    highlightcolors,highlightbonds = mychem.color_ftn(mol,ftn_group_ids)
    im = mychem.highlight_draw(mol,highlightcolors,highlightbonds,format="svg")
else:
    molid,matchidx = st.session_state["data_index"]

    tmp_rxn_name = multi_filtered.query("molid == @molid").query("matchidx == @matchidx").index.get_level_values("rxn_name")[0] #should be the same characterization features for every rxn name, so any will do
    match_specific_data = multi_filtered.loc[(molid,matchidx,tmp_rxn_name)]
    if isinstance(match_specific_data,pd.Series):
        match_specific_data = match_specific_data.to_frame().transpose()
        match_specific_data.index.names = ["molid","matchidx","rxn_name"]

    #evaluation,description = characterize_substituents(match_specific_data) HERE
    evaluation,description = characterize_substituents(match_specific_data)

    # characterize and draw
    compound_smiles = data_rxn.loc[molid].smiles
    mol = Chem.MolFromSmiles(compound_smiles)

    ftn_group_ids = ast.literal_eval(match_specific_data.index.get_level_values("matchidx")[0])
    highlightcolors,highlightbonds = mychem.color_ftn(mol,ftn_group_ids)
    im = mychem.highlight_draw(mol,highlightcolors,highlightbonds,format="svg")


# ===== Display
st.markdown(st.session_state["data_index_description"])
st.image(im)
st.markdown(f"{evaluation}")
checkbox_details = st.checkbox('Show substituent analysis below')
if checkbox_details:
    st.markdown(description)


# ===== Misc

