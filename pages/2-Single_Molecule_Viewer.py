# Main imports
import time, os, ast, copy

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

# Streamlit and web imports
import streamlit as st
import st_utils, app_utils
st_utils.set_font()
st_utils.set_sidebar(max_width=650, min_width=550)

# States 
st.session_state.state["reload_batch_evaluation"] = True

# Preamble
if "settings_initialized" not in st.session_state:
    st.markdown("# Upon browser refresh, please revisit Settings page first.")
    st.stop()
rating_on = False

# Database connection
import gspread_pdlite as gspdl
sheet = st.session_state["backend_sheet"]

# ===== Settings and initialization
quality_dict = {0:"poor",1:"decent",2:"promising"} # unused
index = {"molid":0, "rxn_name":1, "ftn_id":2, "matchid":3} # unused

# State tracking
if "b_update_data_single" not in st.session_state:
    st.session_state["b_update_data_single"] = False
def set_update_data_flag(flag):
    st.session_state["b_update_data_single"] = flag

# ===== Loading
# Load data
if "base_data" not in st.session_state:
    multi,data_rxn = app_utils.load_data()
    st.session_state["base_data"] = (multi,data_rxn)
else:
    multi,data_rxn = st.session_state["base_data"]

# ===== Chemical Utilities
def random_molecule(data,rxn_name = None):
    """return a random molecule from the set"""
    molids = data.index.levels[0].unique().values
    num_mols = len(molids)

    # choose random mol
    molid = np.random.randint(0,num_mols)

    subdata = multi.loc[(molid,slice(None),slice(None),slice(None))]

    return molid, num_mols, subdata

def generate_index(multi_filtered,trial_molid=None):
    """Generate next molecule's index (and data) if needed

    Args:
        multi (_type_): _description_
        data_rxn (_type_): _description_
        rxn_selection (_type_): _description_

    Returns:
        tuple: (index tuple)

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
            molid = trial_molid # will throw an error, which will get caught, then return to original behavior. probably not idea that it's so segmented.
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
        substituents.append( (tmp_data.subid, mychem.eprop_dict[tmp_data.electronics], tmp_data.bulkiness))

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
        substituents.append( (tmp_data.subid, mychem.eprop_dict[tmp_data.electronics], tmp_data.bulkiness))

    ftn_group_ids = ast.literal_eval(match_specific_data.iloc[0].matchidx)

    n_substituents = len(substituents)
    evaluation = f"ID'd functional group (colored in `green`) w/ atoms `{ftn_group_ids}`\n"
    evaluation +=f"    with `{n_substituents}` substituents.\n\n"
    description = ""
    for entry in substituents:
        description += f"- root atom `{entry[0]}`\n  - electronic properties: `{entry[1]}`\n  - bulkiness: `{entry[2]:0.2}`\n"
    return evaluation,description


# ===== Load Data
if "prev_data" not in st.session_state \
    or "data_index" not in st.session_state \
    or "data_index_description" not in st.session_state: #first time through
    if "rxn_selection" in st.session_state:
        multi_filtered = app_utils.filter_rxn(multi,data_rxn,st.session_state.rxn_selection)
    else:
        multi_filtered = app_utils.filter_rxn(multi,data_rxn,None)
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

multi_filtered0 = app_utils.filter_rxn(multi,data_rxn,None) #default data set if filters fail


# ===== Sidebar and submission logic
def update_filters():
    if st.session_state["b_update_data_single"]:
        tmp_multi_filtered = app_utils.filter_rxn(multi,data_rxn,st.session_state.rxn_selection)
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
        st.session_state["b_update_data_single"] = False
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
    if rating_on:
        comment_dict["rating_mol"] = st.session_state["rating_mol"]
    else:
        comment_dict["rating_mol"] = ""
    if "skip" not in comment_dict["rating_mol"] \
        and comment_dict["comments_ftn"].strip() != "" \
        and comment_dict["comments_ftn"].strip() != "":
        comment_dict_df = pd.DataFrame([comment_dict])
        st.session_state["eval_mol"] = pd.concat([ st.session_state["eval_mol"], 
                                                    comment_dict_df],
                                                    ignore_index=True)

    # rxn-specific data
    if rating_on:
        rxn_ratings = []
        for rxn_name_entry_in_session_state in st.session_state["rxns_for_this_ftn_group"]:
            rxn_rating_dict = copy.copy(root_dict)
            rxn_name = rxn_name_entry_in_session_state.removeprefix("rating_")
            if rxn_name == "other":
                rxn_name = st.session_state["rating_other_name"]
            rxn_rating_dict["rxn_name"] = rxn_name
            rxn_rating_dict["rating"] = st.session_state[rxn_name_entry_in_session_state]

            if "skip" not in rxn_rating_dict["rating"]:
                rxn_ratings.append(rxn_rating_dict)
        rxn_ratings_dict_df = pd.DataFrame(rxn_ratings)
        st.session_state["eval_details"] = pd.concat([ st.session_state["eval_details"], 
                                                    rxn_ratings_dict_df ],
                                                    ignore_index=True)

    # saving
    st.session_state["eval_mol"].to_csv("eval_mol.csv",index=False)
    st.session_state["eval_details"].to_csv("eval_details.csv",index=False)

    if "skip" not in comment_dict["rating_mol"] \
        and comment_dict["comments_ftn"].strip() != "" \
        and comment_dict["comments_ftn"].strip() != "":
        ws = sheet.worksheet(st.secrets.data.name_mol)
        gspdl.worksheet_append( ws, comment_dict_df )

    if rating_on and rxn_ratings_dict_df.size > 0:
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
    """Action upon entering a specific input
    """
    if "settings_initialized" in st.session_state:
        submit_update(log=False)
        try:
            trial_molid = int(st.session_state.molid_input)
            update_filters()
            generate_index(st.session_state.prev_data,trial_molid)
        except:
            pass
            #st.markdown(f"Input molid `{st.session_state.molid_input}` invaild, ignoring.")
        st.session_state["molid_input"] = "" #callbacks are executed before everything, so the session_state can be set *before* the input field is declared/defined.


if st.session_state["b_update_data_single"]: #in multipage form, make sure we always update when we come back from the settings page, IF THE SETTINGS WERE CHANGED
    submit_update(log=False)
    # note that sidebar loads data again, multi_filtered = st.session_state["prev_data"]


with st.sidebar:
    molid_input = st.text_input("specify molecule ID (optional)",key="molid_input",on_change=clear_input)

    # Submission
    if "userinfo" not in st.session_state \
        or st.session_state["userinfo"] in ["",None] \
        or "@" not in st.session_state["userinfo"]:
        st.markdown("**Enter a valid e-mail on Settings page to submit evaluations.**")
    else:
        with st.form("evaluation",clear_on_submit = True):
            st.markdown("**Functional group quality for the following polymerizations...**")
            #st.markdown("(scoring scale: `0` to skip, `1`: bad, `3`: interesting, `5`: good)")
            if rating_on:
                st.markdown("(scoring scale: `1`: unworkable to `5`: almost definitely works)")
            radio_quality_list = []

            multi_filtered = st.session_state["prev_data"] 
            # data might have been updated by batch_evaluation, 
            # in which case need to ensure data_index is still in the new data set
            try:
                match_specific_data = multi_filtered.loc[ st.session_state["data_index"] ]
            except:
                st.session_state["data_index"] = multi_filtered.index[0][:2]
                match_specific_data = multi_filtered.loc[ st.session_state["data_index"] ]
            

            rxn_types = match_specific_data.index.unique("rxn_name")

            st.session_state["rxns_for_this_ftn_group"] = []
            if rating_on:
                for rxn_name in rxn_types:  
                    if rxn_name in app_utils.rxn_name_alias:
                        rxn_name = app_utils.rxn_name_alias[rxn_name]
                    keyname = "rating_" + rxn_name
                    st.session_state["rxns_for_this_ftn_group"].append(keyname)
                    #radio_quality_list.append( st.radio(rxn_name + " polymerization",
                    #                                    ("skip","bad","interesting","good"),
                    #                                    horizontal=True,
                    #                                    key = keyname)
                    #)
                    radio_quality_list.append( st.selectbox("**" + rxn_name + " polymerization**",
                                                                app_utils.rating_scale,
                                                                key=keyname))

                    #radio_quality = st.radio("**Ftnl group quality**",("skip","bad","interesting","good"),horizontal=True)
                #radio_quality_list.append( st.radio("other polymerization",
                #                                    ("skip","bad","interesting","good"),
                #                                    horizontal=True,
                #                                    key="rating_other") 
                #)
                with st.expander("manually enter other polymerization motif for this functional group",expanded=False):
                    rxn_types_with_other = ["other (write in comments)"]
                    rxn_types_with_other.extend(st.session_state.rxn_types)
                    other_polymerization_motif = st.selectbox("**choose polymerization motif**",rxn_types_with_other,key="rating_other_name")
                    radio_quality_list.append( st.selectbox("**rating**",
                                                            app_utils.rating_scale,
                                                            key="rating_other"))
                    st.session_state["rxns_for_this_ftn_group"].append("rating_other")

                #radio_quality = st.radio("**Overall monomer quality**",("no comment","bad","interesting","good"),horizontal=True,key="rating_mol")
                radio_quality = st.selectbox("**Is the molecule interesting overall?**",
                                                            app_utils.rating_scale,
                                                            key="rating_mol")
            

            tmp_df = st.session_state.eval_mol
            current_molid, current_molidx = st.session_state["data_index"]
            tmp_slice = tmp_df[ (tmp_df.molid == current_molid) 
                                & (tmp_df.molidx == current_molidx) ]
            if tmp_slice.size > 0:
                comment_0 = tmp_slice.iloc[-1].comments_mol
            else:
                comment_0 = ""

            text_form = st.text_area("comments on the highlighted functional group: (use atom indices if needed)",
                                     "",key="comments_ftn")
            text_form = st.text_area("comments on the monomer, overall: (use atom indices if needed)",
                                     "",key="comments_mol")

            submitted = st.form_submit_button("submit and next",on_click=submit_update)
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

