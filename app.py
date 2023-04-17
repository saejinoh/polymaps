# Main imports
import streamlit as st
import time
import os
import ast
homedir = os.path.dirname(__file__)

# Analysis imports
import random
import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem import rdChemReactions
from rdkit.Chem import Draw
import mychem
iterate_by_matchidx = False

@st.cache_data
def load_data():
    # functional group data
    filename = homedir + "/../data/" + "fg_analysis_2023-04-10.csv"
    data_df = pd.read_csv(filename,index_col=False)
    #def toset(mystr):
    #   return frozenset( ast.literal_eval(mystr))
    #data_df["matchidx"] = data_df["matchidx"].map(toset)
    if iterate_by_matchidx:
        #multi = data_df.set_index(["molid","matchidx","rxn_name"])
        #st.write(multi.loc[(0,'(2, 3)','simple')])
        multi   = data_df.set_index(["molid","rxn_name","ftn_id","matchid"])
        #st.write(multi.loc[(0,'simple',0,1)])
    else:
        multi = data_df.set_index(["molid","matchidx","rxn_name"])

    # functional group count data; makes initial filtering easier
    filename = homedir + "/../data/" + "rxntypes_2023-04-10.csv"
    data_rxn = pd.read_csv(filename,index_col=False)

    # evaluations
    # evaluate #ftnl groups identified
    if "num_ftn" not in data_rxn.index:
        data_rxn["num_ftn"] = 0.

        num_ftn = multi.reset_index().set_index("molid").groupby(["molid"]).agg("nunique").matchidx
        data_rxn.loc[multi.index.unique("molid"),"num_ftn"] = num_ftn

    # evaluate MW
    if "MW" not in data_rxn.index:
        data_rxn["MW"] = 0.
        def get_MW(row):
            mol = Chem.MolFromSmiles(row.smiles)
            row.MW = Chem.Descriptors.MolWt(mol)
            return row
        data_rxn = data_rxn.apply(get_MW,axis=1)
    
    return multi,data_rxn

multi,data_rxn = load_data()

if "max_MW" not in st.session_state:
    st.session_state["max_MW"] = data_rxn.MW.max()
if "max_numftn" not in st.session_state:
    st.session_state["max_numftn"] = int(data_rxn.num_ftn.max())


# Key variables
quality_dict = {0:"poor",1:"decent",2:"promising"}
eprop_dict = {-1:"withdrawing", 1:"donating", 0:"undetermined"}
index = {"molid":0, "rxn_name":1, "ftn_id":2, "matchid":3}
rxn_types = data_rxn.columns
# TMP: eliminate step reactions
rxn_types = [x for x in rxn_types if not x.startswith("step") and x != "smiles"]


# Formatting
st.markdown("""
<style>
.big-font {
    font-size:25px !important;
}
</style>
<style>
.med-font {
    font-size:20px !important;
}
</style>
""", unsafe_allow_html=True)


def set_update_data_flag(flag):
    st.session_state["b_update_data"] = flag


# ===== Chemical Utilities
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
    #TMP, remove step reactions
    rxn_names = data.index.get_level_values(1)
    step_rxn_names = [x for x in rxn_names if x.startswith("step")]

    # filter by reaction name. TODO: step reactions needs adjustment
    if rxn_name is None or rxn_name == "choose for me!":
        # TMP, remove step reactions
        sub_data = data.iloc[ ~data.index.get_level_values(1).isin(step_rxn_names) ]
        inds = sub_data.index.get_level_values(0).unique().values
    else:
        # filter by rxn_name
        inds = np.argwhere( (data_rxn[rxn_name]>0).values ).flatten() #alternative way to filter
        #sub_data = data.query("molid in @inds")
        rxns_of_interest = [rxn_name]
        sub_data = data.query("rxn_name in @rxns_of_interest")
        #sub_data = sub_data[ sub_data.index.get_level_values(rxn_name).isin([rxn_name])] 

    # filter by MW
    inds_where_MW_range = data_rxn.loc[ 
        (data_rxn.MW>= slider_MW[0]) & 
        (data_rxn.MW <= slider_MW[1]) ].index.values

    sub_data = sub_data.query( "molid in @inds_where_MW_range")

    # filter by number of functional groups
    inds_where_num_ftn_range = data_rxn.loc[ 
        (data_rxn.num_ftn>= slider_num_ftn[0]) & 
        (data_rxn.num_ftn <= slider_num_ftn[1]) ].index.values
    sub_data = sub_data.query( "molid in @inds_where_num_ftn_range")

    # return
    return sub_data

def generate_index(multi,data_rxn,rxn_selection):
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
            - subid
        And also saves to st.session_state, both the above, and:
            - data_index
            - data_index_description
        (may refactor to clean up in the future)

        Global variables used:
            user inputs
    """
    
    if "initialized" not in st.session_state:
        multi_filtered = filter_rxn(multi,data_rxn,rxn_selection)
        molids = multi_filtered.index.get_level_values(0).unique()
        molnum = molids.values.size
        mol_subid = np.random.randint(0,molnum) #default to random molecule
        molid = molids[mol_subid]

        mol_specific_data = multi_filtered.query("molid == @molid")
        num_ftnl_groups = mol_specific_data.reset_index().agg("nunique").matchidx
        ftn_subid = 0 #default init
        ftn_group_ids = mol_specific_data.index.get_level_values("matchidx").unique()[ftn_subid]
        
        description_base = ""
        st.session_state["initialized"] = True
    else:
        multi_filtered = filter_rxn(multi,data_rxn,rxn_selection)
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

        multi_filtered = filter_rxn(multi,data_rxn,rxn_selection)

        if b_mol_out_of_set:
            molids = multi_filtered.index.get_level_values(0).unique()
            molnum = molids.values.size

            description_base = "##### previous molecule out of filtered set, returning random molecule  \n"
            mol_subid = np.random.randint(0,molnum) 
            molid = molids[mol_subid]
            ftn_subid = 0 #always default to first ftnl group
        elif b_finished_mol:
            molids = multi_filtered.index.get_level_values(0).unique()
            molnum = molids.values.size

            description_base = "##### new molecule  \n"
            current_mol_subid = np.argwhere( molids == prev_molid )[0][0]
            molid = prev_molid

            reached_last_mol = True if (current_mol_subid >= molnum - 1) else False

            if iteration_selection == "random":
                mol_subid = np.random.randint(0,molnum) #IMPORTANT
                molid = molids[mol_subid]
            elif iteration_selection == "sequential":
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
            molids = multi_filtered.index.get_level_values(0).unique()
            molnum = molids.values.size
            mol_subid = np.argwhere( molids == prev_molid )[0][0] #prev_mol_subid
            molid = molids[mol_subid]
            description_base = "##### continuing with molecule  \n"

            ftn_subid = prev_ftn_subid + 1
            
    mol_specific_data = multi_filtered.query("molid == @molid")
    num_ftnl_groups = mol_specific_data.reset_index().agg("nunique").matchidx
    ftn_group_ids = mol_specific_data.index.get_level_values("matchidx").unique()[ftn_subid]
            
    # Save and return
    numdigits=len(str(molnum))
    numdigits_str = "0"+str(numdigits)

    description = description_base + f"**Molecule ID:**\t\t `{molid:{numdigits_str}}` ({mol_subid+1}/{molnum} monomers identified for the chosen reaction type)  \n**Showing:**\t\t potential functional group `{ftn_subid+1}`/`{num_ftnl_groups}` for rxn type `{rxn_selection}`"
    st.session_state["prev_data"] = multi_filtered
    st.session_state["data_index"] = (molid, ftn_group_ids)
    st.session_state["ftn_tracking"] = (ftn_subid,num_ftnl_groups)
    #st.session_state["mol_subid"] = mol_subid
    st.session_state["data_index_description"] = description

    return (molid,ftn_group_ids),multi_filtered
    
def generate_index_by_matchid(multi,data_rxn,rxn_selection):
    multi = filter_rxn(multi,data_rxn,rxn_selection)
    st.session_state["prev_data"] = multi

    # utility
    molids = multi.index.get_level_values(0).unique()
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

    # figure out generation mode
    #if radio2 == "random" and radio1 == "new molecule":
    #    mode_new_mol = "random"
    #   mode_new_rxn = "random"
    #    mode_new_ftnid = "random"
    #    mode_new_matchid = "true"
    #elif radio2 == "random" and radio1.startswith("same molecule"):
        # only do new (random) molecule once matchids are exhausted
    #    mode_new_mol = "same"
    #    mode_new_rxn = "random"
    #    mode_new_matchid = "random"
    #    mode_new_ftnid = "random"
    #elif radio2 == "sequential" and radio1 == "new molecule":
    #    pass

    # actually select index
    #if radio1 == "new molecule" and radio2 == "random" or first_time or prev_molid not in molids:
    #    subid = np.random.randint(0,molnum) #IMPORTANT
    #elif radio1 == "new molecule":
    #    pass
    #else:
    #    prev_subid = np.argwhere(molids == prev_molid)[0][0]
    #    prev_molid = molids[prev_subid]
    #    subid = np.minimum(prev_subid+1,molnum-1)

    subid = np.random.randint(0,molnum) #IMPORTANT
    molid = molids[subid]
    #molid = 1818 # TICKET: rop-amine, bad rule

    mol_specific_data = multi.query(f"molid == {molid}")

    # if "choose for me!", choose a random reaction
    if rxn_selection == "choose for me!":
        rxn_types = mol_specific_data.index.get_level_values(1).unique().values
        #rxn_types = [ x for x in rxn_types if not x.startswith("step")] #TMP, remove step reactions for now
        #rxn_types
        rxn_name = random.choice(rxn_types)
    else:
        rxn_name = mol_specific_data.index[0][1] # in my scheme, should all be the same reaction name

    if not rxn_name.startswith("step"): #TMP filter out step
        ftn_id = 0
        match_totals = len(mol_specific_data.index.get_level_values(3).unique().values)
        match_id = np.random.randint(0,match_totals) #chooses a random match

    # store and return
    #st.session_state["subid"] = subid
    st.session_state["data_index_description"] = f"**Molecule:**\t\t `{molid:{numdigits_str}}` ({subid}/{molnum} monomers identified for the chosen reaction type)  \n**Rxn type:**\t\t `{rxn_name}`  \n**Showing:**\t\t `{match_id+1}`/`{match_totals}` reactive sites identified"

    return (molid,rxn_name,ftn_id,match_id),multi

def characterize_substituents(match_specific_data):
    ftn_group_ids = ast.literal_eval(match_specific_data.iloc[0].matchidx)

    n_substituents = len(substituents)
    evaluation = f"Functional group (colored in `green`) w/ atoms `{ftn_group_ids}`\n"
    evaluation +=f"    with `{n_substituents}` substituents.\n\n"
    return evaluation, ""

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
    #evaluation += f"\nEstimated quality: {quality_dict[quality]}"
    return evaluation,description


# ===== Navigation Sidebar
with st.sidebar:
    #st.write("Navigation")
    #st.markdown('<u><p class="big-font">Navigation</p><u>',unsafe_allow_html=True)
    #st.markdown('<p class="med-font">On submission present...</p>',unsafe_allow_html=True)
    st.markdown("# Navigation")
    st.markdown("### On next molecule, show...")
    
    rxn_selection = st.selectbox("reaction type",("choose for me!",*rxn_types),
                                 on_change=lambda: set_update_data_flag(True))

    #radio1 = st.radio("iteration mode",("same molecule, new prediction","new molecule"),horizontal=True)

    #radio2 = st.radio("randomization mode:",("sequential","random"),horizontal=True)

    iteration_selection = st.selectbox("molecule iteration mode:",
        ("random","sequential"),index=1, on_change = lambda x:set_update_data_flag(True))

    # MW
    slider_MW = st.slider("MW range",0.,st.session_state.max_MW,(10.,st.session_state.max_MW),
                          on_change = lambda x: set_update_data_flag(True))
    
    # Simplicity/Complexity
    # (# of polymerizations identified, # functional groups, # subsitutents)
    slider_num_ftn = st.slider("number functional groups",1,st.session_state.max_numftn,(1,st.session_state.max_numftn),
                               on_change = lambda x: set_update_data_flag(True))

    # Bulkiness

    
    # Submission
    with st.form("evaluation",clear_on_submit = True):
        st.markdown("##### Quality of this identified polymerization site (in green)?")
        radio_quality = st.radio("rxn quality",("skip","bad","interesting","good"),label_visibility="hidden",horizontal=True)
        text_form = st.text_area("Additional comments (use atom indices to refer to specific atoms):","")

        submitted = st.form_submit_button("submit")

        if submitted:
            # only update data after submission
            if iterate_by_matchidx:
                #b = filter_rxn(multi,data_rxn,rxn_selection)
                #st.session_state["prev_data"] = b
                st.session_state["data_index"],st.session_state["prev_data"] = generate_index_by_matchid(multi,data_rxn,rxn_selection)
            else:
                st.session_state["data_index"],st.session_state["prev_data"] = generate_index(multi,data_rxn,rxn_selection)
            set_update_data_flag(False)

    # Other comments
    st.markdown("---")
    with st.form("general comments"):
        comment_area = st.text_area("General comments?","")
        submitted = st.form_submit_button("submit")


# ===== Test
#rxn_name = "rop-amine"
#a,b = filter_rxn(multi,data_rxn,rxn_selection)
# initialize data
if "prev_data" not in st.session_state: #first time through
    if iterate_by_matchidx:
        st.session_state["data_index"],multi_filtered = generate_index_by_matchid(multi,data_rxn,None)
        st.session_state["prev_data"] = multi_filtered
    else:
        st.session_state["data_index"],multi_filtered = generate_index(multi,data_rxn,None)
        st.session_state["prev_data"] = multi_filtered
else:
    multi_filtered = st.session_state["prev_data"]


molids = multi_filtered.index.get_level_values(0).unique()
molnum = molids.values.size


# ===== Actual main content
# Other
#compound_smiles = 'c1cc(C(=O)O)c(OC(=O)C)cc1'
#compound_smiles = 'CC[C@H](C)C=CC1=CC2=C(C(=O)[C@@]3(C(=C(C(=O)O3)C(=O)C(=CC)C)C2=CN1CCCC(=O)O)C)Cl'

# retrieve index, molecule, match
if iterate_by_matchidx:
    molid,rxn_name,ftn_id,match_id = st.session_state["data_index"]
    match_specific_data = multi_filtered.loc[(molid,rxn_name,ftn_id,match_id)]
    evaluation,description = characterize_substituents_by_matchid(match_specific_data)
    if isinstance(match_specific_data,pd.Series):
        match_specific_data = match_specific_data.to_frame().transpose()
        match_specific_data.index.names = ["molid","rxn_name","ftn_id","matchid"]

    # characterize and draw
    compound_smiles = data_rxn.loc[molid].smiles
    mol = Chem.MolFromSmiles(compound_smiles)

    ftn_group_ids = ast.literal_eval(match_specific_data.iloc[0].matchidx)
    #core_bool = mychem.get_max_scaffold(mol,ftn_group_ids)
    #highlightcolors,highlightbonds = mychem.color_scaffold(mol,core_bool)
    highlightcolors,highlightbonds = mychem.color_ftn(mol,ftn_group_ids)
    im = mychem.highlight_draw(mol,highlightcolors,highlightbonds)
else:
    molid,matchidx = st.session_state["data_index"]

    tmp_rxn_name = multi_filtered.query("molid == @molid").query("matchidx == @matchidx").index.get_level_values("rxn_name")[0]

    match_specific_data = multi_filtered.loc[(molid,matchidx,tmp_rxn_name)]
    if isinstance(match_specific_data,pd.Series):
        match_specific_data = match_specific_data.to_frame().transpose()
        match_specific_data.index.names = ["molid","matchidx","rxn_name"]

    #evaluation,description = characterize_substituents(match_specific_data) HERE
    evaluation = "Test"
    description = "dest_description"

    compound_smiles = data_rxn.loc[molid].smiles
    mol = Chem.MolFromSmiles(compound_smiles)

    # characterize and draw
    compound_smiles = data_rxn.loc[molid].smiles
    mol = Chem.MolFromSmiles(compound_smiles)

    ftn_group_ids = ast.literal_eval(match_specific_data.index.get_level_values("matchidx")[0])
    highlightcolors,highlightbonds = mychem.color_ftn(mol,ftn_group_ids)
    im = mychem.highlight_draw(mol,highlightcolors,highlightbonds)



# ===== Display
st.image(im)
st.markdown(st.session_state["data_index_description"])
st.markdown(f"{evaluation}")
checkbox_details = st.checkbox('Show substituent analysis below')
if checkbox_details:
    st.markdown(description)


# ===== Misc

#st.markdown("# Navigation")
#st.markdown("### On submission, show...")
#selection = st.selectbox("reaction type2",("choose for me!","next"))
#radio1 = st.radio("",("same molecule2, new prediction","new molecule"),horizontal=True)
#radio2 = st.radio("",("sequential2","random"),horizontal=True)
#st.markdown("### Polymerization monomer browser")


#col1,col2 = st.columns(2)
#with col1:
#    st.markdown(f"**Molecule:**\t\t `{molid:{numdigits_str}}` ({subid}/{molnum} monomers identified for the chosen reaction type)  \n**Rxn type:**\t\t `{rxn_name}`  \n**Showing:**\t\t `{match_id+1}`/`{match_totals}` reactive sites identified")

#with col2:
    #rxn_id   = 1
    #rxn_totals = 2
    #({rxn_id+1}/{rxn_totals} reaction types for this molecule)
    #st.write("Molecule: \t",molid,"/",molnum)
    #st.write("Rxn type: \t",rxn_name)

    #st.text(f"{evaluation}")
#    st.markdown(f"{evaluation}")
    #previous = st.form_submit_button("previous")

