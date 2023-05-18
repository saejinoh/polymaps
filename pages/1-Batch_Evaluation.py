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

# Streamlit imports
import streamlit as st
st.set_page_config(layout="wide")
import st_utils, app_utils
st_utils.set_font()

# Preamble
if "settings_initialized" not in st.session_state:
    st.markdown("# Upon browser refresh, please revisit Settings page first.")
    st.stop()

# Database connection
import gspread_pdlite as gspdl
sheet = st.session_state["backend_sheet"] #alias for convenience

# ===== Settings and initialization
# Constants
from app_utils import rating_scale

# General settings
if "entry_other_reset" not in st.session_state:
    st.session_state["entry_other_reset"] = True

# State tracking
if "b_update_data_batch" not in st.session_state:
    st.session_state["b_update_data_batch"] = False
def set_update_data_flag(flag):
    st.session_state["b_update_data_batch"] = flag

if "batch_page" not in st.session_state:
    st.session_state["batch_page"] = 1 #change to 0

if "mols_per_page" not in st.session_state:
    st.session_state["mols_per_page"] = 24

# ===== Loading
# ----- Raw data
app_utils.first_load() #complete loading of unfiltered data if it hasn't been loaded yet.

# ----- Filtered data
# Filter if not completed on the settings page.
# Filtering is not too slow, is ok to redo?

# update filters, data, and proxy index
if st.session_state["b_update_data_batch"]:
    #since updating proxies is pretty cheap, always update proxy if filter is updated
    #don't do special multi_filtered checks
    app_utils.update_filters()
    st.session_state["batch_page"] = 1

# retrieve multi_filtered and prune (keep only one entry per unique multi index)
multi_filtered = st.session_state["prev_data"]
if multi_filtered.size == 0:
     st.markdown("##### No molecules left to analyze. Adjust filters on settings page.")
     st.stop()
multi_filtered = multi_filtered[ ~multi_filtered.index.duplicated(keep='first') ]

# update proxy index
if st.session_state["b_update_data_batch"] or "proxy_indices" not in st.session_state:
    proxy_indices = np.arange(multi_filtered.shape[0])
    if st.session_state.iteration_selection == "random":
        np.random.shuffle(proxy_indices)
    st.session_state["proxy_indices"] = proxy_indices
    st.session_state["b_update_data_batch"] = False
else:
    proxy_indices = st.session_state["proxy_indices"]


# ===== Page number
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


# ===== generate images of molecules, all the information
data_slice = multi_filtered.iloc[proxy_indices[ind0:ind1]]
#data_slice
#data_rxn.iloc[data_slice.index.unique("molid")]

mols = []
svgs = []
molids = []
molidxs = []
rxn_names = []
for row in data_slice.iterrows():
    mols.append( Chem.MolFromSmiles(row[1].smiles) )
    molids.append( row[0][0] )
    rxn_name = row[0][2]
    rxn_names.append( app_utils.rxn_name_alias[rxn_name] if rxn_name in app_utils.rxn_name_alias else rxn_name )
    molidxs.append( row[0][1] ) #leave as string
    ftn_group_ids = ast.literal_eval(row[0][1]) #depends on ordering of index; should get ids of atoms in ftnl group.
    

    highlightcolors,highlightbonds = mychem.color_ftn(mols[-1],ftn_group_ids)
    im = mychem.highlight_draw(mols[-1],highlightcolors,highlightbonds,format="svg",
                               wd=250,ht=250)
    svgs.append(im)

# ===== Display
# --- actions
def log():
    """save results of this page
    Notes:
    - uses the following global variables:
        - entries (in session state)
        - entries_general (in session state)
        - eval_mol and eval_details (in session state from 0-Settings)
        - rxn_names
        - molids
    """
    timestamp = pd.Timestamp(time.time(),unit="s")

    mol_ratings = []
    rxn_ratings = []
    for ii in range(data_slice.shape[0]):
        root_dict = {}
        molid, molidx = molids[ii],molidxs[ii]
        rxn_name = rxn_names[ii]

        # base information in dictionary
        root_dict["molid"],root_dict["molidx"] = molid, molidx
        root_dict["smiles"] = st.session_state["data_rxn"].loc[molid].smiles
        root_dict["userinfo"] = st.session_state["userinfo"]
        root_dict["timestamp"] = timestamp

        #general ratings
        # --> comments first
        comment_dict = copy.copy(root_dict)
        if "skip" not in st.session_state[f"entry_other_{ii}"]:
            comment_dict["comments_ftn"] = \
                st.session_state[f"entry_other_name_{ii}"] \
                + "; " + st.session_state[f"entry_other_{ii}"]
        else:
            comment_dict["comments_ftn"] = ""
        comment_dict["comments_mol"] = st.session_state[f"entry_comments_{ii}"].strip()
        # --> ratings
        comment_dict["rating_mol"] = st.session_state[f"entry_general_{ii}"]
        # --> check if we should store or not
        tmp_df = st.session_state.eval_mol
        tmp_slice = tmp_df[ (tmp_df.molid == molid) 
                            & (tmp_df.molidx == molidx) ]
        if "skip" in comment_dict["rating_mol"]\
            and comment_dict["comments_ftn"] == ""\
            and comment_dict["comments_mol"].strip() == "":
            #always skip if blank result
            pass
        else: #check if need to save
            if tmp_slice.size == 0:
                #always save if no previous results (and already know that it's not blank)
                mol_ratings.append(comment_dict)
            else:
                if comment_dict["comments_ftn"] not in tmp_slice.comments_ftn.values \
                    or comment_dict["comments_mol"] not in tmp_slice.comments_mol.values \
                    or comment_dict["rating_mol"] != tmp_slice.iloc[-1].rating_mol:
                    #only save if any of the values are different from previous results
                    mol_ratings.append(comment_dict)
        # --> reset
        st.session_state[f"entry_general_{ii}"] = rating_scale[0] 

        #detailed ratings
        rxn_rating_dict = copy.copy(root_dict)
        rxn_rating_dict["rxn_name"] = rxn_name
        rxn_rating_dict["rating"] = st.session_state[f"entry_{ii}"]

        if "skip" not in rxn_rating_dict["rating"]:
            tmp_df = st.session_state.eval_details
            tmp_slice = tmp_df[ (tmp_df.molid == molid)
                            & (tmp_df.molidx == molidx)
                            & (tmp_df.rxn_name == rxn_name)   ]
            if tmp_slice.size > 0:
                if rxn_rating_dict["rating"] != tmp_slice.iloc[-1].rating:
                    #only save if different from previous result
                    rxn_ratings.append(rxn_rating_dict)
            else: #not in previous results
                rxn_ratings.append(rxn_rating_dict)
        st.session_state[f"entry_{ii}"] = rating_scale[0] #reset

        #other ratings
        rxn_rating_dict = copy.copy(root_dict)
        rxn_rating_dict["rxn_name"] = st.session_state[f"entry_other_name_{ii}"]
        rxn_rating_dict["rating"] = st.session_state[f"entry_other_{ii}"]
        if "skip" not in rxn_rating_dict["rating"]:
            tmp_df = st.session_state.eval_details
            tmp_slice = tmp_df[ (tmp_df.molid == molid)
                            & (tmp_df.molidx == molidx)
                            & (tmp_df.rxn_name == rxn_name)   ]
            if tmp_slice.size > 0:
                if rxn_rating_dict["rating"] != tmp_slice.iloc[-1].rating:
                    #only save if different from previous result
                    rxn_ratings.append(rxn_rating_dict)
            else: #not in previous results
                rxn_ratings.append(rxn_rating_dict)

    # save
    mol_ratings_df = pd.DataFrame(mol_ratings)
    st.session_state["eval_mol"] = pd.concat([ st.session_state["eval_mol"], 
                                                 mol_ratings_df],
                                                ignore_index=True)
    rxn_ratings_df = pd.DataFrame(rxn_ratings)
    st.session_state["eval_details"] = pd.concat([ st.session_state["eval_details"], 
                                                rxn_ratings_df ],
                                                ignore_index=True)
    
    st.session_state["eval_mol"].to_csv("eval_mol.csv",index=False)
    st.session_state["eval_details"].to_csv("eval_details.csv",index=False)

    if mol_ratings_df.size > 0:
        ws = sheet.worksheet(st.secrets.data.name_mol)
        gspdl.worksheet_append( ws, mol_ratings_df )
    if rxn_ratings_df.size > 0:
        ws = sheet.worksheet(st.secrets.data.name_details)
        gspdl.worksheet_append( ws, rxn_ratings_df )
# end log

valid_email = not ("userinfo" not in st.session_state \
                or st.session_state["userinfo"] in ["",None] \
                or "@" not in st.session_state["userinfo"])

def get_closest_page():
    st.session_state.mols_per_page = st.session_state.mols_per_page_widget
    st.session_state["batch_page"] = 1 + int( np.floor( ind0/st.session_state.mols_per_page ) )
    if valid_email:
        log()
    st.session_state["entry_other_reset"] = True
def page_forward():
    st.session_state.batch_page = min(last_page_index+1, st.session_state.batch_page + 1)
    if valid_email:
        log()
    st.session_state["entry_other_reset"] = True
def page_backward():
    st.session_state.batch_page = max(1, st.session_state.batch_page - 1)
    if valid_email:
        log()
    st.session_state["entry_other_reset"] = True
def store_page():
    st.session_state.batch_page = int(st.session_state.batch_page_text)
    st.session_state.batch_page_text = ""
    if valid_email:
        log()
    st.session_state["entry_other_reset"] = True
def save_only():
    if valid_email:
        log()
    
# --- Sidebar
with st.sidebar:
    if not valid_email:
        st.markdown("**Enter a valid e-mail on Settings page to submit evaluations.**")

    st.markdown(f"## Page `{st.session_state.batch_page}`/`{last_page_index+1}`")

    st.select_slider("**Results per page**",(12,24,36,48),value = st.session_state.mols_per_page,
                     key="mols_per_page_widget",on_change=get_closest_page)
    
    st.text_input("**Save & Jump to page**",value = "", key="batch_page_text",
                  on_change=store_page)
    st.button("Save & Next page",on_click = page_forward)
    st.button("Save & Previous page",on_click = page_backward)
    st.button("Save",on_click = save_only)
    with st.expander("""Rating Guide"""):
        st.markdown(""" 
- skip (don't answer)
- incorrectly ID'd
- `1`: impossible
- `3`: potentially works
- `5`: definitely works
""")
    st_utils.change_widget_fontsize("""Rating Guide""","18px")

# --- Main Text
# Entry Area
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
                    if valid_email:
                        current_molid = molids[index_abs]
                        current_rxn = rxn_names[index_abs]
                        current_molidx = molidxs[index_abs]

                        # DUE to streamlit behavior,
                        # We need to manually figure out when to and not to load and persist 
                        # widget values. Here, we only refresh
                        if st.session_state[f"entry_other_reset"]:
                            st.session_state[f"entry_other_name_{index_abs}"] = "other (write in comments)"
                            st.session_state[f"entry_other_{index_abs}"] = rating_scale[0]
                            st.session_state[f"entry_comments_{index_abs}"] = ""

                        tmp_df = st.session_state.eval_details
                        tmp_slice = tmp_df[ (tmp_df.molid == current_molid)
                                            & (tmp_df.molidx == current_molidx)
                                            & (tmp_df.rxn_name == current_rxn)   ]
                        if f"entry_{index_abs}_load" not in st.session_state:
                            st.session_state[f"entry_{index_abs}_load"] = True
                        if st.session_state[f"entry_{index_abs}_load"]: 
                            if tmp_slice.size > 0:
                                st.session_state[f"entry_{index_abs}"] = tmp_slice.iloc[-1].rating
                            else:
                                st.session_state[f"entry_{index_abs}"] = rating_scale[0]
                            st.session_state[f"entry_{index_abs}_load"] = False #False after first reload
                        #otherwise, just keep result

                        tmp_df = st.session_state.eval_mol
                        tmp_slice = tmp_df[ (tmp_df.molid == current_molid) 
                                            & (tmp_df.molidx == current_molidx) ]
                        if f"entry_general_{index_abs}_load" not in st.session_state:
                            st.session_state[f"entry_general_{index_abs}_load"] = True
                        if st.session_state[f"entry_general_{index_abs}_load"]:
                            if tmp_slice.size > 0:
                                #if previous result exists, and in reload mode
                                st.session_state[f"entry_general_{index_abs}"] = tmp_slice.iloc[-1].rating_mol
                            else:
                                #new molecule, go to default value
                                st.session_state[f"entry_general_{index_abs}"] = rating_scale[0]
                            st.session_state[f"entry_general_{index_abs}_load"] = False #False after first reload
                        #otherwise, just keep result


                        # Actually create the entries
                        entries.append( st.radio(f"**functional group (highlighted) quality for polymerization: `{rxn_names[index_abs]}`**",
                                                rating_scale,
                                                key=f"entry_{index_abs}",horizontal=True
                                                ))
                        st.multiselect("**comments on ftn group** (all that apply)",["radical","anionic","cationic","metathesis",
                                        "too bulky","too electron poor","too electron rich",
                                        "not nucleophilic enough","not electrophilic enough",
                                        "aromatic too stable","other"],key=f"select_comments_{index_abs}")
                        entries_general.append( st.radio("**overall monomer quality for polymerization**",
                                                rating_scale,
                                                key=f"entry_general_{index_abs}",horizontal=True
                                                ))
                        st.multiselect("**comments on monomer overall** (all that apply)",[
                                        "too bulky","too many interfering functional groups","repeat unit will look very different",
                                        "aromatic too stable","other"],key=f"select_comments_general_{index_abs}")

                        #if f"container_state_{index_abs}" not in st.session_state:
                        #    st.session_state[f"container_state_{index_abs}"] = False
                        #if st.session_state[f"entry_other_reset"]:
                        #    st.session_state[f"container_state_{index_abs}"] = False
                        #    container =  st.expander("**Manually enter comments and other reactions**",
                        #        expanded=True) #force open, experimental re run will then close
                        #else:
                        #    container =  st.expander("**Manually enter comments and other reactions**",
                        #        expanded = False)
                        if st.session_state["entry_other_reset"]:
                            container =  st.expander("**Manually enter comments and other reactions**",
                                expanded=True) #force open on reset. Experimental rerun and close on next refresh.
                        else:
                            container =  st.expander("**Manually enter comments and other reactions**",
                                expanded=False)

                        with container:
                            st.text_area("other comments: (use atom indices if needed)",
                                         "",key=f"entry_comments_{index_abs}")
                            rxn_types_with_other = ["other (write in comments)"]
                            rxn_types_with_other.extend(st.session_state.rxn_types)
                            st.multiselect("**suggest another polymerization motif**",
                                         rxn_types_with_other,key=f"entry_other_name_{index_abs}")
                            st.radio(f"**quality for selected (above) polymerization(s)**",
                                                rating_scale,
                                                key=f"entry_other_{index_abs}",horizontal=True
                                                )
                            
                    else:
                        st.write(f"identified for: `{rxn_names[index_abs]}`")
                #data_slice.iloc[ia*3 + ic]

# final state handling
if st.session_state[f"entry_other_reset"]:
    st.session_state[f"entry_other_reset"] = False
    js = '''
    f"""
        <p>{st.session_state.batch_page}</p>
        <script>
            window.parent.document.querySelector('section.main').scrollTo(0, 0);
        </script>
    """,
    height=0
    '''
    st.components.v1.html(js)
    time.sleep(st.session_state.mols_per_page/12 * 1) 
    #roughly need 0.2s per 12 results to open/close all expanders
    #or 0.4s per 12 results if want to properly jump to top of page
    st.experimental_rerun()
