import streamlit as st
st.set_page_config(layout="wide")
from rdkit import Chem
from rdkit.Chem import rdChemReactions
from rdkit.Chem import Draw
import numpy as np
import st_utils, app_utils
st_utils.set_font()

st.markdown("# Results")

if "prev_data" not in st.session_state or not "settings_initialized" in st.session_state:
    st.markdown("No molecules evaluated yet.")
else:
    mols = []
    smiles = []
    good_molids = []
    legends = []

    # process data
    data = st.session_state["eval_mol"]
    molids = data.molid.unique()
    for molid in molids:
        tmp_data = data.loc[ data.molid == molid ]
        ratings = tmp_data.rating_mol.values
        if app_utils.rating_scale[-3] in ratings \
                or app_utils.rating_scale[-2] in ratings \
                or app_utils.rating_scale[-1] in ratings: # rating 3~5
            smiles.append( tmp_data.iloc[0].smiles )
            mols.append( Chem.MolFromSmiles(smiles[-1]) )
            legends.append( f"{str(molid)}" )
            good_molids.append(molid)
    
    data_details = st.session_state["eval_details"]
    molids = data_details.molid.unique()
    for molid in molids:
        if molid not in good_molids: #only check if it's not already included
            tmp_data = data_details.loc[ data_details.molid == molid ]
            ratings = tmp_data.rating.values
            if app_utils.rating_scale[-3] in ratings \
                    or app_utils.rating_scale[-2] in ratings \
                    or app_utils.rating_scale[-1] in ratings: # rating 3~5
                smiles.append( tmp_data.iloc[0].smiles )
                mols.append( Chem.MolFromSmiles(smiles[-1]) )
                legends.append( f"{str(molid)}" )
                good_molids.append(molid)
    
    # get combined, unique, molid list
    all_molids = data.molid.unique().tolist()
    all_molids.extend( data_details.molid.unique() )
    all_molids = np.unique(all_molids).tolist()

    # Draw
    svg = Chem.Draw.MolsToGridImage(mols,molsPerRow=3,subImgSize=(250,250),legends=legends, useSVG=True)

    st.markdown(f"### `{len(mols)}`/`{len(all_molids)}` saved molecules were rated as `> 3`:")
    if len(mols) > 0:
        st.image(svg)

    # Download results
    user_file = f"eval_mol_{st.session_state['userinfo']}.csv"
    #st.session_state["eval_mol"].to_csv(user_file,index=False)
    st.download_button("Download all molecules evaluated",
                       st.session_state["eval_mol"].to_csv(index=False).encode("utf-8"),
                       file_name = user_file
    )


# Diagnostic stuff
st.markdown("-----")
st.markdown("optional diagnostics (for debugging)")
showlog = st.checkbox("see full user log")
if showlog and "settings_initialized" in st.session_state:
    st.session_state["eval_mol"]
    st.session_state["eval_details"]
    st.session_state["eval_general"]


#show_session_state = st.checkbox("see full session state") #for diagnostics
#if show_session_state:
#    st.session_state