import streamlit as st
from rdkit import Chem
from rdkit.Chem import rdChemReactions
from rdkit.Chem import Draw

st.markdown("# Results")

if "prev_data" not in st.session_state or not "settings_initialized" in st.session_state:
    st.markdown("No molecules evaluated yet.")
else:
    data = st.session_state["eval_mol"]

    molids = data.molid.unique()
    mols = []
    smiles = []
    legends = []
    for molid in molids:
        tmp_data = data.loc[ data.molid == molid ]
        ratings = tmp_data.rating_mol.values
        if "interesting" in ratings or "good" in ratings:
            smiles.append( tmp_data.iloc[0].smiles )
            mols.append( Chem.MolFromSmiles(smiles[-1]) )
            legends.append( f"{str(molid)}" )
    svg = Chem.Draw.MolsToGridImage(mols,molsPerRow=3,subImgSize=(250,250),legends=legends, useSVG=True)

    st.markdown(f"### `{len(mols)}`/`{molids.size}` molecules were rated as `interesting` or `good`:")
    st.image(svg)

# Diagnostic stuff
st.markdown("-----")
st.markdown("optional diagnostics (for debugging)")
showlog = st.checkbox("see full user log")
if showlog and "settings_initialized" in st.session_state:
    st.session_state["eval_mol"]
    st.session_state["eval_details"]
    st.session_state["eval_general"]

#show_session_state = st.checkbox("see full session state")
#if show_session_state:
#    st.session_state