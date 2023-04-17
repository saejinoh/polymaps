# Imports
import datamol as dm
from rdkit import Chem
from rdkit.Chem import rdChemReactions
from rdkit.Chem import Draw
from rdkit.Chem.Draw import rdMolDraw2D

import networkx as nx
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# "Tol" colormap from https://davidmathlogic.com/colorblind
colors_int = [(51,34,136),(17,119,51),(68,170,153),(136,204,238),(221,204,119),(204,102,119),(170,68,153),(136,34,85)]
colors = []
for c in colors_int:
    colors.append( (c[0]/256.,c[1]/256.,c[2]/256.) )


# ===== Processing tools
def mol_to_nx(mol):
    G = nx.Graph()

    for atom in mol.GetAtoms():
        G.add_node(atom.GetIdx(),
                   atomic_num=atom.GetAtomicNum(),
                   formal_charge=atom.GetFormalCharge(),
                   chiral_tag=atom.GetChiralTag(),
                   hybridization=atom.GetHybridization(),
                   num_explicit_hs=atom.GetNumExplicitHs(),
                   is_aromatic=atom.GetIsAromatic())
    for bond in mol.GetBonds():
        G.add_edge(bond.GetBeginAtomIdx(),
                   bond.GetEndAtomIdx(),
                   bond_type=bond.GetBondType())
    return G


# ===== Selection Helpers and Filters
def get_max_scaffold(mol,selection,start=None,end=None):
    """Get maximum scaffold (atoms on simple paths connecting start and end), and return coloring

    Args:
        mol (rdkit Chem Mol): 
        selection (listlike): base selection guaranteed to be in scaffold
        start (int): start index
        end (int): end index

    Returns:
        listlike: array, 0 if atom is pendant, > 0 if in scaffold
    """
    n_atoms = mol.GetNumAtoms()
    g = mol_to_nx(mol)
    
    core_bool = [0] * n_atoms
    for idx in selection:
        core_bool[idx] = 2

    if start is None:
        start = selection[0]
    if end is None:
        end = selection[-1]
    paths = nx.all_simple_paths(g,start,end)

    for p in paths:
        # first check that no node in path is in selection
        # not necessary?  

        # then add nodes to core nodes
        for n in p:
            if core_bool[n] != 2:
                core_bool[n] = 1

    return core_bool

def get_substituent_roots(mol,selection,core=None):
    """Returns atom indices of atom roots attached to core.

    Args:
        mol (rdkit Molecule): 
        selection (listlike): indices
        core (listlike, optional): actual core for which we care about substituents. Defaults to None (core will be set to selection).

    Returns:
        (listlike)
    """
    if core is None:
        core = selection

    seen = [0] * mol.GetNumAtoms()
    substituents = []

    for idx in selection:
        seen[idx] = 1

    for idx in core:
        atom = mol.GetAtomWithIdx(idx)
        for nbr in atom.GetNeighbors():
            if not seen[nbr.GetIdx()]:
                substituents.append(nbr.GetIdx())
                seen[nbr.GetIdx()] = 1

    return substituents

def characterize_substituents(mol,selection,core=None):
    """Returns dictionary with basic characterization

    Args:
        mol (rdkit Molecule): 
        selection (listlike): indices
        core (listlike, optional): actual core for which we care about substituents. Defaults to None (core will be set to selection).

    Returns:
        (dict)
    """
    n_atoms = mol.GetNumAtoms()
    #seen = [0] * n_atoms

    if core is None:
        core = selection

    rgroup_roots = get_substituent_roots(mol,selection,core)

    #for idx in selection:
    #    seen[idx] = 1
    #for idx in rgroup_roots:
    #    seen[idx] = 1

    props = {}
    # === Check different chemical environments
    for root in rgroup_roots:
        props[root] = {}
        atom = mol.GetAtomWithIdx(root)
        # check element
        props[root]["el"] = atom.GetSymbol() #AtomicNum() 

        # check bonding environment
        # maybe not necessary, since the electronics matching does look for = and # bonds
        props[root]["bonds"] = []
        for bond in atom.GetBonds():
            #print(bond)
            idx1 = bond.GetBeginAtomIdx()
            idx2 = bond.GetEndAtomIdx()
            #if seen[idx1]==1 or seen[idx2]==1:
            #    continue
            if idx1 in selection or idx2 in selection:
                continue
            props[root]["bonds"].append(bond.GetBondTypeAsDouble())
        

        # check other functional group environments
        # right now, just do crude EDG/EWG binary classification
        props[root]["electronics"] = check_electron_direction(mol,root,selection)

        # === Check size of pendants, or how to assess bulk of pendants? 
        # esp. for "intramolecular" polymerization
        # I think this is the trickiest...
        # e.g. if it's a long linear pendant piece, then that isn't so bad
        # i.e. MW alone isn't enough?
        # aspherity and eccentricity?
        # or only look up to the 6 closest atoms? go up to 9? (e.g. tert-butyl)
        # and use a measure, e.g. of average distance
        props[root]["bulkiness"],_,_ = check_substituent_bulkiness(mol,root,selection)

    return props

# Evaluate Functional Groups
def eval_functional_group(mol,smarts,base):
    """evaluate quality of functional group matches 

    Args:
        mol (rdchem mol): 
        smarts (rdchem mol, smarts): smarts string or mol initiated from smarts
        base (listlike): list of indices (to the smarts match) for which we care about substituents
    """
    if isinstance(smarts,str):
        pattern = Chem.MolFromSmarts(smarts)
    else:
        pattern = smarts

    matches = mol.GetSubstructMatches(pattern)

    res = []
    for im,m in enumerate(matches):
        core = [ m[ix] for ix in base ] #core = m[base(..)]
        characterization = characterize_substituents(mol,m,core)

        substituents = []
        for k,v in characterization.items():
            v["subid"] = k
            v["matchid"] = im
            v["matchidx"] = m
            substituents.append(v)

        #tmp = {}
        #tmp["matchid"] = im
        #tmp["matchidx"] = m 
        #tmp["match"] = characterization

        res.extend(substituents)
    

    return matches,res



# ===== Visualization Helpers
def color_scaffold(mol,core_bool):
    n_atoms = mol.GetNumAtoms()
    scaffold = []
    free = []
    rgroup_roots = []
    highlightcolors = {}


    color_og = colors[2] #green
    color_rgroup = colors[5] #red
    color_scaffold = colors[3] #blue

    for idx in range(n_atoms):
        if core_bool[idx] == 0:
            free.append(idx)
            highlightcolors[idx] = color_rgroup
        elif core_bool[idx] == 2: #the original selection
            scaffold.append(idx)
            highlightcolors[idx] = color_og
        elif core_bool[idx] == 1:
            highlightcolors[idx] = color_scaffold

    highlightbonds = {}
    for ib,b in enumerate(mol.GetBonds()):
        idx1,idx2 = b.GetBeginAtomIdx(), b.GetEndAtomIdx()
        if core_bool[idx1] == 0 and core_bool[idx2] == 0:
            highlightbonds[ib] = color_rgroup
        elif core_bool[idx1] > 0 and core_bool[idx2] > 0:
            highlightbonds[ib] = color_scaffold
        else: #is bond connecting scaffold with Rgroup
            rgroup_roots.append([ib,(idx1,idx2)])

    #return scaffold, free, rgroup_roots
    return highlightcolors,highlightbonds

def color_ftn(mol,ftn_group_ids):
    n_atoms = mol.GetNumAtoms()
    scaffold = []
    free = []
    rgroup_roots = []
    highlightcolors = {}

    color_og = colors[2] #green
    color_rgroup = colors[5] #blue
    color_scaffold = colors[3] #red

    for idx in range(n_atoms):
        if idx in ftn_group_ids:
            scaffold.append(idx)
            highlightcolors[idx] = color_og
        else:
            highlightcolors[idx] = (1,1,1)

    highlightbonds = {}
    for ib,b in enumerate(mol.GetBonds()):
        idx1,idx2 = b.GetBeginAtomIdx(), b.GetEndAtomIdx()
        if idx1 in ftn_group_ids and idx2 in ftn_group_ids:
            highlightbonds[ib] = color_og

    #return scaffold, free, rgroup_roots
    return highlightcolors,highlightbonds


def highlight_draw(mol,highlight_colors,highlight_bonds={},wd=550,ht=350,format="png",atomIndices=True):
    num = mol.GetNumAtoms()
    
    # prepare drawing
    if format.lower() == "svg":
        d = Chem.Draw.rdMolDraw2D.MolDraw2DSVG(300,300)
    else:
        d = Chem.Draw.rdMolDraw2D.MolDraw2DCairo(wd,ht)
    d.drawOptions().addAtomIndices = atomIndices
    d.drawOptions().annotationFontScale = 0.9

    # prepare highlight options
    highlight = list(range(num))
    highlightRadii = dict( [(ii,0.4) for ii in highlight])

    #draw
    d.DrawMolecule(mol,highlightAtoms=highlight,
                   highlightAtomColors=highlight_colors,
                   highlightBonds=highlight_bonds,
                   highlightBondColors=highlight_bonds,
                   highlightAtomRadii = highlightRadii)
    #d.DrawMoleculeWithHighlights(mol,highlight_atom_map = highlight_colors,
    #                             highlight_bond_map = highlight_bonds,
    #                             highlight_radii = highlightRadii)
    d.FinishDrawing()

    #return
    if format.lower() == "svg":
        svg = d.GetDrawingText().replace('svg:','')
        return svg
    else:
        png = d.GetDrawingText()
        return png
    
