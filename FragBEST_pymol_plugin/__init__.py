"""
Modified by Yu-Yuan Yang 2025

Pablo Gainza Cirauqui 2016 LPDI IBI STI EPFL
This pymol plugin for Masif just enables the load ply functions.
"""

from pymol import cmd

from .loadFRAG import load_pdb_with_frags
from .loadPLY import load_ply
from .loadPREDICT import load_grids
from .superPLY import super_ply


def __init_plugin__(app):
    cmd.extend("loadfrag", load_pdb_with_frags)
    cmd.extend("loadpredict", load_grids)
    cmd.extend("loadply", load_ply)
    cmd.extend("superply", super_ply)
