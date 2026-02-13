"""
This file is modified from masif_pymol_plugin.
masif_pymol_plugin: https://github.com/LPDI-EPFL/masif/tree/master/source/masif_pymol_plugin

Original work Copyright (c) 2019 Gainza P, Sverrisson F, Monti F, Rodola,
Bronstein MM, Correia BE
FragBEST_pymol_plugin Copyright (c) 2024-2026 Yu-Yuan (Stuart) Yang /
Arianna Fornili's Lab

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at
    http://www.apache.org/licenses/LICENSE-2.0

Original header:
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
