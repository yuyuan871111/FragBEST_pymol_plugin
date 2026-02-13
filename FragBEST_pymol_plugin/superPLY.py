"""
This file is a newly-added file (compared to masif_pymol_plugin).

FragBEST_pymol_plugin Copyright (c) 2024-2026 Yu-Yuan (Stuart) Yang /
Arianna Fornili's Lab

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at
    http://www.apache.org/licenses/LICENSE-2.0
"""

from pymol import cmd


def super_ply(
    pdb_query: str,
    ply_query: str,
    pdb_ref: str,
):
    """
    pdb_query: str (object name of the query pdb)
    ply_query: str (object name, the surface ply to describe pdb_query)
    pdb_ref: str (object name of the reference pdb to superimpose)
    """

    cmd.super(pdb_query, pdb_ref)
    transformation_matrix = cmd.get_object_matrix(pdb_query)
    cmd.set_object_ttt(ply_query, transformation_matrix)
