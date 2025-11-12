"""
    loadFRAG.py: This pymol function loads PDB file with fragments into pymol.
    It is modified by Yu-Yuan (Stuart) Yang - Arianna group QMUL 2024-2027
    This is part of probinsense_pymol_plugin.
"""

import numpy as np
from pymol import cmd

from .color_palette import colorDict_for_labels
from .loadPREDICT import get_from


def load_pdb_with_frags(filename, ligand_name):
    # get the object name
    obj_name = filename.replace(".pdb", "")

    # load the object
    cmd.load(filename, obj_name)
    cmd.select("ligand", f"resn {ligand_name}")

    # set the color
    fragment_num = max(get_from("ligand", what="chain"))
    frag_spectrum_list = []
    for chain_id in range(1, int(fragment_num) + 1):
        cmd.set_color(f"frag{chain_id}", colorDict_for_labels[str(chain_id)][1:])
        frag_spectrum_list.append(f"frag{chain_id}")
    frag_spectrum = " ".join(frag_spectrum_list)

    # color the fragments
    cmd.spectrum(
        "chain", frag_spectrum, selection="ligand", minimum=1, maximum=fragment_num
    )

    # zoom to the ligand
    cmd.zoom("ligand")

    # deselect the ligand
    cmd.deselect()
