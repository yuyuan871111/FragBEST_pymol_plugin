"""
    loadPREDICT.py: This pymol function loads prediction files into pymol.
    It is modified by Yu-Yuan (Stuart) Yang - Arianna group QMUL 2024-2027
    This is part of probinsense_pymol_plugin.
"""

from pymol import cmd
from pymol.cgo import SPHERE

from .loadPLY import label_color


def get_from(obj_name, what="b"):
    """
    Get the values of a property from an object.

    Input:
        - `obj_name`: str
        - `where`: str (`b`: b-factor; `q`: occupancy)
    """
    values = []
    cmd.iterate_state(1, obj_name, f"values.append({what})", space=locals(), atomic=0)
    return values


def draw_points(verts, label, dotSizes, name="_grid", scaler=0.3):
    # Plot points
    obj = []
    color_array_surf = label_color(label)
    for v_ix in range(len(verts)):
        vert = verts[v_ix]
        dotSize = dotSizes[v_ix]
        colorToAdd = color_array_surf[v_ix]
        obj.extend(colorToAdd)
        obj.extend([SPHERE, vert[0], vert[1], vert[2], dotSize * scaler])

    cmd.load_cgo(obj, name, 1.0)


def load_grids(filename, scaler=0.2):
    # get the object name
    obj_name = f'raw_{filename.replace(".pdb", "")}'

    # load the object
    cmd.load(filename, obj_name)
    xyz = cmd.get_coordset(obj_name)
    b_factors = get_from(obj_name, what="b")
    probs = get_from(obj_name, what="q")

    # draw the points
    grid_name = obj_name.replace("raw_", "grids_")
    draw_points(xyz, b_factors, probs, name=grid_name, scaler=scaler)

    # group the objects
    group_names = f"{obj_name} {grid_name}"
    cmd.group(filename, group_names)
