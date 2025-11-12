"""
    loadPLY.py: This pymol function loads ply files into pymol.
    It is modified by Yu-Yuan (Stuart) Yang - Arianna group QMUL 2024-2027
    This is part of probinsense_pymol_plugin.

    Pablo Gainza - LPDI STI EPFL 2016-2019
    This file is part of MaSIF.
    Released under an Apache License 2.0
"""

from pathlib import Path
from typing import Callable, Union

import numpy as np
from pymol import cmd
from pymol.cgo import BEGIN, COLOR, END, LINES, NORMAL, SPHERE, TRIANGLES, VERTEX

from .color_palette import colorDict, colorDict_for_labels

# def color_gradient(vals, color1, color2):
#     """
#     Create a gradient color from color 1 to whitish, to color 2. val goes from 0 (color1) to 1 (color2).
#     """
#     c1 = Color("white")
#     c2 = Color("orange")
#     ix = np.floor(vals * 100).astype(int)
#     crange = list(c1.range_to(c2, 100))
#     mycolor = []
#     print(crange[0].get_rgb())
#     for x in ix:
#         myc = crange[x].get_rgb()
#         mycolor.append([COLOR, myc[0], myc[1], myc[2]])
#     return mycolor


def iface_color(iface):
    # max value is 1, min values is 0
    hp = iface.copy()
    hp = hp * 2 - 1
    mycolor = charge_color(-hp)
    return mycolor


def apbs_color(apbs):
    # max value is 3, min values is -3
    apbs = apbs / 3
    return charge_color(apbs)


def hphob_color(hphob):
    """
    Returns the color of each vertex according to the charge.
    The most purple colors are the most hydrophilic values, and the most
    white colors are the most positive colors.
    """
    # max value is 4.5, min values is -4.5
    hp = hphob.copy()
    # normalize
    hp = hp + 4.5
    hp = hp / 9.0
    # mycolor = [ [COLOR, 1.0, hp[i], 1.0]  for i in range(len(hp)) ]
    mycolor = [[COLOR, 1.0, 1.0 - hp[i], 1.0] for i in range(len(hp))]
    return mycolor


def label_color(label):
    mycolor = []
    for i in range(len(label)):
        mycolor.append(colorDict_for_labels[str(int(label[i]))])
    return mycolor


def true_false_label_color(label):
    mycolor = []
    for i in range(len(label)):
        if label[i] == 1:
            mycolor.append(colorDict_for_labels["correct"])
        else:
            mycolor.append(colorDict_for_labels["incorrect"])
    return mycolor


def single_color(data, color="green"):
    assert color in colorDict.keys(), "Color not found in colorDict"
    return [colorDict[color]] * len(data)


def charge_color(charges):
    """
    Returns the color of each vertex according to the charge.
    The most red colors are the most negative values, and the most
    blue colors are the most positive colors.
    """
    # Assume a std deviation equal for all proteins....
    max_val = 1.0
    min_val = -1.0

    norm_charges = charges
    blue_charges = np.array(norm_charges)
    red_charges = np.array(norm_charges)
    blue_charges[blue_charges < 0] = 0
    red_charges[red_charges > 0] = 0
    red_charges = abs(red_charges)
    red_charges[red_charges > max_val] = max_val
    blue_charges[blue_charges < min_val] = min_val
    red_charges = red_charges / max_val
    blue_charges = blue_charges / max_val
    # red_charges[red_charges>1.0] = 1.0
    # blue_charges[blue_charges>1.0] = 1.0
    # green_color = np.array([0.0] * len(charges))
    mycolor = [
        [
            COLOR,
            0.9999 - blue_charges[i],
            0.9999 - (blue_charges[i] + red_charges[i]),
            0.9999 - red_charges[i],
        ]
        for i in range(len(charges))
    ]
    for i in range(len(mycolor)):
        for k in range(1, 4):
            if mycolor[i][k] < 0:
                mycolor[i][k] = 0

    return mycolor


def draw_on(
    data,
    verts,
    color_style: Callable,
    interest=None,
    faces=None,
    normals=None,
    where: str = "vertex",
    dotSize: float = 0.2,
):
    """
    Draw data on the mesh.

    Intput:
        - `data`: list of values to draw
        - `verts`: list of vertices
        - `color_style`: function that returns a color array
        - `interest`: list of interest values
        - `faces`: list of faces
        - `normals`: list of normals
        - `where`: ['surface', 'vertex']
        - `dotSize`: size of the dots

    Output:
        - `obj`: list of cgo commands
    """
    obj = []
    color_array_surf = color_style(data)

    if where == "vertex":
        # Check if dotSize is a float or a list
        if type(dotSize) is float or type(dotSize) is int:
            dotSize = [dotSize] * len(verts)
            dotSize = np.array(dotSize)

        # Plot vertices
        for v_ix, each_dotSize in enumerate(dotSize):
            vert = verts[v_ix]
            colorToAdd = color_array_surf[v_ix]

            if interest is not None and interest[v_ix] == 0:
                continue
            obj.extend(colorToAdd)
            obj.extend([SPHERE, vert[0], vert[1], vert[2], each_dotSize])

    elif where == "surface":
        # Plot faces
        for tri in faces:
            vert1 = verts[int(tri[0])]
            vert2 = verts[int(tri[1])]
            vert3 = verts[int(tri[2])]
            na = normals[int(tri[0])]
            nb = normals[int(tri[1])]
            nc = normals[int(tri[2])]
            if (
                interest is not None
                and (interest[int(tri[0])] == 0)
                and (interest[int(tri[1])] == 0)
                and (interest[int(tri[2])] == 0)
            ):
                continue
            obj.extend([BEGIN, TRIANGLES])
            # obj.extend([ALPHA, 0.6])
            obj.extend(color_array_surf[int(tri[0])])
            obj.extend([NORMAL, (na[0]), (na[1]), (na[2])])
            obj.extend([VERTEX, (vert1[0]), (vert1[1]), (vert1[2])])
            obj.extend(color_array_surf[int(tri[1])])
            obj.extend([NORMAL, (nb[0]), (nb[1]), (nb[2])])
            obj.extend([VERTEX, (vert2[0]), (vert2[1]), (vert2[2])])
            obj.extend(color_array_surf[int(tri[2])])
            obj.extend([NORMAL, (nc[0]), (nc[1]), (nc[2])])
            obj.extend([VERTEX, (vert3[0]), (vert3[1]), (vert3[2])])
            obj.append(END)

    else:
        raise ValueError("where should be 'vertex' or 'surface'")

    return obj


def draw_mesh(verts, faces, interest=None):
    obj = []
    # Plot mesh
    for tri in faces:
        pairs = [[tri[0], tri[1]], [tri[0], tri[2]], [tri[1], tri[2]]]
        colorToAdd = colorDict["gray"]
        for pair in pairs:
            vert1 = verts[pair[0]]
            vert2 = verts[pair[1]]
            if (
                interest is not None
                and (interest[pair[0]] == 0)
                and (interest[pair[1]] == 0)
            ):
                continue
            obj.extend([BEGIN, LINES])
            obj.extend(colorToAdd)
            obj.extend([VERTEX, (vert1[0]), (vert1[1]), (vert1[2])])
            obj.extend([VERTEX, (vert2[0]), (vert2[1]), (vert2[2])])
            obj.append(END)
    return obj


def draw_features(**kwargs):
    """
    Draw features on the mesh.
    """
    obj_list = []
    # Draw on vertices
    obj = draw_on(
        data=kwargs["data"],
        verts=kwargs["verts"],
        color_style=kwargs["color_style"],
        interest=kwargs["interest"],
        dotSize=kwargs["dotSize"],
        where="vertex",
    )
    name = kwargs["name"] + "_vert_" + kwargs["custom_name"]
    cmd.load_cgo(obj, name, 1.0)
    obj_list.append(name)

    # Draw on surface
    if "normals" in kwargs and kwargs["ignore_surface"] == 0:
        obj = draw_on(
            data=kwargs["data"],
            verts=kwargs["verts"],
            color_style=kwargs["color_style"],
            faces=kwargs["faces"],
            normals=kwargs["normals"],
            interest=kwargs["interest"],
            where="surface",
        )
        name = kwargs["name"] + "_surf_" + kwargs["custom_name"]
        cmd.load_cgo(obj, name, 1.0)
        obj_list.append(name)

    return obj_list


def load_ply(
    filename: Union[str, Path],
    interest_pt: int = 1,
    ignore_surface: int = 0,
    custom_name: str = None,
    dotSize: float = 0.2,
):
    """
    Main funcion to load a ply file into pymol.
    """
    # Check
    ignore_surface = int(ignore_surface)
    interest_pt = int(interest_pt)
    dotSize = float(dotSize)
    assert ignore_surface == 0 or ignore_surface == 1, "ignore_surface should be 0 or 1"
    assert interest_pt == 0 or interest_pt == 1, "interest_pt should be 0 or 1"

    # Load the mesh
    ## Pymesh should be faster and supports binary ply files. However it is difficult to install with pymol...
    #        import pymesh
    #        mesh = pymesh.load_mesh(filename)
    from .simple_mesh import Simple_mesh

    mesh = Simple_mesh()
    mesh.load_mesh(filename)
    verts = mesh.vertices
    faces = mesh.faces

    if "vertex_nx" in mesh.get_attribute_names():
        nx = mesh.get_attribute("vertex_nx")
        ny = mesh.get_attribute("vertex_ny")
        nz = mesh.get_attribute("vertex_nz")
        normals = np.vstack([nx, ny, nz]).T
        # print(normals.shape)

    # Read interest points
    print(f"interest_pt: {interest_pt}")
    if interest_pt != 0 and "vertex_interest" in mesh.get_attribute_names():
        interest = mesh.get_attribute("vertex_interest")
    else:
        interest = None
    print(f"ignore_surface: {ignore_surface}")

    # Initialisation
    obj_list = []
    if custom_name is None:
        custom_name = filename

    # Draw APBS charges
    if "vertex_charge" in mesh.get_attribute_names():
        apbs = {
            "name": "pb",
            "data": mesh.get_attribute("vertex_charge"),
            "verts": verts,
            "faces": faces,
            "normals": normals,
            "color_style": apbs_color,
            "interest": interest,
            "ignore_surface": ignore_surface,
            "dotSize": dotSize,
            "custom_name": custom_name,
        }
        # Draw features
        objs = draw_features(**apbs)
        obj_list.extend(objs)

    # Draw hydrophobicity
    if "vertex_hphob" in mesh.get_attribute_names():
        hphob = {
            "name": "hphob",
            "data": mesh.get_attribute("vertex_hphob"),
            "verts": verts,
            "faces": faces,
            "normals": normals,
            "color_style": hphob_color,
            "interest": interest,
            "ignore_surface": ignore_surface,
            "dotSize": dotSize,
            "custom_name": custom_name,
        }
        # Draw features
        objs = draw_features(**hphob)
        obj_list.extend(objs)

    # Draw hbond
    if "vertex_hbond" in mesh.get_attribute_names():
        hbond = {
            "name": "hbond",
            "data": mesh.get_attribute("vertex_hbond"),
            "verts": verts,
            "faces": faces,
            "normals": normals,
            "color_style": charge_color,
            "interest": interest,
            "ignore_surface": ignore_surface,
            "dotSize": dotSize,
            "custom_name": custom_name,
        }
        # Draw features
        objs = draw_features(**hbond)
        obj_list.extend(objs)

    # Draw label
    if "vertex_label" in mesh.get_attribute_names():
        label = {
            "name": "label",
            "data": mesh.get_attribute("vertex_label"),
            "verts": verts,
            "faces": faces,
            "normals": normals,
            "color_style": label_color,
            "interest": interest,
            "ignore_surface": ignore_surface,
            "dotSize": dotSize,
            "custom_name": custom_name,
        }
        # Draw features
        objs = draw_features(**label)
        obj_list.extend(objs)

    # Draw prediction
    if "vertex_pred" in mesh.get_attribute_names():
        pred = {
            "name": "pred",
            "data": mesh.get_attribute("vertex_pred"),
            "verts": verts,
            "faces": faces,
            "normals": normals,
            "color_style": label_color,
            "interest": interest,
            "ignore_surface": ignore_surface,
            "dotSize": 0.2 * mesh.get_attribute("vertex_predprobs"),
            "custom_name": custom_name,
        }
        # Draw features
        objs = draw_features(**pred)
        obj_list.extend(objs)

    # Draw comparison of the prediction and the label
    if (
        "vertex_pred" in mesh.get_attribute_names()
        and "vertex_label" in mesh.get_attribute_names()
    ):
        comp_pred_label = [
            1 if _pred == _label else 0
            for _pred, _label in zip(
                mesh.get_attribute("vertex_pred"), mesh.get_attribute("vertex_label")
            )
        ]
        comp = {
            "name": "compare",
            "data": comp_pred_label,
            "verts": verts,
            "faces": faces,
            "normals": normals,
            "color_style": true_false_label_color,
            "interest": interest,
            "ignore_surface": ignore_surface,
            "dotSize": 0.2 * mesh.get_attribute("vertex_predprobs"),
            "custom_name": custom_name,
        }
        # Draw features
        objs = draw_features(**comp)
        obj_list.extend(objs)

    # Draw triangles (faces)
    obj = draw_mesh(verts=verts, faces=faces, interest=interest)
    name = "mesh_" + custom_name
    cmd.load_cgo(obj, name, 1.0)
    obj_list.append(name)

    # Grouping all objects
    group_names = " ".join(obj_list)
    print(group_names)
    cmd.group(filename, group_names)
