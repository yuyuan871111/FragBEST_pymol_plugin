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
loadPLY.py: This pymol function loads ply files into pymol.
Pablo Gainza - LPDI STI EPFL 2016-2019
This file is part of MaSIF.
Released under an Apache License 2.0
"""

import re
from pathlib import Path
from typing import Callable, Union

import numpy as np
from pymol import cmd
from pymol.cgo import (
    BEGIN,
    COLOR,
    END,
    LINES,
    NORMAL,
    SPHERE,
    TRIANGLES,
    VERTEX,
)

from .color_palette import colorDict, colorDict_for_labels


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


def gradient_label_color(label):
    """Color class 0 gray and all other classes along a red-to-purple rainbow.

    The gradient is assigned to the sorted, non-zero class IDs present in the
    mesh. This supports an arbitrary number of classes while keeping adjacent
    class IDs visually ordered. Class 0 is deliberately excluded from the
    rainbow because it represents the background.
    """
    labels = np.asarray(label, dtype=int)
    colors = np.empty((len(labels), 3), dtype=float)
    background = labels == 0
    colors[background] = colorDict["gray"][1:]

    classes, class_indices = np.unique(labels[~background], return_inverse=True)
    if len(classes):
        # HSV hue 0.0 is red; 0.78 is a purple/magenta. Traversing the hue
        # range gives red, yellow, green, cyan, blue, and finally purple.
        hues = np.linspace(0.0, 0.78, len(classes))
        hue = hues[class_indices]
        sector = np.floor(hue * 6).astype(int) % 6
        fraction = hue * 6 - np.floor(hue * 6)
        rainbow = np.empty((len(hue), 3), dtype=float)
        for value in range(6):
            mask = sector == value
            fraction_at_sector = fraction[mask]
            zeros = np.zeros_like(fraction_at_sector)
            ones = np.ones_like(fraction_at_sector)
            if value == 0:
                rainbow[mask] = np.column_stack((ones, fraction_at_sector, zeros))
            elif value == 1:
                rainbow[mask] = np.column_stack(
                    (ones - fraction_at_sector, ones, zeros)
                )
            elif value == 2:
                rainbow[mask] = np.column_stack((zeros, ones, fraction_at_sector))
            elif value == 3:
                rainbow[mask] = np.column_stack(
                    (zeros, ones - fraction_at_sector, ones)
                )
            elif value == 4:
                rainbow[mask] = np.column_stack((fraction_at_sector, zeros, ones))
            else:
                rainbow[mask] = np.column_stack(
                    (ones, zeros, ones - fraction_at_sector)
                )
        colors[~background] = rainbow

    return [[COLOR, red, green, blue] for red, green, blue in colors]


def label_color(label):
    """Legacy fixed palette for labels 0 through 10."""
    return [colorDict_for_labels[str(int(value))] for value in label]


def pymol_object_name(name):
    """Return a PyMOL-safe object name without path separators or punctuation."""
    safe_name = re.sub(r"[^A-Za-z0-9_]", "_", Path(str(name)).name)
    return safe_name.strip("_") or "ply"


def show_vertex_text_labels(
    verts, values, name, interest=None, normals=None, offset=0.3
):
    """Show screen-facing integer labels offset from selected mesh vertices."""
    indices = (
        np.flatnonzero(np.asarray(interest) != 0)
        if interest is not None
        else np.arange(len(verts))
    )
    positions = np.asarray(verts, dtype=float)[indices].copy()
    directions = np.tile(np.array([1.0, 1.0, 1.0]) / np.sqrt(3), (len(indices), 1))
    if normals is not None:
        candidate_directions = np.asarray(normals, dtype=float)[indices]
        lengths = np.linalg.norm(candidate_directions, axis=1)
        valid = lengths > 0
        directions[valid] = candidate_directions[valid] / lengths[valid, None]
    positions += float(offset) * directions

    cmd.delete(name)
    for position, value in zip(positions, np.asarray(values)[indices]):
        cmd.pseudoatom(
            name,
            pos=position.tolist(),
            label=str(int(value)),
            quiet=1,
        )
    cmd.hide("everything", name)
    cmd.show("labels", name)
    cmd.set("label_color", "white", name)
    cmd.set("label_outline_color", "black", name)
    cmd.set("label_size", 14, name)
    return name


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
    if where == "vertex":
        # Select first: an interest-point view should not spend time creating
        # colors or CGO commands for the (usually much larger) background.
        indices = (
            np.flatnonzero(np.asarray(interest) != 0)
            if interest is not None
            else np.arange(len(verts))
        )
        selected_data = np.asarray(data)[indices]
        selected_verts = np.asarray(verts)[indices]
        color_array = color_style(selected_data)
        selected_sizes = None if np.isscalar(dotSize) else np.asarray(dotSize)[indices]

        obj = []
        for position, (vert, color_to_add) in enumerate(
            zip(selected_verts, color_array)
        ):
            each_dotSize = (
                dotSize if selected_sizes is None else selected_sizes[position]
            )
            obj.extend(color_to_add)
            obj.extend([SPHERE, vert[0], vert[1], vert[2], each_dotSize])
        return obj

    if where == "surface":
        obj = []
        color_array_surf = color_style(data)
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

        return obj

    raise ValueError("where should be 'vertex' or 'surface'")


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
    if kwargs.get("normals") is not None and kwargs["ignore_surface"] == 0:
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
    label_palette: str = "gradient",
    show_vertex_labels: int = 0,
    custom_name: str = None,
    dotSize: float = 0.2,
    text_offset: float = 0.3,
):
    """
    Load a PLY file into PyMOL.

    ``label_palette`` may be ``"gradient"`` (the default, supports any
    number of labels) or ``"legacy"`` (the original fixed palette for labels
    0 through 10). Set integer ``show_vertex_labels`` to ``1`` to display
    available ``vertex_label`` and ``vertex_pred`` class numbers next to each
    selected vertex.
    """
    # Check
    ignore_surface = int(ignore_surface)
    interest_pt = int(interest_pt)
    show_vertex_labels = int(show_vertex_labels)
    dotSize = float(dotSize)
    text_offset = float(text_offset)
    label_palette = str(label_palette).lower()
    assert ignore_surface == 0 or ignore_surface == 1, "ignore_surface should be 0 or 1"
    assert interest_pt == 0 or interest_pt == 1, "interest_pt should be 0 or 1"
    assert show_vertex_labels in {0, 1}, "show_vertex_labels should be 0 or 1"
    assert text_offset >= 0, "text_offset should be non-negative"
    if label_palette not in {"gradient", "legacy"}:
        raise ValueError("label_palette should be 'gradient' or 'legacy'")
    label_color_style = (
        gradient_label_color if label_palette == "gradient" else label_color
    )

    # Load the mesh
    ## Pymesh should be faster and supports binary ply files. However it is difficult to install with pymol...
    #        import pymesh
    #        mesh = pymesh.load_mesh(filename)
    from .simple_mesh import Simple_mesh

    mesh = Simple_mesh()
    mesh.load_mesh(filename)
    verts = mesh.vertices
    faces = mesh.faces

    normals = None
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
        vertex_labels = mesh.get_attribute("vertex_label")
        label = {
            "name": "label",
            "data": vertex_labels,
            "verts": verts,
            "faces": faces,
            "normals": normals,
            "color_style": label_color_style,
            "interest": interest,
            "ignore_surface": ignore_surface,
            "dotSize": dotSize,
            "custom_name": custom_name,
        }
        # Draw features
        objs = draw_features(**label)
        obj_list.extend(objs)

        if show_vertex_labels:
            name = "labelname_vert_" + pymol_object_name(custom_name)
            show_vertex_text_labels(
                verts,
                vertex_labels,
                name,
                interest=interest,
                normals=normals,
                offset=text_offset,
            )
            obj_list.append(name)

    # Draw prediction
    if "vertex_pred" in mesh.get_attribute_names():
        vertex_predictions = mesh.get_attribute("vertex_pred")
        pred = {
            "name": "pred",
            "data": vertex_predictions,
            "verts": verts,
            "faces": faces,
            "normals": normals,
            "color_style": label_color_style,
            "interest": interest,
            "ignore_surface": ignore_surface,
            "dotSize": 0.2 * mesh.get_attribute("vertex_predprobs"),
            "custom_name": custom_name,
        }
        # Draw features
        objs = draw_features(**pred)
        obj_list.extend(objs)

        if show_vertex_labels:
            name = "predname_vert_" + pymol_object_name(custom_name)
            show_vertex_text_labels(
                verts,
                vertex_predictions,
                name,
                interest=interest,
                normals=normals,
                offset=text_offset,
            )
            obj_list.append(name)

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
    cmd.group(custom_name, group_names)
