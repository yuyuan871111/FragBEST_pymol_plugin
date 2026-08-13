"""A small, dependency-free PLY reader used by the PyMOL plugin.

This file is modified from masif_pymol_plugin.
Copyright (c) 2019 Gainza P, Sverrisson F, Monti F, Rodola, Bronstein MM,
Correia BE; FragBEST_pymol_plugin Copyright (c) 2024-2026 Yu-Yuan (Stuart)
Yang / Arianna Fornili's Lab. Licensed under the Apache License, Version 2.0.

It reads the scalar vertex properties emitted by MaSIF/FragBEST and triangle
faces. Both ASCII and standard binary little-/big-endian PLY files are
supported.
"""

from pathlib import Path

import numpy as np

_PLY_DTYPES = {
    "char": "i1",
    "int8": "i1",
    "uchar": "u1",
    "uint8": "u1",
    "short": "i2",
    "int16": "i2",
    "ushort": "u2",
    "uint16": "u2",
    "int": "i4",
    "int32": "i4",
    "uint": "u4",
    "uint32": "u4",
    "float": "f4",
    "float32": "f4",
    "double": "f8",
    "float64": "f8",
}


class Simple_mesh:
    """Load PLY meshes without requiring PyMesh or another external package."""

    def __init__(self):
        self.vertices = np.empty((0, 3), dtype=np.float32)
        self.faces = np.empty((0, 3), dtype=np.int32)
        self.attributes = {}
        self.attribute_names = []

    @staticmethod
    def _dtype(type_name, byte_order="="):
        try:
            return np.dtype(byte_order + _PLY_DTYPES[type_name])
        except KeyError as exc:
            raise ValueError(f"Unsupported PLY property type: {type_name}") from exc

    @staticmethod
    def _read_header(handle):
        first_line = handle.readline()
        if first_line.strip() != b"ply":
            raise ValueError("Not a PLY file (missing 'ply' header).")

        fmt = None
        elements = []
        current = None
        while True:
            raw_line = handle.readline()
            if not raw_line:
                raise ValueError("PLY header ended before 'end_header'.")
            fields = raw_line.decode("ascii").strip().split()
            if not fields or fields[0] in {"comment", "obj_info"}:
                continue
            if fields[0] == "format":
                fmt = fields[1]
            elif fields[0] == "element":
                current = {"name": fields[1], "count": int(fields[2]), "properties": []}
                elements.append(current)
            elif fields[0] == "property":
                if current is None:
                    raise ValueError(
                        "PLY property found before an element declaration."
                    )
                if fields[1] == "list":
                    current["properties"].append(
                        ("list", fields[2], fields[3], fields[4])
                    )
                else:
                    current["properties"].append(("scalar", fields[1], fields[2]))
            elif fields[0] == "end_header":
                break

        if fmt not in {"ascii", "binary_little_endian", "binary_big_endian"}:
            raise ValueError(f"Unsupported PLY format: {fmt!r}")
        return fmt, elements

    def load_mesh(self, filename):
        """Read mesh data in one pass, avoiding per-vertex Python work."""
        filename = Path(filename)
        with filename.open("rb") as handle:
            fmt, elements = self._read_header(handle)
            vertex = next((item for item in elements if item["name"] == "vertex"), None)
            face = next((item for item in elements if item["name"] == "face"), None)
            if vertex is None:
                raise ValueError("PLY file has no vertex element.")
            if face is None:
                raise ValueError("PLY file has no face element.")
            if any(prop[0] != "scalar" for prop in vertex["properties"]):
                raise ValueError("List-valued vertex properties are not supported.")

            self.num_verts = vertex["count"]
            self.num_faces = face["count"]
            self.attribute_names = [
                "vertex_" + prop[2] for prop in vertex["properties"]
            ]
            scalar_names = [prop[2] for prop in vertex["properties"]]

            if fmt == "ascii":
                vertex_data = np.loadtxt(
                    handle, dtype=np.float32, max_rows=self.num_verts, ndmin=2
                )
                faces = self._read_ascii_faces(handle, self.num_faces)
            else:
                order = "<" if fmt == "binary_little_endian" else ">"
                vertex_dtype = np.dtype(
                    [
                        (name, self._dtype(type_name, order))
                        for _, type_name, name in vertex["properties"]
                    ]
                )
                vertex_data = np.fromfile(
                    handle, dtype=vertex_dtype, count=self.num_verts
                )
                faces = self._read_binary_faces(handle, face, order)

        if vertex_data.shape[0] != self.num_verts:
            raise ValueError("PLY file ended before all vertices were read.")
        self.attributes = {
            "vertex_"
            + name: np.asarray(
                vertex_data[:, index] if fmt == "ascii" else vertex_data[name]
            )
            for index, name in enumerate(scalar_names)
        }
        try:
            self.vertices = np.column_stack(
                (
                    self.attributes["vertex_x"],
                    self.attributes["vertex_y"],
                    self.attributes["vertex_z"],
                )
            )
        except KeyError as exc:
            raise ValueError("PLY vertex properties must include x, y, and z.") from exc
        self.faces = faces

    @staticmethod
    def _read_ascii_faces(handle, count):
        faces = []
        for _ in range(count):
            fields = handle.readline().split()
            if not fields:
                raise ValueError("PLY file ended before all faces were read.")
            size = int(fields[0])
            if len(fields) < size + 1:
                raise ValueError("Malformed PLY face record.")
            faces.extend(
                Simple_mesh._triangulate([int(value) for value in fields[1 : size + 1]])
            )
        return np.asarray(faces, dtype=np.int32).reshape((-1, 3))

    def _read_binary_faces(self, handle, face, order):
        if not face["properties"] or face["properties"][0][0] != "list":
            raise ValueError(
                "PLY face element must start with a list of vertex indices."
            )
        _, count_type, index_type, _ = face["properties"][0]
        count_dtype, index_dtype = self._dtype(count_type, order), self._dtype(
            index_type, order
        )
        faces = []
        for _ in range(self.num_faces):
            size_data = np.fromfile(handle, dtype=count_dtype, count=1)
            if not size_data.size:
                raise ValueError("PLY file ended before all faces were read.")
            indices = np.fromfile(handle, dtype=index_dtype, count=int(size_data[0]))
            if len(indices) != int(size_data[0]):
                raise ValueError("Malformed binary PLY face record.")
            faces.extend(self._triangulate(indices))
        return np.asarray(faces, dtype=np.int32).reshape((-1, 3))

    @staticmethod
    def _triangulate(indices):
        if len(indices) < 3:
            return []
        return [
            (indices[0], indices[index], indices[index + 1])
            for index in range(1, len(indices) - 1)
        ]

    def get_attribute_names(self):
        return list(self.attribute_names)

    def get_attribute(self, attribute_name):
        return self.attributes[attribute_name]
