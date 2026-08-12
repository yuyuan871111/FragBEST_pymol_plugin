import struct
import tempfile
import unittest
from pathlib import Path

import numpy as np

from FragBEST_pymol_plugin.simple_mesh import Simple_mesh


class SimpleMeshTests(unittest.TestCase):
    def test_loads_ascii_ply_and_triangulates_polygons(self):
        with tempfile.TemporaryDirectory() as directory:
            ply = Path(directory) / "mesh.ply"
            ply.write_text(
                """ply
format ascii 1.0
element vertex 4
property float x
property float y
property float z
property float charge
element face 1
property list uchar int vertex_indices
end_header
0 0 0 -1
1 0 0 0
1 1 0 1
0 1 0 2
4 0 1 2 3
"""
            )

            mesh = Simple_mesh()
            mesh.load_mesh(ply)

            np.testing.assert_allclose(
                mesh.vertices, [[0, 0, 0], [1, 0, 0], [1, 1, 0], [0, 1, 0]]
            )
            np.testing.assert_array_equal(mesh.faces, [[0, 1, 2], [0, 2, 3]])
            np.testing.assert_allclose(
                mesh.get_attribute("vertex_charge"), [-1, 0, 1, 2]
            )

    def test_loads_binary_little_endian_ply(self):
        with tempfile.TemporaryDirectory() as directory:
            ply = Path(directory) / "mesh_binary.ply"
            header = b"""ply
format binary_little_endian 1.0
element vertex 3
property float x
property float y
property float z
property float interest
element face 1
property list uchar int vertex_indices
end_header
"""
            vertices = struct.pack("<12f", 0, 0, 0, 1, 1, 0, 0, 0, 0, 1, 0, 1)
            face = struct.pack("<B3i", 3, 0, 1, 2)
            ply.write_bytes(header + vertices + face)

            mesh = Simple_mesh()
            mesh.load_mesh(ply)

            np.testing.assert_allclose(mesh.vertices, [[0, 0, 0], [1, 0, 0], [0, 1, 0]])
            np.testing.assert_array_equal(mesh.faces, [[0, 1, 2]])
            np.testing.assert_allclose(mesh.get_attribute("vertex_interest"), [1, 0, 1])
