import unittest

import numpy as np
from pymol import cmd

from FragBEST_pymol_plugin.loadPLY import pymol_object_name, show_vertex_text_labels


class VertexLabelTests(unittest.TestCase):
    def tearDown(self):
        cmd.delete("test_vertex_labels")

    def test_creates_screen_facing_integer_labels_at_offset_positions(self):
        name = show_vertex_text_labels(
            np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]]),
            np.array([1.9, 2.1]),
            "test_vertex_labels",
            interest=np.array([0, 1]),
            normals=np.array([[0.0, 0.0, 1.0], [0.0, 0.0, 1.0]]),
            offset=0.5,
        )

        self.assertEqual(name, "test_vertex_labels")
        self.assertEqual(cmd.count_atoms(name), 1)
        np.testing.assert_allclose(cmd.get_atom_coords(name), [1.0, 0.0, 0.5])
        labels = []
        cmd.iterate(name, "labels.append(label)", space={"labels": labels})
        self.assertEqual(labels, ["2"])

    def test_sanitizes_file_paths_for_pymol_object_names(self):
        self.assertEqual(pymol_object_name("/tmp/example mesh.ply"), "example_mesh_ply")
