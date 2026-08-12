import unittest

import numpy as np

from FragBEST_pymol_plugin.color_palette import colorDict
from FragBEST_pymol_plugin.loadPLY import gradient_label_color, label_color


class LabelColorTests(unittest.TestCase):
    def test_background_is_gray_and_labels_form_a_rainbow(self):
        colors = np.asarray(gradient_label_color(np.arange(13)))

        np.testing.assert_allclose(colors[0, 1:], colorDict["gray"][1:])
        np.testing.assert_allclose(colors[1, 1:], [1.0, 0.0, 0.0])
        np.testing.assert_allclose(colors[-1, 1:], [0.68, 0.0, 1.0])
        self.assertEqual(len({tuple(color[1:]) for color in colors[1:]}), 12)

    def test_legacy_palette_remains_available(self):
        colors = np.asarray(label_color([0, 1, 10]))

        np.testing.assert_allclose(colors[0, 1:], colorDict["gray"][1:])
        np.testing.assert_allclose(colors[1, 1:], colorDict["blue"][1:])
        np.testing.assert_allclose(colors[2, 1:], colorDict["greentint"][1:])
