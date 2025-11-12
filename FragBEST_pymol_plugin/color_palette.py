"""
Yu-Yuan Yang 2024
This is a color palette for PyMOL.
"""

from pymol.cgo import COLOR

colorDict = {
    "sky": [COLOR, 0.0, 0.76, 1.0],  # 00C3FF
    "sea": [COLOR, 0.0, 0.90, 0.5],  # 00E67F
    "yellowtint": [COLOR, 0.88, 0.97, 0.02],  # E0F700
    "hotpink": [COLOR, 0.90, 0.40, 0.70],  # FF69B4
    "greentint": [COLOR, 0.50, 0.90, 0.40],  # 80E6CC
    "blue": [COLOR, 0.0, 0.0, 1.0],  # 0000FF
    "green": [COLOR, 0.0, 1.0, 0.0],  # 00FF00
    "yellow": [COLOR, 1.0, 1.0, 0.0],  # FFFF00
    "orange": [COLOR, 1.0, 0.5, 0.0],  # FF8000
    "red": [COLOR, 1.0, 0.0, 0.0],  # FF0000
    "black": [COLOR, 0.0, 0.0, 0.0],  # 000000
    "white": [COLOR, 1.0, 1.0, 1.0],  # FFFFFF
    "gray": [COLOR, 0.9, 0.9, 0.9],  # E6E6E6
}

colorDict_for_labels = {
    "0": colorDict["gray"],
    "1": colorDict["blue"],
    "2": colorDict["green"],
    "3": colorDict["yellow"],
    "4": colorDict["orange"],
    "5": colorDict["red"],
    "6": colorDict["sky"],
    "7": colorDict["sea"],
    "8": colorDict["yellowtint"],
    "9": colorDict["hotpink"],
    "10": colorDict["greentint"],
    # For the comparison of the prediction and the ground truth
    "correct": colorDict["green"],
    "incorrect": colorDict["red"],
}
