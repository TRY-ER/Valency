import pytest
from engine.explorers.mol_explorer.mol_visualizer import MolVisualizer

def test_get_image_valid_smiles():
    vis = MolVisualizer()
    result = vis.get_image("CCC")
    assert result is not None
    assert isinstance(result, str)

def test_get_image_invalid_smiles():
    vis = MolVisualizer()
    result = vis.get_image("invalid_smiles")
    assert result is None

def test_visualize_valid_smiles():
    vis = MolVisualizer()
    result = vis.visualize("CCO")
    assert result is not None
    assert isinstance(result, str)