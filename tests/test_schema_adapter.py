"""Test schema adapter for workflow editing."""

import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent.parent))

import pytest
import numpy as np
from src.python.generators.workflow.schema_adapter import (
    workflow_to_core,
    core_to_workflow,
    validate_workflow_structure,
    validate_core_structure,
    _cell_matrix_to_params
)


def test_workflow_to_core_simple():
    """Test conversion from workflow to core schema."""
    workflow = {
        "species": ["Si", "Si"],
        "positions": [[0.0, 0.0, 0.0], [2.715, 2.715, 2.715]],
        "cell": [[5.43, 0, 0], [0, 5.43, 0], [0, 0, 5.43]],
        "n_atoms": 2
    }

    core = workflow_to_core(workflow)

    # Check structure
    assert "atoms" in core
    assert "lattice" in core
    assert len(core["atoms"]) == 2

    # Check atoms
    assert core["atoms"][0]["element"] == "Si"
    assert core["atoms"][1]["element"] == "Si"

    # Check fractional coords (should be [0, 0, 0] and [0.5, 0.5, 0.5])
    coords0 = core["atoms"][0]["coords"]
    coords1 = core["atoms"][1]["coords"]
    assert np.allclose(coords0, [0, 0, 0])
    assert np.allclose(coords1, [2.715/5.43, 2.715/5.43, 2.715/5.43], atol=1e-6)

    # Check lattice
    assert core["lattice"]["a"] == pytest.approx(5.43)
    assert core["lattice"]["b"] == pytest.approx(5.43)
    assert core["lattice"]["c"] == pytest.approx(5.43)
    assert core["lattice"]["alpha"] == pytest.approx(90.0)
    assert core["lattice"]["beta"] == pytest.approx(90.0)
    assert core["lattice"]["gamma"] == pytest.approx(90.0)


def test_core_to_workflow_simple():
    """Test conversion from core to workflow schema."""
    core = {
        "atoms": [
            {"element": "Si", "coords": [0.0, 0.0, 0.0]},
            {"element": "Si", "coords": [0.5, 0.5, 0.5]}
        ],
        "lattice": {
            "matrix": [[5.43, 0, 0], [0, 5.43, 0], [0, 0, 5.43]],
            "a": 5.43, "b": 5.43, "c": 5.43,
            "alpha": 90.0, "beta": 90.0, "gamma": 90.0,
            "volume": 160.10
        }
    }

    workflow = core_to_workflow(core)

    # Check structure
    assert "species" in workflow
    assert "positions" in workflow
    assert "cell" in workflow
    assert workflow["n_atoms"] == 2

    # Check species
    assert workflow["species"] == ["Si", "Si"]

    # Check Cartesian positions
    positions = workflow["positions"]
    assert np.allclose(positions[0], [0, 0, 0])
    assert np.allclose(positions[1], [2.715, 2.715, 2.715], atol=1e-6)

    # Check cell
    assert np.allclose(workflow["cell"], [[5.43, 0, 0], [0, 5.43, 0], [0, 0, 5.43]])


def test_roundtrip_conversion():
    """Test that workflow -> core -> workflow preserves data."""
    original_workflow = {
        "species": ["Na", "Cl"],
        "positions": [[0.0, 0.0, 0.0], [2.82, 2.82, 2.82]],
        "cell": [[5.64, 0, 0], [0, 5.64, 0], [0, 0, 5.64]],
        "n_atoms": 2
    }

    # Convert workflow -> core -> workflow
    core = workflow_to_core(original_workflow)
    recovered_workflow = core_to_workflow(core)

    # Check preservation
    assert recovered_workflow["species"] == original_workflow["species"]
    assert recovered_workflow["n_atoms"] == original_workflow["n_atoms"]
    assert np.allclose(recovered_workflow["positions"], original_workflow["positions"], atol=1e-6)
    assert np.allclose(recovered_workflow["cell"], original_workflow["cell"], atol=1e-6)


def test_cell_matrix_to_params_cubic():
    """Test lattice parameter calculation for cubic cell."""
    cell = np.array([[5.0, 0, 0], [0, 5.0, 0], [0, 0, 5.0]])
    params = _cell_matrix_to_params(cell)

    assert params["a"] == pytest.approx(5.0)
    assert params["b"] == pytest.approx(5.0)
    assert params["c"] == pytest.approx(5.0)
    assert params["alpha"] == pytest.approx(90.0)
    assert params["beta"] == pytest.approx(90.0)
    assert params["gamma"] == pytest.approx(90.0)
    assert params["volume"] == pytest.approx(125.0)


def test_cell_matrix_to_params_triclinic():
    """Test lattice parameter calculation for triclinic cell."""
    # Non-orthogonal cell
    cell = np.array([
        [4.0, 0.0, 0.0],
        [1.0, 5.0, 0.0],
        [0.5, 1.0, 6.0]
    ])
    params = _cell_matrix_to_params(cell)

    assert params["a"] == pytest.approx(4.0)
    assert params["b"] == pytest.approx(np.sqrt(1**2 + 5**2))
    assert params["c"] == pytest.approx(np.sqrt(0.5**2 + 1**2 + 6**2))
    assert params["volume"] == pytest.approx(120.0)  # 4 * 5 * 6


def test_validate_workflow_structure_valid():
    """Test validation of valid workflow structure."""
    valid_workflow = {
        "species": ["Si"],
        "positions": [[0, 0, 0]],
        "cell": [[5, 0, 0], [0, 5, 0], [0, 0, 5]]
    }

    is_valid, error = validate_workflow_structure(valid_workflow)
    assert is_valid
    assert error == ""


def test_validate_workflow_structure_missing_field():
    """Test validation fails for missing required field."""
    invalid_workflow = {
        "species": ["Si"],
        "positions": [[0, 0, 0]]
        # Missing "cell"
    }

    is_valid, error = validate_workflow_structure(invalid_workflow)
    assert not is_valid
    assert "cell" in error


def test_validate_workflow_structure_length_mismatch():
    """Test validation fails for species/positions length mismatch."""
    invalid_workflow = {
        "species": ["Si", "Si"],
        "positions": [[0, 0, 0]],  # Only 1 position for 2 species
        "cell": [[5, 0, 0], [0, 5, 0], [0, 0, 5]]
    }

    is_valid, error = validate_workflow_structure(invalid_workflow)
    assert not is_valid
    assert "mismatch" in error.lower()


def test_validate_workflow_structure_singular_cell():
    """Test validation fails for singular cell matrix."""
    invalid_workflow = {
        "species": ["Si"],
        "positions": [[0, 0, 0]],
        "cell": [[1, 0, 0], [2, 0, 0], [3, 0, 0]]  # Singular matrix
    }

    is_valid, error = validate_workflow_structure(invalid_workflow)
    assert not is_valid
    assert "singular" in error.lower()


def test_validate_core_structure_valid():
    """Test validation of valid core structure."""
    valid_core = {
        "atoms": [
            {"element": "Si", "coords": [0, 0, 0]}
        ],
        "lattice": {
            "matrix": [[5, 0, 0], [0, 5, 0], [0, 0, 5]]
        }
    }

    is_valid, error = validate_core_structure(valid_core)
    assert is_valid
    assert error == ""


def test_validate_core_structure_missing_element():
    """Test validation fails for atom missing element."""
    invalid_core = {
        "atoms": [
            {"coords": [0, 0, 0]}  # Missing "element"
        ],
        "lattice": {
            "matrix": [[5, 0, 0], [0, 5, 0], [0, 0, 5]]
        }
    }

    is_valid, error = validate_core_structure(invalid_core)
    assert not is_valid
    assert "element" in error.lower()


def test_non_orthogonal_conversion():
    """Test conversion with non-orthogonal cell."""
    # Hexagonal cell
    a = 5.0
    c = 8.0
    workflow = {
        "species": ["C", "C"],
        "positions": [[0, 0, 0], [a/2, a*np.sqrt(3)/6, c/2]],
        "cell": [
            [a, 0, 0],
            [a/2, a*np.sqrt(3)/2, 0],
            [0, 0, c]
        ],
        "n_atoms": 2
    }

    # Convert to core and back
    core = workflow_to_core(workflow)
    recovered = core_to_workflow(core)

    # Check preservation
    assert recovered["species"] == workflow["species"]
    assert np.allclose(recovered["positions"], workflow["positions"], atol=1e-6)
    assert np.allclose(recovered["cell"], workflow["cell"], atol=1e-6)

    # Check lattice parameters
    assert core["lattice"]["a"] == pytest.approx(a)
    assert core["lattice"]["b"] == pytest.approx(a)
    assert core["lattice"]["c"] == pytest.approx(c)
    assert core["lattice"]["alpha"] == pytest.approx(90.0)
    assert core["lattice"]["beta"] == pytest.approx(90.0)
    # gamma is 60° for cell vectors [a,0,0] and [a/2, a*sqrt(3)/2, 0]
    assert core["lattice"]["gamma"] == pytest.approx(60.0)
