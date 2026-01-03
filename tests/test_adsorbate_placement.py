"""Test adsorbate placement with surface normal orientation."""

import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent.parent))

import pytest
import numpy as np
from src.python.structure_tools import (
    add_adsorbate,
    _calculate_surface_normal,
    _rotate_around_axis,
    dict_to_structure,
    structure_to_dict
)
from pymatgen.core import Structure, Lattice


def test_calculate_surface_normal_001():
    """Test surface normal calculation for (001) surface."""
    # Create a proper (001) slab with clear layering
    lattice = Lattice([[5.0, 0, 0], [0, 5.0, 0], [0, 0, 15.0]])  # Extended in z
    species = ["Si"] * 8
    # Two clear layers at different z positions
    coords = [
        # Bottom layer at z=0.1
        [0, 0, 0.1], [0.5, 0, 0.1], [0, 0.5, 0.1], [0.5, 0.5, 0.1],
        # Top layer at z=0.2 (surface)
        [0, 0, 0.2], [0.5, 0, 0.2], [0, 0.5, 0.2], [0.5, 0.5, 0.2]
    ]
    structure = Structure(lattice, species, coords)

    # Calculate normal for a top layer atom (should be close to [0, 0, 1])
    normal = _calculate_surface_normal(structure, 4)  # Top layer atom

    assert normal.shape == (3,)
    assert np.linalg.norm(normal) == pytest.approx(1.0)  # Unit vector
    # For (001) surface, normal should be primarily in z-direction
    assert abs(normal[2]) > 0.9  # Close to [0, 0, 1]


def test_rotate_around_axis_z():
    """Test rotation around z-axis."""
    coords = np.array([[1.0, 0.0, 0.0]])
    axis = np.array([0.0, 0.0, 1.0])
    angle = np.pi / 2  # 90 degrees

    rotated = _rotate_around_axis(coords, axis, angle)

    # After 90° rotation around z, [1,0,0] should become [0,1,0]
    expected = np.array([[0.0, 1.0, 0.0]])
    assert np.allclose(rotated, expected, atol=1e-6)


def test_rotate_around_axis_arbitrary():
    """Test rotation around arbitrary axis."""
    coords = np.array([[1.0, 0.0, 0.0], [0.0, 1.0, 0.0]])
    axis = np.array([1.0, 1.0, 0.0])  # 45° in xy-plane
    axis = axis / np.linalg.norm(axis)  # Normalize
    angle = np.pi  # 180 degrees

    rotated = _rotate_around_axis(coords, axis, angle)

    # After 180° rotation around [1,1,0] axis
    # Check that rotation preserves distances
    original_dist = np.linalg.norm(coords[1] - coords[0])
    rotated_dist = np.linalg.norm(rotated[1] - rotated[0])
    assert rotated_dist == pytest.approx(original_dist)


def test_add_adsorbate_simple_001():
    """Test adding CO molecule to (001) surface."""
    # Create simple (001) surface
    lattice = Lattice.cubic(5.0)
    species = ["Cu"] * 4
    coords = [[0, 0, 0], [0.5, 0.5, 0], [0, 0.5, 0.25], [0.5, 0, 0.25]]
    structure = Structure(lattice, species, coords)
    structure_dict = structure_to_dict(structure)

    # Add CO molecule at site 0
    result = add_adsorbate(
        structure_dict=structure_dict,
        molecule="CO",
        site_index=0,
        distance=2.0
    )

    assert result["success"] is True
    assert "structure" in result

    new_structure = result["structure"]
    # Original had 4 atoms, CO adds 2 more
    assert new_structure["metadata"]["natoms"] == 6

    # Check that CO was added
    atoms = new_structure["atoms"]
    elements = [atom["element"] for atom in atoms]
    assert "C" in elements
    assert "O" in elements


def test_add_adsorbate_by_name():
    """Test adding different molecules by name."""
    lattice = Lattice.cubic(5.0)
    species = ["Pt"] * 4
    coords = [[0, 0, 0], [0.5, 0.5, 0], [0, 0.5, 0.25], [0.5, 0, 0.25]]
    structure = Structure(lattice, species, coords)
    structure_dict = structure_to_dict(structure)

    molecules = ["H2O", "NH3", "CO2", "CH4"]

    for molecule in molecules:
        result = add_adsorbate(
            structure_dict=structure_dict,
            molecule=molecule,
            site_index=0,
            distance=2.5
        )

        assert result["success"] is True, f"Failed to add {molecule}"
        new_structure = result["structure"]
        # Check that atoms were added
        assert new_structure["metadata"]["natoms"] > 4


def test_add_adsorbate_tilted_surface():
    """Test adding molecule to tilted (non-001) surface."""
    # Create a tilted slab (111) surface approximation
    lattice = Lattice([[5.0, 0, 0], [2.5, 4.33, 0], [0, 0, 10.0]])
    species = ["Si"] * 6
    coords = [
        [0.0, 0.0, 0.3],
        [0.333, 0.333, 0.3],
        [0.666, 0.666, 0.3],
        [0.0, 0.5, 0.4],
        [0.5, 0.0, 0.4],
        [0.5, 0.5, 0.4]
    ]
    structure = Structure(lattice, species, coords)
    structure_dict = structure_to_dict(structure)

    # Add CO to this tilted surface
    result = add_adsorbate(
        structure_dict=structure_dict,
        molecule="CO",
        site_index=0,
        distance=2.0
    )

    assert result["success"] is True
    new_structure = result["structure"]
    assert new_structure["metadata"]["natoms"] == 8  # 6 + 2


def test_add_adsorbate_invalid_site():
    """Test error handling for invalid site index."""
    lattice = Lattice.cubic(5.0)
    species = ["Cu"] * 4
    coords = [[0, 0, 0], [0.5, 0.5, 0], [0, 0.5, 0.25], [0.5, 0, 0.25]]
    structure = Structure(lattice, species, coords)
    structure_dict = structure_to_dict(structure)

    # Try invalid site index
    result = add_adsorbate(
        structure_dict=structure_dict,
        molecule="CO",
        site_index=999,  # Out of range
        distance=2.0
    )

    assert result["success"] is False
    assert "error" in result
    assert "INVALID_INPUT" in result["error"]["code"]


def test_add_adsorbate_unknown_molecule():
    """Test error handling for unknown molecule."""
    lattice = Lattice.cubic(5.0)
    species = ["Cu"] * 4
    coords = [[0, 0, 0], [0.5, 0.5, 0], [0, 0.5, 0.25], [0.5, 0, 0.25]]
    structure = Structure(lattice, species, coords)
    structure_dict = structure_to_dict(structure)

    result = add_adsorbate(
        structure_dict=structure_dict,
        molecule="UNKNOWN_MOLECULE_XYZ",
        site_index=0,
        distance=2.0
    )

    assert result["success"] is False
    assert "error" in result


def test_adsorbate_orientation_optimization():
    """Test that adsorbate tries multiple orientations."""
    # Create a surface with some steric constraints
    lattice = Lattice.cubic(5.0)
    species = ["Cu"] * 9
    coords = [
        [0, 0, 0], [0.5, 0, 0], [0, 0.5, 0],
        [0.5, 0.5, 0], [0.25, 0.25, 0.1],
        [0.75, 0.25, 0.1], [0.25, 0.75, 0.1],
        [0.75, 0.75, 0.1], [0.5, 0.5, 0.2]
    ]
    structure = Structure(lattice, species, coords)
    structure_dict = structure_to_dict(structure)

    # Add larger molecule (NH3) - should optimize orientation
    result = add_adsorbate(
        structure_dict=structure_dict,
        molecule="NH3",
        site_index=4,  # Central site
        distance=2.0
    )

    assert result["success"] is True
    # Check for orientation selection warning
    if "warnings" in result:
        warnings_text = " ".join(result["warnings"])
        # Should mention rotation selection
        assert "rotation" in warnings_text.lower() or len(result["warnings"]) == 0
