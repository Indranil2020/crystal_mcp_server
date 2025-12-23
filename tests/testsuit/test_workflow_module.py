"""
Tests for Workflow Module: Parametric Scans and Iterative Editing

This test file covers:
- Z-scan generation for molecule-surface systems
- Distance and rotation scans
- Parameter sweeps
- Iterative structure editing
- Trajectory/animation generation
"""

import pytest
import numpy as np
import sys
import os

# Add source path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../../src/python'))


# =============================================================================
# Fixtures
# =============================================================================

@pytest.fixture
def si_structure():
    """Create a simple Si structure for testing."""
    return {
        "species": ["Si", "Si"],
        "positions": [[0.0, 0.0, 0.0], [1.36, 1.36, 1.36]],
        "cell": [[5.43, 0, 0], [0, 5.43, 0], [0, 0, 5.43]],
        "pbc": [True, True, True],
        "n_atoms": 2,
        "formula": "Si2"
    }


@pytest.fixture
def nacl_slab():
    """Create a simple NaCl slab structure for testing."""
    a = 5.64
    positions = []
    species = []

    # Create a 2x2x3 slab
    for i in range(2):
        for j in range(2):
            for k in range(3):
                x = i * a/2
                y = j * a/2
                z = k * a/2
                positions.append([x, y, z])
                species.append("Na" if (i + j + k) % 2 == 0 else "Cl")

    return {
        "species": species,
        "positions": positions,
        "cell": [[a, 0, 0], [0, a, 0], [0, 0, 20.0]],  # 20 Å vacuum
        "pbc": [True, True, False],
        "n_atoms": len(species),
        "formula": "Na12Cl12"
    }


@pytest.fixture
def co_molecule():
    """Create a CO molecule structure."""
    return {
        "species": ["C", "O"],
        "positions": [[0.0, 0.0, 0.0], [0.0, 0.0, 1.13]],
        "cell": [[10, 0, 0], [0, 10, 0], [0, 0, 10]],
        "pbc": [False, False, False],
        "n_atoms": 2,
        "formula": "CO"
    }


# =============================================================================
# Parametric Scans Tests
# =============================================================================

class TestZScan:
    """Tests for z-direction scanning."""

    def test_z_scan_basic(self, nacl_slab, co_molecule):
        """Test basic z-scan generation."""
        from generators.workflow.parametric_scans import generate_z_scan

        result = generate_z_scan(
            surface_structure=nacl_slab,
            molecule=co_molecule,
            z_start=2.0,
            z_end=5.0,
            z_step=1.0
        )

        assert result["success"] is True
        assert "structures" in result
        assert len(result["structures"]) == 4  # 2.0, 3.0, 4.0, 5.0
        assert "trajectory_xyz" in result
        assert result["n_frames"] == 4

    def test_z_scan_with_material_name(self):
        """Test z-scan with material name instead of structure."""
        pytest.importorskip("ase")
        from generators.workflow.parametric_scans import generate_z_scan

        # This uses internal surface generation
        result = generate_z_scan(
            surface_material="NaCl",
            surface_miller=(1, 0, 0),
            molecule="CO",
            z_start=3.0,
            z_end=6.0,
            z_step=1.0,
            surface_layers=2,
            supercell=(2, 2, 1)
        )

        # May fail if ASE has issues, but should at least not raise
        assert "success" in result
        if result["success"]:
            assert len(result["structures"]) > 0

    def test_z_scan_z_values(self, nacl_slab, co_molecule):
        """Test that z_values are correctly generated."""
        from generators.workflow.parametric_scans import generate_z_scan

        result = generate_z_scan(
            surface_structure=nacl_slab,
            molecule=co_molecule,
            z_start=3.0,
            z_end=5.0,
            z_step=0.5
        )

        assert result["success"] is True
        expected_z = [3.0, 3.5, 4.0, 4.5, 5.0]
        assert len(result["z_values"]) == len(expected_z)
        for z_got, z_exp in zip(result["z_values"], expected_z):
            assert abs(z_got - z_exp) < 0.01

    def test_z_scan_energy_template(self, nacl_slab, co_molecule):
        """Test that energy calculation template is included."""
        from generators.workflow.parametric_scans import generate_z_scan

        result = generate_z_scan(
            surface_structure=nacl_slab,
            molecule=co_molecule,
            z_start=3.0,
            z_end=4.0,
            z_step=0.5,
            include_energy_template=True
        )

        assert result["success"] is True
        assert "energy_template" in result
        assert "calculation_steps" in result["energy_template"]


class TestDistanceScan:
    """Tests for distance scanning."""

    def test_distance_scan_basic(self, si_structure):
        """Test basic distance scan between two structures."""
        from generators.workflow.parametric_scans import generate_distance_scan

        result = generate_distance_scan(
            structure1=si_structure,
            structure2=si_structure,
            distance_start=3.0,
            distance_end=6.0,
            distance_step=1.0
        )

        assert result["success"] is True
        assert len(result["structures"]) == 4
        assert len(result["distances"]) == 4
        assert result["n_frames"] == 4

    def test_distance_scan_approach_axis(self, si_structure):
        """Test distance scan with different approach axes."""
        from generators.workflow.parametric_scans import generate_distance_scan

        for axis in ["x", "y", "z"]:
            result = generate_distance_scan(
                structure1=si_structure,
                structure2=si_structure,
                distance_start=4.0,
                distance_end=5.0,
                distance_step=0.5,
                approach_axis=axis
            )
            assert result["success"] is True
            assert result["approach_axis"] == axis


class TestRotationScan:
    """Tests for rotation scanning."""

    def test_rotation_scan_basic(self, si_structure):
        """Test basic rotation scan."""
        pytest.importorskip("scipy")
        from generators.workflow.parametric_scans import generate_rotation_scan

        result = generate_rotation_scan(
            structure=si_structure,
            rotation_axis="z",
            angle_start=0.0,
            angle_end=90.0,
            angle_step=30.0
        )

        assert result["success"] is True
        assert len(result["structures"]) == 4  # 0, 30, 60, 90
        assert len(result["angles"]) == 4

    def test_rotation_scan_full_rotation(self, si_structure):
        """Test full 360 degree rotation."""
        pytest.importorskip("scipy")
        from generators.workflow.parametric_scans import generate_rotation_scan

        result = generate_rotation_scan(
            structure=si_structure,
            rotation_axis="z",
            angle_start=0.0,
            angle_end=360.0,
            angle_step=90.0
        )

        assert result["success"] is True
        assert len(result["structures"]) == 5  # 0, 90, 180, 270, 360


class TestParameterSweep:
    """Tests for parameter sweep functionality."""

    def test_parameter_sweep_info(self):
        """Test that parameter sweep returns proper info."""
        from generators.workflow.parametric_scans import get_available_scans

        scans = get_available_scans()

        assert "z_scan" in scans
        assert "distance_scan" in scans
        assert "rotation_scan" in scans
        assert "parameter_sweep" in scans


class TestTrajectoryAnimation:
    """Tests for trajectory/animation generation."""

    def test_xyz_trajectory(self, si_structure):
        """Test XYZ trajectory generation."""
        from generators.workflow.parametric_scans import generate_trajectory_animation

        structures = [si_structure, si_structure, si_structure]
        result = generate_trajectory_animation(
            structures=structures,
            output_format="xyz"
        )

        assert result["success"] is True
        assert result["n_frames"] == 3
        assert "content" in result
        assert "Si" in result["content"]

    def test_xsf_animation(self, si_structure):
        """Test XSF animation generation."""
        from generators.workflow.parametric_scans import generate_trajectory_animation

        structures = [si_structure, si_structure]
        result = generate_trajectory_animation(
            structures=structures,
            output_format="xsf"
        )

        assert result["success"] is True
        assert "ANIMSTEPS" in result["content"]

    def test_extxyz_trajectory(self, si_structure):
        """Test extended XYZ trajectory generation."""
        from generators.workflow.parametric_scans import generate_trajectory_animation

        structures = [si_structure]
        result = generate_trajectory_animation(
            structures=structures,
            output_format="extxyz"
        )

        assert result["success"] is True
        assert "Lattice=" in result["content"]

    def test_html_animation(self, si_structure):
        """Test HTML animation viewer generation."""
        from generators.workflow.parametric_scans import generate_trajectory_animation

        structures = [si_structure, si_structure]
        result = generate_trajectory_animation(
            structures=structures,
            output_format="html"
        )

        assert result["success"] is True
        assert "<!DOCTYPE html>" in result["content"]
        assert "three.js" in result["content"].lower() or "Three" in result["content"]


# =============================================================================
# Iterative Editing Tests
# =============================================================================

class TestStructureEditor:
    """Tests for the StructureEditor class."""

    def test_editor_initialization(self, si_structure):
        """Test editor initialization."""
        from generators.workflow.iterative_editing import StructureEditor

        editor = StructureEditor(si_structure)
        struct = editor.get_structure()

        assert struct["n_atoms"] == 2
        assert struct["species"] == ["Si", "Si"]

    def test_add_atom(self, si_structure):
        """Test adding an atom."""
        from generators.workflow.iterative_editing import StructureEditor

        editor = StructureEditor(si_structure)
        editor.add_atom("C", [2.0, 2.0, 2.0])

        struct = editor.get_structure()
        assert len(struct["species"]) == 3
        assert "C" in struct["species"]

    def test_remove_atom(self, si_structure):
        """Test removing an atom."""
        from generators.workflow.iterative_editing import StructureEditor

        editor = StructureEditor(si_structure)
        editor.remove_atom(0)

        struct = editor.get_structure()
        assert len(struct["species"]) == 1

    def test_substitute_atom(self, si_structure):
        """Test substituting an atom."""
        from generators.workflow.iterative_editing import StructureEditor

        editor = StructureEditor(si_structure)
        editor.substitute_atom(0, "Ge")

        struct = editor.get_structure()
        assert "Ge" in struct["species"]
        assert struct["species"].count("Si") == 1

    def test_translate_atoms(self, si_structure):
        """Test translating atoms."""
        from generators.workflow.iterative_editing import StructureEditor

        editor = StructureEditor(si_structure)
        original = editor.get_structure()
        original_pos = np.array(original["positions"])

        editor.translate_atoms([1.0, 0.0, 0.0])

        struct = editor.get_structure()
        new_pos = np.array(struct["positions"])

        # All atoms should be shifted by 1.0 in x
        assert np.allclose(new_pos[:, 0], original_pos[:, 0] + 1.0)

    def test_rotate_atoms(self, si_structure):
        """Test rotating atoms."""
        pytest.importorskip("scipy")
        from generators.workflow.iterative_editing import StructureEditor

        editor = StructureEditor(si_structure)
        editor.rotate_atoms(angle=90.0, axis="z")

        # Should complete without error
        struct = editor.get_structure()
        assert len(struct["species"]) == 2

    def test_make_supercell(self, si_structure):
        """Test creating a supercell."""
        from generators.workflow.iterative_editing import StructureEditor

        editor = StructureEditor(si_structure)
        editor.make_supercell([2, 2, 2])

        struct = editor.get_structure()
        assert struct["n_atoms"] == 2 * 8  # 2 atoms * 2^3

    def test_apply_strain(self, si_structure):
        """Test applying strain."""
        from generators.workflow.iterative_editing import StructureEditor

        editor = StructureEditor(si_structure)
        original_cell = np.array(si_structure["cell"])

        editor.apply_strain(0.02, "biaxial")

        struct = editor.get_structure()
        new_cell = np.array(struct["cell"])

        # x and y should be strained, z unchanged
        assert new_cell[0, 0] > original_cell[0, 0]
        assert new_cell[1, 1] > original_cell[1, 1]
        assert abs(new_cell[2, 2] - original_cell[2, 2]) < 0.01

    def test_undo_redo(self, si_structure):
        """Test undo/redo functionality."""
        from generators.workflow.iterative_editing import StructureEditor

        editor = StructureEditor(si_structure)

        # Add an atom
        editor.add_atom("C", [2.0, 2.0, 2.0])
        assert len(editor.get_structure()["species"]) == 3

        # Undo
        success = editor.undo()
        assert success is True
        assert len(editor.get_structure()["species"]) == 2

        # Redo
        success = editor.redo()
        assert success is True
        assert len(editor.get_structure()["species"]) == 3

    def test_chained_operations(self, si_structure):
        """Test chaining multiple operations."""
        from generators.workflow.iterative_editing import StructureEditor

        editor = StructureEditor(si_structure)

        # Chain operations
        result = (editor
                  .make_supercell(2)
                  .substitute_atom(0, "Ge")
                  .translate_atoms([0.5, 0.5, 0.5])
                  .get_structure())

        assert result["n_atoms"] == 2 * 8
        assert "Ge" in result["species"]

    def test_history_tracking(self, si_structure):
        """Test that edit history is tracked."""
        from generators.workflow.iterative_editing import StructureEditor

        editor = StructureEditor(si_structure, track_history=True)

        editor.add_atom("C", [2.0, 2.0, 2.0])
        editor.remove_atom(0)

        history = editor.get_history()
        assert len(history) >= 2


class TestEditStructureFunction:
    """Tests for the edit_structure function."""

    def test_edit_structure_sequence(self, si_structure):
        """Test applying a sequence of operations."""
        from generators.workflow.iterative_editing import edit_structure

        result = edit_structure(
            structure=si_structure,
            operations=[
                {"operation": "make_supercell", "scaling": 2},
                {"operation": "substitute_atom", "index": 0, "new_element": "Ge"},
            ]
        )

        assert result["success"] is True
        assert result["n_operations"] == 2
        assert result["structure"]["n_atoms"] == 16
        assert "Ge" in result["structure"]["species"]

    def test_edit_structure_invalid_operation(self, si_structure):
        """Test handling of invalid operation."""
        from generators.workflow.iterative_editing import edit_structure

        result = edit_structure(
            structure=si_structure,
            operations=[
                {"operation": "nonexistent_operation"}
            ]
        )

        assert result["success"] is False
        assert "error" in result


class TestCompareStructures:
    """Tests for structure comparison."""

    def test_compare_identical(self, si_structure):
        """Test comparing identical structures."""
        from generators.workflow.iterative_editing import compare_structures

        result = compare_structures(si_structure, si_structure)

        assert result["success"] is True
        assert result["compositions_match"] is True

    def test_compare_different_composition(self, si_structure, co_molecule):
        """Test comparing structures with different composition."""
        from generators.workflow.iterative_editing import compare_structures

        result = compare_structures(si_structure, co_molecule)

        assert result["success"] is True
        assert result["compositions_match"] is False


class TestAvailableOperations:
    """Tests for operation listing functions."""

    def test_get_available_edit_operations(self):
        """Test listing available edit operations."""
        from generators.workflow.iterative_editing import get_available_edit_operations

        ops = get_available_edit_operations()

        assert "atom_operations" in ops
        assert "defect_operations" in ops
        assert "cell_operations" in ops
        assert "add_atom" in ops["atom_operations"]


# =============================================================================
# Integration Tests
# =============================================================================

class TestWorkflowIntegration:
    """Integration tests for workflow functionality."""

    def test_z_scan_to_animation(self, nacl_slab, co_molecule):
        """Test generating z-scan and converting to animation."""
        from generators.workflow.parametric_scans import (
            generate_z_scan,
            generate_trajectory_animation
        )

        # Generate z-scan
        scan_result = generate_z_scan(
            surface_structure=nacl_slab,
            molecule=co_molecule,
            z_start=3.0,
            z_end=5.0,
            z_step=1.0
        )

        assert scan_result["success"] is True

        # Convert to HTML animation
        anim_result = generate_trajectory_animation(
            structures=scan_result["structures"],
            output_format="html"
        )

        assert anim_result["success"] is True
        assert "<!DOCTYPE html>" in anim_result["content"]

    def test_edit_then_scan(self, si_structure):
        """Test editing structure then performing a scan."""
        from generators.workflow.iterative_editing import StructureEditor
        from generators.workflow.parametric_scans import generate_rotation_scan

        pytest.importorskip("scipy")

        # Edit the structure
        editor = StructureEditor(si_structure)
        editor.make_supercell(2)
        edited = editor.get_structure()

        # Perform rotation scan
        result = generate_rotation_scan(
            structure=edited,
            angle_start=0,
            angle_end=90,
            angle_step=30
        )

        assert result["success"] is True
        assert result["n_frames"] == 4


# =============================================================================
# Registry Integration Tests
# =============================================================================

class TestRegistryIntegration:
    """Test that workflow operations are properly registered."""

    def test_workflow_category_exists(self):
        """Test that workflow category is in registry."""
        from generators import GENERATOR_REGISTRY

        assert "workflow" in GENERATOR_REGISTRY

    def test_workflow_operations_registered(self):
        """Test that workflow operations are registered."""
        from generators import GENERATOR_REGISTRY

        workflow_ops = GENERATOR_REGISTRY["workflow"]["operations"]

        assert "generate_z_scan" in workflow_ops
        assert "generate_distance_scan" in workflow_ops
        assert "generate_rotation_scan" in workflow_ops
        assert "edit_structure" in workflow_ops
        assert "compare_structures" in workflow_ops

    def test_get_operation_info(self):
        """Test getting operation info for workflow operations."""
        from generators import get_operation_info

        info = get_operation_info("generate_z_scan")

        assert info["success"] is True
        assert info["info"]["category"] == "workflow"
        assert "surface_material" in info["info"]["params"]


# =============================================================================
# Edge Cases and Error Handling
# =============================================================================

# =============================================================================
# New Comprehensive Scan Tests
# =============================================================================

class TestDirectionalScan:
    """Tests for arbitrary direction scanning."""

    def test_directional_scan_basic(self, si_structure):
        """Test basic directional scan."""
        from generators.workflow.parametric_scans import generate_directional_scan

        result = generate_directional_scan(
            structure=si_structure,
            direction=[1, 0, 0],  # X direction
            distance_start=0.0,
            distance_end=2.0,
            distance_step=0.5
        )

        assert result["success"] is True
        assert len(result["structures"]) == 5  # 0, 0.5, 1.0, 1.5, 2.0
        assert "trajectory_xyz" in result

    def test_directional_scan_arbitrary_direction(self, si_structure):
        """Test directional scan with arbitrary direction vector."""
        from generators.workflow.parametric_scans import generate_directional_scan

        result = generate_directional_scan(
            structure=si_structure,
            direction=[1, 1, 1],  # Diagonal
            distance_start=0.0,
            distance_end=3.0,
            distance_step=1.0
        )

        assert result["success"] is True
        # Direction should be normalized
        import numpy as np
        norm = np.linalg.norm(result["direction"])
        assert abs(norm - 1.0) < 0.001


class TestPlanarScan:
    """Tests for 2D planar scanning."""

    def test_planar_scan_xy(self, si_structure):
        """Test planar scan in XY plane."""
        from generators.workflow.parametric_scans import generate_planar_scan

        result = generate_planar_scan(
            structure=si_structure,
            plane_normal=[0, 0, 1],  # XY plane
            x_range=(-1.0, 1.0, 1.0),  # 3 points
            y_range=(-1.0, 1.0, 1.0)   # 3 points
        )

        assert result["success"] is True
        assert result["grid_shape"] == (3, 3)
        assert result["n_structures"] == 9


class TestRadialScan:
    """Tests for radial scanning."""

    def test_radial_scan_basic(self, si_structure):
        """Test basic radial scan."""
        pytest.importorskip("scipy")
        from generators.workflow.parametric_scans import generate_radial_scan

        result = generate_radial_scan(
            structure=si_structure,
            center=[5.0, 5.0, 5.0],
            radius_start=1.0,
            radius_end=3.0,
            radius_step=1.0,
            n_angles=4
        )

        assert result["success"] is True
        # 3 radii * 4 angles = 12 structures
        assert result["n_structures"] == 12


class TestLatticeScan:
    """Tests for lattice parameter scanning."""

    def test_volume_scan(self, si_structure):
        """Test volume scaling scan."""
        from generators.workflow.parametric_scans import generate_lattice_scan

        result = generate_lattice_scan(
            structure=si_structure,
            volume_range=(0.9, 1.1, 0.1)  # 3 points
        )

        assert result["success"] is True
        assert result["n_structures"] == 3

    def test_lattice_a_scan(self, si_structure):
        """Test lattice parameter a scan."""
        from generators.workflow.parametric_scans import generate_lattice_scan

        result = generate_lattice_scan(
            structure=si_structure,
            a_range=(5.0, 6.0, 0.5)  # 3 points
        )

        assert result["success"] is True
        assert result["n_structures"] == 3


class TestStrainTensorScan:
    """Tests for strain tensor scanning."""

    def test_strain_xx_scan(self, si_structure):
        """Test uniaxial strain scan."""
        from generators.workflow.parametric_scans import generate_strain_tensor_scan

        result = generate_strain_tensor_scan(
            structure=si_structure,
            strain_component="xx",
            strain_range=(-0.02, 0.02, 0.02)  # 3 points
        )

        assert result["success"] is True
        assert result["n_structures"] == 3
        assert result["strain_component"] == "xx"

    def test_biaxial_strain_scan(self, si_structure):
        """Test biaxial strain scan."""
        from generators.workflow.parametric_scans import generate_strain_tensor_scan

        result = generate_strain_tensor_scan(
            structure=si_structure,
            strain_component="biaxial",
            strain_range=(-0.05, 0.05, 0.05)
        )

        assert result["success"] is True
        assert result["strain_component"] == "biaxial"


class TestCoverageScan:
    """Tests for coverage scanning."""

    def test_coverage_scan_basic(self, nacl_slab, co_molecule):
        """Test basic coverage scan."""
        from generators.workflow.parametric_scans import generate_coverage_scan

        result = generate_coverage_scan(
            surface_structure=nacl_slab,
            adsorbate=co_molecule,
            coverage_range=(0.25, 0.75, 0.25)
        )

        assert result["success"] is True
        assert result["n_structures"] == 3


class TestNEBPath:
    """Tests for NEB path generation."""

    def test_neb_path_basic(self, si_structure):
        """Test basic NEB path generation."""
        from generators.workflow.parametric_scans import generate_neb_path
        import copy

        # Create final structure with shifted positions
        final_structure = copy.deepcopy(si_structure)
        final_structure["positions"][0][0] += 1.0

        result = generate_neb_path(
            initial_structure=si_structure,
            final_structure=final_structure,
            n_images=5
        )

        assert result["success"] is True
        # 5 images + 2 endpoints = 7 structures
        assert result["n_structures"] == 7
        assert "neb_info" in result

    def test_neb_path_climbing_image(self, si_structure):
        """Test NEB path with climbing image."""
        from generators.workflow.parametric_scans import generate_neb_path
        import copy

        final_structure = copy.deepcopy(si_structure)
        final_structure["positions"][0][0] += 1.0

        result = generate_neb_path(
            initial_structure=si_structure,
            final_structure=final_structure,
            n_images=5,
            climbing_image=True
        )

        assert result["success"] is True
        # Check that one image is marked as climbing
        climbing_count = sum(1 for s in result["structures"] if s.get("is_climbing_image"))
        assert climbing_count == 1


class TestBondScan:
    """Tests for bond distance scanning."""

    def test_bond_scan_basic(self, co_molecule):
        """Test basic bond scan."""
        from generators.workflow.parametric_scans import generate_bond_scan

        result = generate_bond_scan(
            structure=co_molecule,
            atom_indices=(0, 1),
            bond_range=(1.0, 2.0, 0.5)
        )

        assert result["success"] is True
        assert result["n_structures"] == 3
        assert result["atom_indices"] == [0, 1]


class TestDihedralScan:
    """Tests for dihedral angle scanning."""

    def test_dihedral_scan_basic(self):
        """Test basic dihedral scan."""
        pytest.importorskip("scipy")
        from generators.workflow.parametric_scans import generate_dihedral_scan

        # Create a simple 4-atom structure
        structure = {
            "species": ["C", "C", "C", "C"],
            "positions": [
                [0.0, 0.0, 0.0],
                [1.5, 0.0, 0.0],
                [1.5, 1.5, 0.0],
                [3.0, 1.5, 0.0]
            ],
            "cell": [[10, 0, 0], [0, 10, 0], [0, 0, 10]],
            "n_atoms": 4
        }

        result = generate_dihedral_scan(
            structure=structure,
            atom_indices=(0, 1, 2, 3),
            angle_range=(0.0, 180.0, 60.0)
        )

        assert result["success"] is True
        assert result["n_structures"] == 4  # 0, 60, 120, 180


class TestNewScansInRegistry:
    """Tests for new scan operations in registry."""

    def test_new_scans_registered(self):
        """Test that new scan operations are in registry."""
        from generators import GENERATOR_REGISTRY

        workflow_ops = GENERATOR_REGISTRY["workflow"]["operations"]

        new_scans = [
            "generate_directional_scan",
            "generate_planar_scan",
            "generate_radial_scan",
            "generate_lattice_scan",
            "generate_strain_tensor_scan",
            "generate_coverage_scan",
            "generate_neb_path",
            "generate_bond_scan",
            "generate_dihedral_scan"
        ]

        for scan in new_scans:
            assert scan in workflow_ops, f"{scan} not in registry"

    def test_get_available_scans_complete(self):
        """Test that get_available_scans returns all scan types."""
        from generators.workflow.parametric_scans import get_available_scans

        scans = get_available_scans()

        expected_scans = [
            "z_scan", "distance_scan", "rotation_scan",
            "directional_scan", "planar_scan", "radial_scan",
            "lattice_scan", "strain_tensor_scan", "coverage_scan",
            "neb_path", "bond_scan", "dihedral_scan",
            "parameter_sweep", "multi_parameter_scan", "trajectory_animation"
        ]

        for scan in expected_scans:
            assert scan in scans, f"{scan} missing from available scans"


class TestEdgeCases:
    """Test edge cases and error handling."""

    def test_empty_structure(self):
        """Test handling of empty structure."""
        from generators.workflow.iterative_editing import StructureEditor

        empty = {"species": [], "positions": [], "cell": [[1, 0, 0], [0, 1, 0], [0, 0, 1]]}
        editor = StructureEditor(empty)

        # Adding to empty should work
        editor.add_atom("H", [0.5, 0.5, 0.5])
        assert len(editor.get_structure()["species"]) == 1

    def test_no_structure_set(self):
        """Test operations on uninitialized editor."""
        from generators.workflow.iterative_editing import StructureEditor

        editor = StructureEditor()

        with pytest.raises(ValueError):
            editor.add_atom("H", [0, 0, 0])

    def test_invalid_atom_index(self, si_structure):
        """Test removing atom with invalid index."""
        from generators.workflow.iterative_editing import StructureEditor

        editor = StructureEditor(si_structure)

        with pytest.raises(IndexError):
            editor.remove_atom(999)

    def test_single_frame_trajectory(self, si_structure):
        """Test trajectory with single frame."""
        from generators.workflow.parametric_scans import generate_trajectory_animation

        result = generate_trajectory_animation(
            structures=[si_structure],
            output_format="xyz"
        )

        assert result["success"] is True
        assert result["n_frames"] == 1


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
