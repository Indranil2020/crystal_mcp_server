"""
test_phase3_implementations.py - Tests for Phase 3 implementations

Tests for:
1. Visualization exports (XSF, CUBE, VESTA, Three.js, HTML)
2. Cell settings (Primitive, Conventional, Niggli, Standardized)
3. Advanced magnetic orderings (multi-q, DMI cycloid, conical, skyrmion)
4. Defects (trivacancy, gradient, vacancy cluster)
5. Ferroelectrics (Dion-Jacobson, Aurivillius)
"""

import pytest
import sys
import os
import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../../src/python'))


# =============================================================================
# Test Visualization Exports
# =============================================================================
class TestVisualizationExports:
    """Test visualization export formats."""

    def test_module_importable(self):
        """Test that visualization module is importable."""
        from generators.output_formats.visualization import (
            export_xsf, export_cube, export_vesta, export_threejs_json, export_html_viewer
        )
        assert callable(export_xsf)
        assert callable(export_cube)
        assert callable(export_vesta)
        assert callable(export_threejs_json)
        assert callable(export_html_viewer)

    def test_export_xsf(self):
        """Test XSF export."""
        from generators.output_formats.visualization import export_xsf
        from generators.bulk.spacegroups import generate_prototype

        # Generate a simple structure
        si = generate_prototype("diamond", {"X": "Si"})
        structure = si.get("structure", {})

        result = export_xsf(structure)
        assert result["success"] is True
        assert "content" in result
        assert "CRYSTAL" in result["content"]
        assert "PRIMVEC" in result["content"]
        assert "PRIMCOORD" in result["content"]

    def test_export_cube(self):
        """Test Gaussian CUBE export."""
        from generators.output_formats.visualization import export_cube
        from generators.bulk.spacegroups import generate_prototype

        si = generate_prototype("diamond", {"X": "Si"})
        structure = si.get("structure", {})

        result = export_cube(structure)
        assert result["success"] is True
        assert "content" in result
        assert "cube" in result["format"]

    def test_export_vesta(self):
        """Test VESTA CIF export."""
        from generators.output_formats.visualization import export_vesta
        from generators.bulk.spacegroups import generate_prototype

        si = generate_prototype("diamond", {"X": "Si"})
        structure = si.get("structure", {})

        result = export_vesta(structure)
        assert result["success"] is True
        assert "content" in result
        assert "_cell_length_a" in result["content"]

    def test_export_threejs_json(self):
        """Test Three.js JSON export."""
        from generators.output_formats.visualization import export_threejs_json
        from generators.bulk.spacegroups import generate_prototype

        si = generate_prototype("diamond", {"X": "Si"})
        structure = si.get("structure", {})

        result = export_threejs_json(structure)
        assert result["success"] is True
        # Check for expected structure - might be "json_data" or similar
        assert "json_data" in result or "data" in result

    def test_export_html_viewer(self):
        """Test HTML viewer export."""
        from generators.output_formats.visualization import export_html_viewer
        from generators.bulk.spacegroups import generate_prototype

        si = generate_prototype("diamond", {"X": "Si"})
        structure = si.get("structure", {})

        result = export_html_viewer(structure, title="Silicon")
        assert result["success"] is True
        assert "content" in result
        assert "<html>" in result["content"].lower() or "<!doctype" in result["content"].lower()


# =============================================================================
# Test Cell Settings
# =============================================================================
class TestCellSettings:
    """Test cell setting transformations."""

    def test_module_importable(self):
        """Test that cell_settings module is importable."""
        from generators.bulk.cell_settings import (
            get_primitive_cell, get_conventional_cell, get_niggli_cell,
            get_standardized_cell, transform_cell, compare_cell_settings
        )
        assert callable(get_primitive_cell)
        assert callable(get_conventional_cell)
        assert callable(get_niggli_cell)

    def test_get_primitive_cell(self):
        """Test primitive cell transformation."""
        from generators.bulk.cell_settings import get_primitive_cell
        from generators.bulk.spacegroups import generate_prototype

        # FCC structure has primitive cell with 1 atom
        fcc = generate_prototype("fcc", {"X": "Cu"})
        structure = fcc.get("structure", {})

        result = get_primitive_cell(structure)
        assert result["success"] is True
        # Check for structure output
        assert "structure" in result

    def test_get_conventional_cell(self):
        """Test conventional cell transformation."""
        from generators.bulk.cell_settings import get_conventional_cell
        from generators.bulk.spacegroups import generate_prototype

        bcc = generate_prototype("bcc", {"X": "Fe"})
        structure = bcc.get("structure", {})

        result = get_conventional_cell(structure)
        assert result["success"] is True
        assert "structure" in result

    def test_get_niggli_cell(self):
        """Test Niggli cell transformation."""
        from generators.bulk.cell_settings import get_niggli_cell
        from generators.bulk.spacegroups import generate_prototype

        si = generate_prototype("diamond", {"X": "Si"})
        structure = si.get("structure", {})

        result = get_niggli_cell(structure)
        assert result["success"] is True
        assert "niggli_parameters" in result

    def test_transform_cell(self):
        """Test transform_cell convenience function."""
        from generators.bulk.cell_settings import transform_cell
        from generators.bulk.spacegroups import generate_prototype

        si = generate_prototype("diamond", {"X": "Si"})
        structure = si.get("structure", {})

        for setting in ["primitive", "conventional", "niggli"]:
            result = transform_cell(structure, setting=setting)
            assert result["success"] is True, f"Failed for {setting}"

    def test_compare_cell_settings(self):
        """Test comparison of all cell settings."""
        from generators.bulk.cell_settings import compare_cell_settings
        from generators.bulk.spacegroups import generate_prototype

        si = generate_prototype("diamond", {"X": "Si"})
        structure = si.get("structure", {})

        result = compare_cell_settings(structure)
        assert result["success"] is True
        # Check for comparison data - might have different structure
        assert "comparison" in result or "settings" in result or "primitive" in result


# =============================================================================
# Test Advanced Magnetic Orderings
# =============================================================================
class TestAdvancedMagneticOrderings:
    """Test multi-q, DMI cycloid, conical, and skyrmion structures."""

    def test_module_importable(self):
        """Test that advanced_orderings module is importable."""
        from generators.magnetic.advanced_orderings import (
            generate_multi_q_structure, generate_dmi_cycloid,
            generate_conical_phase, generate_skyrmion_lattice
        )
        assert callable(generate_multi_q_structure)
        assert callable(generate_dmi_cycloid)
        assert callable(generate_conical_phase)
        assert callable(generate_skyrmion_lattice)

    def test_single_q_structure(self):
        """Test single-q (simple spiral) structure."""
        from generators.magnetic.advanced_orderings import generate_multi_q_structure

        result = generate_multi_q_structure("Fe", q_type="single_q", supercell=[2, 2, 2])
        assert result["success"] is True
        assert result["n_q_vectors"] == 1
        assert len(result["magnetic_moments"]) == result["n_atoms"]

    def test_2q_structure(self):
        """Test 2-q magnetic structure."""
        from generators.magnetic.advanced_orderings import generate_multi_q_structure

        result = generate_multi_q_structure("Fe", q_type="2q", supercell=[2, 2, 2])
        assert result["success"] is True
        assert result["n_q_vectors"] == 2

    def test_3q_structure(self):
        """Test 3-q (triple-q/tetrahedral) structure."""
        from generators.magnetic.advanced_orderings import generate_multi_q_structure

        result = generate_multi_q_structure("Fe", q_type="3q", supercell=[2, 2, 2])
        assert result["success"] is True
        assert result["n_q_vectors"] == 3
        assert "tetrahedral" in result.get("description", "").lower() or "triple" in result.get("description", "").lower()

    def test_4q_structure(self):
        """Test 4-q magnetic structure."""
        from generators.magnetic.advanced_orderings import generate_multi_q_structure

        result = generate_multi_q_structure("Fe", q_type="4q", supercell=[2, 2, 2])
        assert result["success"] is True
        assert result["n_q_vectors"] == 4

    def test_dmi_cycloid_bifeo3(self):
        """Test DMI cycloid for BiFeO3."""
        from generators.magnetic.advanced_orderings import generate_dmi_cycloid

        result = generate_dmi_cycloid("BiFeO3")
        assert result["success"] is True
        assert result["cycloid_type"] == "spin_cycloid"
        assert result["cycloid_period_nm"] == 62.0  # Known BiFeO3 period
        assert result["is_multiferroic"] is True

    def test_dmi_cycloid_mnsi(self):
        """Test DMI cycloid for MnSi."""
        from generators.magnetic.advanced_orderings import generate_dmi_cycloid

        result = generate_dmi_cycloid("MnSi")
        assert result["success"] is True
        assert result["can_host_skyrmions"] is True

    def test_conical_phase(self):
        """Test conical magnetic phase."""
        from generators.magnetic.advanced_orderings import generate_conical_phase

        result = generate_conical_phase("MnSi", cone_angle_deg=45.0)
        assert result["success"] is True
        assert result["phase"] == "conical"
        assert result["cone_angle_deg"] == 45.0

    def test_skyrmion_lattice_bloch(self):
        """Test Bloch skyrmion lattice."""
        from generators.magnetic.advanced_orderings import generate_skyrmion_lattice

        result = generate_skyrmion_lattice("MnSi", skyrmion_type="bloch", supercell_size=10)
        assert result["success"] is True
        assert result["skyrmion_type"] == "bloch"
        assert result["topological_charge_per_skyrmion"] == -1

    def test_skyrmion_lattice_neel(self):
        """Test Neel skyrmion lattice."""
        from generators.magnetic.advanced_orderings import generate_skyrmion_lattice

        result = generate_skyrmion_lattice("FeGe", skyrmion_type="neel", supercell_size=10)
        assert result["success"] is True
        assert result["skyrmion_type"] == "neel"


# =============================================================================
# Test Advanced Defects
# =============================================================================
class TestAdvancedDefects:
    """Test trivacancy, gradient, and vacancy cluster defects."""

    def test_module_importable(self):
        """Test that point_defects module has new functions."""
        from generators.defect.point_defects import (
            generate_trivacancy, generate_defect_gradient, generate_vacancy_cluster
        )
        assert callable(generate_trivacancy)
        assert callable(generate_defect_gradient)
        assert callable(generate_vacancy_cluster)

    def test_trivacancy_triangular(self):
        """Test triangular trivacancy."""
        from generators.defect.point_defects import generate_trivacancy
        from generators.bulk.spacegroups import generate_prototype

        # Create supercell
        si = generate_prototype("diamond", {"X": "Si"})
        from pymatgen.core import Structure, Lattice
        lattice = Lattice.cubic(5.43)
        structure = Structure(lattice, ["Si"]*8, [
            [0, 0, 0], [0.5, 0.5, 0], [0.5, 0, 0.5], [0, 0.5, 0.5],
            [0.25, 0.25, 0.25], [0.75, 0.75, 0.25], [0.75, 0.25, 0.75], [0.25, 0.75, 0.75]
        ])
        host = {
            "lattice": {"matrix": lattice.matrix.tolist()},
            "atoms": [{"element": "Si", "coords": list(s.frac_coords)} for s in structure]
        }

        result = generate_trivacancy(host, center_site=0, trivacancy_type="triangular")
        assert result["success"] is True
        assert result["n_vacancies"] == 3
        assert result["trivacancy_type"] == "triangular"

    def test_trivacancy_linear(self):
        """Test linear trivacancy."""
        from generators.defect.point_defects import generate_trivacancy
        from pymatgen.core import Structure, Lattice

        lattice = Lattice.cubic(5.43)
        structure = Structure(lattice, ["Si"]*8, [
            [0, 0, 0], [0.5, 0.5, 0], [0.5, 0, 0.5], [0, 0.5, 0.5],
            [0.25, 0.25, 0.25], [0.75, 0.75, 0.25], [0.75, 0.25, 0.75], [0.25, 0.75, 0.75]
        ])
        host = {
            "lattice": {"matrix": lattice.matrix.tolist()},
            "atoms": [{"element": "Si", "coords": list(s.frac_coords)} for s in structure]
        }

        result = generate_trivacancy(host, center_site=0, trivacancy_type="linear")
        assert result["success"] is True
        assert result["n_vacancies"] == 3

    def test_defect_gradient_linear(self):
        """Test linear defect gradient."""
        from generators.defect.point_defects import generate_defect_gradient
        from pymatgen.core import Structure, Lattice

        # Create GaAs supercell
        lattice = Lattice.cubic(5.65 * 3)
        coords = []
        species = []
        for i in range(3):
            for j in range(3):
                for k in range(3):
                    species.extend(["Ga", "As"])
                    coords.append([(i + 0)/3, (j + 0)/3, (k + 0)/3])
                    coords.append([(i + 0.25)/3, (j + 0.25)/3, (k + 0.25)/3])

        host = {
            "lattice": {"matrix": lattice.matrix.tolist()},
            "atoms": [{"element": s, "coords": c} for s, c in zip(species, coords)]
        }

        result = generate_defect_gradient(
            host, defect_element="N", gradient_direction=[0, 0, 1],
            min_concentration=0.0, max_concentration=0.5,
            gradient_profile="linear", host_element="As", seed=42
        )
        assert result["success"] is True
        assert result["gradient_profile"] == "linear"
        assert result["n_substituted"] > 0

    def test_defect_gradient_sigmoid(self):
        """Test sigmoid defect gradient."""
        from generators.defect.point_defects import generate_defect_gradient
        from pymatgen.core import Lattice

        lattice = Lattice.cubic(5.65 * 3)
        coords = []
        species = []
        for i in range(3):
            for j in range(3):
                for k in range(3):
                    species.extend(["Ga", "As"])
                    coords.append([(i + 0)/3, (j + 0)/3, (k + 0)/3])
                    coords.append([(i + 0.25)/3, (j + 0.25)/3, (k + 0.25)/3])

        host = {
            "lattice": {"matrix": lattice.matrix.tolist()},
            "atoms": [{"element": s, "coords": c} for s, c in zip(species, coords)]
        }

        result = generate_defect_gradient(
            host, defect_element="N", gradient_profile="sigmoid",
            max_concentration=0.3, seed=42
        )
        assert result["success"] is True
        assert result["gradient_profile"] == "sigmoid"

    def test_vacancy_cluster_compact(self):
        """Test compact vacancy cluster."""
        from generators.defect.point_defects import generate_vacancy_cluster
        from pymatgen.core import Lattice

        lattice = Lattice.cubic(5.43 * 2)
        coords = []
        species = []
        for i in range(2):
            for j in range(2):
                for k in range(2):
                    for base in [[0, 0, 0], [0.5, 0.5, 0], [0.5, 0, 0.5], [0, 0.5, 0.5]]:
                        species.append("Si")
                        coords.append([(i + base[0])/2, (j + base[1])/2, (k + base[2])/2])

        host = {
            "lattice": {"matrix": lattice.matrix.tolist()},
            "atoms": [{"element": s, "coords": c} for s, c in zip(species, coords)]
        }

        result = generate_vacancy_cluster(host, cluster_size=4, cluster_type="compact")
        assert result["success"] is True
        assert result["n_vacancies"] == 4
        assert result["cluster_type"] == "compact"


# =============================================================================
# Test Ferroelectric Phases
# =============================================================================
class TestFerroelectricPhases:
    """Test Dion-Jacobson and Aurivillius phases."""

    def test_module_importable(self):
        """Test that perovskites module has ferroelectric functions."""
        from generators.bulk.perovskites import (
            generate_dion_jacobson, generate_aurivillius, get_layered_perovskite_database
        )
        assert callable(generate_dion_jacobson)
        assert callable(generate_aurivillius)
        assert callable(get_layered_perovskite_database)

    def test_dion_jacobson_klanb2o7(self):
        """Test Dion-Jacobson KLaNb2O7."""
        from generators.bulk.perovskites import generate_dion_jacobson

        result = generate_dion_jacobson("KLaNb2O7")
        assert result["success"] is True
        assert result["family"] == "dion-jacobson"
        assert result["n_layers"] == 2
        assert result["alkali_element"] == "K"
        assert result["has_rock_salt_layer"] is False
        assert result["ion_exchange_capable"] is True

    def test_dion_jacobson_n3(self):
        """Test Dion-Jacobson with n=3 layers."""
        from generators.bulk.perovskites import generate_dion_jacobson

        result = generate_dion_jacobson("KCa2Nb3O10")
        assert result["success"] is True
        assert result["n_layers"] == 3

    def test_aurivillius_bi4ti3o12(self):
        """Test Aurivillius Bi4Ti3O12."""
        from generators.bulk.perovskites import generate_aurivillius

        result = generate_aurivillius("Bi4Ti3O12")
        assert result["success"] is True
        assert result["family"] == "aurivillius"
        assert result["n_layers"] == 3
        assert result["has_bi2o2_layer"] is True
        assert result["is_ferroelectric"] is True

    def test_aurivillius_srbi2ta2o9(self):
        """Test Aurivillius SrBi2Ta2O9 (fatigue-free ferroelectric)."""
        from generators.bulk.perovskites import generate_aurivillius

        result = generate_aurivillius("SrBi2Ta2O9")
        assert result["success"] is True
        assert result["n_layers"] == 2
        assert result["is_fatigue_free"] is True

    def test_aurivillius_n1(self):
        """Test Aurivillius with n=1 layer (Bi2WO6)."""
        from generators.bulk.perovskites import generate_aurivillius

        result = generate_aurivillius("Bi2WO6")
        assert result["success"] is True
        assert result["n_layers"] == 1

    def test_layered_perovskite_database(self):
        """Test layered perovskite database."""
        from generators.bulk.perovskites import get_layered_perovskite_database

        result = get_layered_perovskite_database()
        assert result["success"] is True
        assert "families" in result
        assert "ruddlesden_popper" in result["families"]
        assert "dion_jacobson" in result["families"]
        assert "aurivillius" in result["families"]


# =============================================================================
# Test Registry Integration
# =============================================================================
class TestPhase3RegistryIntegration:
    """Test that new operations are registered in the operation registry."""

    def _find_operation(self, op_name):
        """Helper to find operation in registry."""
        from generators import GENERATOR_REGISTRY
        for cat, cat_data in GENERATOR_REGISTRY.items():
            if isinstance(cat_data, dict) and "operations" in cat_data:
                if op_name in cat_data["operations"]:
                    return True
        return False

    def test_visualization_operations_registered(self):
        """Test visualization operations are in registry."""
        viz_ops = ["export_xsf", "export_cube", "export_vesta", "export_threejs_json", "export_html_viewer"]
        for op in viz_ops:
            assert self._find_operation(op), f"Operation {op} not found in registry"

    def test_cell_settings_operations_registered(self):
        """Test cell settings operations are in registry."""
        cell_ops = ["get_primitive_cell", "get_conventional_cell", "get_niggli_cell",
                    "get_standardized_cell", "transform_cell", "compare_cell_settings"]
        for op in cell_ops:
            assert self._find_operation(op), f"Operation {op} not found in registry"

    def test_magnetic_operations_registered(self):
        """Test advanced magnetic operations are in registry."""
        mag_ops = ["generate_multi_q_structure", "generate_dmi_cycloid",
                   "generate_conical_phase", "generate_skyrmion_lattice"]
        for op in mag_ops:
            assert self._find_operation(op), f"Operation {op} not found in registry"

    def test_defect_operations_registered(self):
        """Test advanced defect operations are in registry."""
        defect_ops = ["generate_trivacancy", "generate_defect_gradient", "generate_vacancy_cluster"]
        for op in defect_ops:
            assert self._find_operation(op), f"Operation {op} not found in registry"

    def test_ferroelectric_operations_registered(self):
        """Test ferroelectric operations are in registry."""
        ferro_ops = ["generate_dion_jacobson", "generate_aurivillius", "get_layered_perovskite_database"]
        for op in ferro_ops:
            assert self._find_operation(op), f"Operation {op} not found in registry"


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
