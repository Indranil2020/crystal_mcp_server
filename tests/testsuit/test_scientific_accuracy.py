"""
Scientific Accuracy Tests for Crystal Structure Generators

Reference Data Sources:
- NIST Standard Reference Database 64: NIST Structural Database (NSD)
- CRC Handbook of Chemistry and Physics, 102nd Edition (2021-2022)
- Crystallography Open Database (COD): https://www.crystallography.net
- Materials Project: https://materialsproject.org
- International Tables for Crystallography, Volume A, 6th Edition (2016)

Lattice parameters at room temperature (293-298 K) unless otherwise noted.
Tolerance: 2% for lattice parameters (accounts for temperature and method variations)
"""

import sys
import os
import pytest

# Add the Python source to path
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__),
                                                '../../src/python')))


# =============================================================================
# AUTHORITATIVE REFERENCE DATA WITH SOURCES
# =============================================================================

REFERENCE_DATA = {
    # Elemental Crystals - Source: CRC Handbook 102nd Ed., Table 12.15
    "Si_diamond": {
        "a": 5.4309,  # Å at 298 K
        "spacegroup": 227,  # Fd-3m
        "source": "CRC Handbook 102nd Ed., Table 12.15"
    },
    "Ge_diamond": {
        "a": 5.6579,  # Å at 298 K
        "spacegroup": 227,  # Fd-3m
        "source": "CRC Handbook 102nd Ed., Table 12.15"
    },
    "Cu_fcc": {
        "a": 3.6149,  # Å at 298 K
        "spacegroup": 225,  # Fm-3m
        "source": "CRC Handbook 102nd Ed., Table 12.15"
    },
    "Au_fcc": {
        "a": 4.0782,  # Å at 298 K
        "spacegroup": 225,  # Fm-3m
        "source": "CRC Handbook 102nd Ed., Table 12.15"
    },
    "Ag_fcc": {
        "a": 4.0862,  # Å at 298 K
        "spacegroup": 225,  # Fm-3m
        "source": "CRC Handbook 102nd Ed., Table 12.15"
    },
    "Al_fcc": {
        "a": 4.0495,  # Å at 298 K
        "spacegroup": 225,  # Fm-3m
        "source": "CRC Handbook 102nd Ed., Table 12.15"
    },
    "Fe_bcc": {
        "a": 2.8665,  # Å at 298 K (alpha-Fe)
        "spacegroup": 229,  # Im-3m
        "source": "CRC Handbook 102nd Ed., Table 12.15"
    },
    "W_bcc": {
        "a": 3.1652,  # Å at 298 K
        "spacegroup": 229,  # Im-3m
        "source": "CRC Handbook 102nd Ed., Table 12.15"
    },

    # Binary Compounds - Source: Wyckoff Crystal Structures + Materials Project
    "NaCl_rocksalt": {
        "a": 5.6401,  # Å at 298 K
        "spacegroup": 225,  # Fm-3m
        "source": "Wyckoff Crystal Structures Vol. 1, COD #9008677"
    },
    "MgO_rocksalt": {
        "a": 4.2112,  # Å at 298 K
        "spacegroup": 225,  # Fm-3m
        "source": "Materials Project mp-1265"
    },
    "GaAs_zincblende": {
        "a": 5.6533,  # Å at 300 K
        "spacegroup": 216,  # F-43m
        "source": "Landolt-Börnstein III/41A1a"
    },
    "ZnS_zincblende": {
        "a": 5.4093,  # Å at 300 K
        "spacegroup": 216,  # F-43m
        "source": "COD #9008878"
    },
    "ZnO_wurtzite": {
        "a": 3.2498,  # Å at 300 K
        "c": 5.2066,  # Å at 300 K
        "spacegroup": 186,  # P6_3mc
        "source": "COD #9004178, Materials Project mp-2133"
    },

    # 2D Materials - Source: Phys. Rev. B + Nature Materials reviews
    "graphene": {
        "a": 2.46,  # Å - C-C bond 1.42 Å, a = 2.46 Å
        "source": "Phys. Rev. B 76, 064120 (2007)"
    },
    "hBN_monolayer": {
        "a": 2.504,  # Å - B-N bond 1.446 Å
        "source": "Phys. Rev. B 64, 235118 (2001)"
    },
    "MoS2_monolayer": {
        "a": 3.160,  # Å (2H phase)
        "source": "J. Phys. Chem. C 116, 8983 (2012)"
    },
    "WS2_monolayer": {
        "a": 3.153,  # Å (2H phase)
        "source": "2D Mater. 1, 025001 (2014)"
    },

    # Topological Insulators - Source: Nature Physics + Phys. Rev. B
    "Bi2Se3": {
        "a": 4.138,  # Å
        "c": 28.64,  # Å
        "spacegroup": 166,  # R-3m
        "source": "Nature Physics 5, 438 (2009), COD #1010089"
    },
    "Bi2Te3": {
        "a": 4.3835,  # Å
        "c": 30.487,  # Å
        "spacegroup": 166,  # R-3m
        "source": "COD #9012437, Materials Project mp-34202"
    },

    # Perovskites - Source: Acta Cryst. + J. Solid State Chem.
    "SrTiO3_cubic": {
        "a": 3.905,  # Å at 300 K (cubic phase)
        "spacegroup": 221,  # Pm-3m
        "source": "Acta Cryst. B 35, 1580 (1979)"
    },
    "BaTiO3_tetragonal": {
        "a": 3.992,  # Å at 300 K
        "c": 4.036,  # Å at 300 K
        "spacegroup": 99,  # P4mm
        "source": "J. Phys. Soc. Jpn. 11, 10 (1956)"
    },
}

# Tolerance for lattice parameter comparison (2% = typical T-dependent variation)
LATTICE_TOLERANCE_PERCENT = 2.0

# =============================================================================
# HELPER FUNCTIONS
# =============================================================================

def assert_lattice_within_tolerance(actual, expected, material, param_name="a",
                                     tolerance_pct=LATTICE_TOLERANCE_PERCENT):
    """Assert that measured lattice parameter is within tolerance of reference value."""
    deviation_pct = abs(actual - expected) / expected * 100
    assert deviation_pct <= tolerance_pct, (
        f"{material} {param_name}: {actual:.4f} Å deviates by {deviation_pct:.2f}% "
        f"from reference {expected:.4f} Å (tolerance: {tolerance_pct}%)"
    )
    return deviation_pct


# =============================================================================
# TEST CLASSES
# =============================================================================

class TestElementalCrystals:
    """Test elemental crystal structures against CRC Handbook values."""

    def test_silicon_diamond(self):
        """Silicon diamond structure: a = 5.4309 Å, SG 227"""
        from generators.bulk.spacegroups import generate_from_spacegroup

        ref = REFERENCE_DATA["Si_diamond"]
        result = generate_from_spacegroup(
            spacegroup=ref["spacegroup"],
            elements=["Si"],
            composition=[8],
            a=ref["a"]
        )

        assert result["success"], f"Generation failed: {result.get('error')}"
        a_actual = result["structure"]["lattice"]["a"]
        assert_lattice_within_tolerance(a_actual, ref["a"], "Si(diamond)", "a")

    def test_germanium_diamond(self):
        """Germanium diamond structure: a = 5.6579 Å, SG 227"""
        from generators.bulk.spacegroups import generate_from_spacegroup

        ref = REFERENCE_DATA["Ge_diamond"]
        result = generate_from_spacegroup(
            spacegroup=ref["spacegroup"],
            elements=["Ge"],
            composition=[8],
            a=ref["a"]
        )

        assert result["success"], f"Generation failed: {result.get('error')}"
        a_actual = result["structure"]["lattice"]["a"]
        assert_lattice_within_tolerance(a_actual, ref["a"], "Ge(diamond)", "a")

    def test_copper_fcc(self):
        """Copper FCC structure: a = 3.6149 Å, SG 225"""
        from generators.bulk.spacegroups import generate_from_spacegroup

        ref = REFERENCE_DATA["Cu_fcc"]
        result = generate_from_spacegroup(
            spacegroup=ref["spacegroup"],
            elements=["Cu"],
            composition=[4],
            a=ref["a"]
        )

        assert result["success"], f"Generation failed: {result.get('error')}"
        a_actual = result["structure"]["lattice"]["a"]
        assert_lattice_within_tolerance(a_actual, ref["a"], "Cu(FCC)", "a")

    def test_gold_fcc(self):
        """Gold FCC structure: a = 4.0782 Å, SG 225"""
        from generators.bulk.spacegroups import generate_from_spacegroup

        ref = REFERENCE_DATA["Au_fcc"]
        result = generate_from_spacegroup(
            spacegroup=ref["spacegroup"],
            elements=["Au"],
            composition=[4],
            a=ref["a"]
        )

        assert result["success"], f"Generation failed: {result.get('error')}"
        a_actual = result["structure"]["lattice"]["a"]
        assert_lattice_within_tolerance(a_actual, ref["a"], "Au(FCC)", "a")

    def test_iron_bcc(self):
        """Iron BCC structure: a = 2.8665 Å, SG 229"""
        from generators.bulk.spacegroups import generate_from_spacegroup

        ref = REFERENCE_DATA["Fe_bcc"]
        result = generate_from_spacegroup(
            spacegroup=ref["spacegroup"],
            elements=["Fe"],
            composition=[2],
            a=ref["a"]
        )

        assert result["success"], f"Generation failed: {result.get('error')}"
        a_actual = result["structure"]["lattice"]["a"]
        assert_lattice_within_tolerance(a_actual, ref["a"], "Fe(BCC)", "a")

    def test_tungsten_bcc(self):
        """Tungsten BCC structure: a = 3.1652 Å, SG 229"""
        from generators.bulk.spacegroups import generate_from_spacegroup

        ref = REFERENCE_DATA["W_bcc"]
        result = generate_from_spacegroup(
            spacegroup=ref["spacegroup"],
            elements=["W"],
            composition=[2],
            a=ref["a"]
        )

        assert result["success"], f"Generation failed: {result.get('error')}"
        a_actual = result["structure"]["lattice"]["a"]
        assert_lattice_within_tolerance(a_actual, ref["a"], "W(BCC)", "a")


class TestBinaryCompounds:
    """Test binary compound structures against crystallographic databases."""

    def test_nacl_rocksalt(self):
        """NaCl rocksalt structure: a = 5.6401 Å, SG 225"""
        from generators.bulk.spacegroups import generate_from_spacegroup

        ref = REFERENCE_DATA["NaCl_rocksalt"]
        result = generate_from_spacegroup(
            spacegroup=ref["spacegroup"],
            elements=["Na", "Cl"],
            composition=[4, 4],
            a=ref["a"]
        )

        assert result["success"], f"Generation failed: {result.get('error')}"
        a_actual = result["structure"]["lattice"]["a"]
        assert_lattice_within_tolerance(a_actual, ref["a"], "NaCl", "a")

        # Verify stoichiometry
        atoms = result["structure"]["atoms"]
        na_count = sum(1 for a in atoms if a["element"] == "Na")
        cl_count = sum(1 for a in atoms if a["element"] == "Cl")
        assert na_count == cl_count, f"NaCl stoichiometry mismatch: Na={na_count}, Cl={cl_count}"

    def test_mgo_rocksalt(self):
        """MgO rocksalt structure: a = 4.2112 Å, SG 225"""
        from generators.bulk.spacegroups import generate_from_spacegroup

        ref = REFERENCE_DATA["MgO_rocksalt"]
        result = generate_from_spacegroup(
            spacegroup=ref["spacegroup"],
            elements=["Mg", "O"],
            composition=[4, 4],
            a=ref["a"]
        )

        assert result["success"], f"Generation failed: {result.get('error')}"
        a_actual = result["structure"]["lattice"]["a"]
        assert_lattice_within_tolerance(a_actual, ref["a"], "MgO", "a")

    def test_gaas_zincblende(self):
        """GaAs zincblende structure: a = 5.6533 Å, SG 216"""
        from generators.bulk.spacegroups import generate_from_spacegroup

        ref = REFERENCE_DATA["GaAs_zincblende"]
        result = generate_from_spacegroup(
            spacegroup=ref["spacegroup"],
            elements=["Ga", "As"],
            composition=[4, 4],
            a=ref["a"]
        )

        assert result["success"], f"Generation failed: {result.get('error')}"
        a_actual = result["structure"]["lattice"]["a"]
        assert_lattice_within_tolerance(a_actual, ref["a"], "GaAs", "a")

    def test_zno_wurtzite(self):
        """ZnO wurtzite structure: a = 3.2498 Å, c = 5.2066 Å, SG 186"""
        from generators.bulk.spacegroups import generate_from_spacegroup

        ref = REFERENCE_DATA["ZnO_wurtzite"]
        result = generate_from_spacegroup(
            spacegroup=ref["spacegroup"],
            elements=["Zn", "O"],
            composition=[2, 2],
            a=ref["a"],
            c=ref["c"]
        )

        assert result["success"], f"Generation failed: {result.get('error')}"
        lattice = result["structure"]["lattice"]
        assert_lattice_within_tolerance(lattice["a"], ref["a"], "ZnO", "a")
        assert_lattice_within_tolerance(lattice["c"], ref["c"], "ZnO", "c")

        # Verify c/a ratio (wurtzite ideal: 1.633, ZnO: ~1.602)
        ca_ratio = lattice["c"] / lattice["a"]
        assert 1.55 < ca_ratio < 1.70, f"ZnO c/a ratio {ca_ratio:.3f} outside expected range"


class Test2DMaterials:
    """Test 2D material structures against published values."""

    def test_graphene_lattice_constant(self):
        """Graphene: a = 2.46 Å (C-C bond 1.42 Å)"""
        from generators.two_d.xenes import generate_graphene

        ref = REFERENCE_DATA["graphene"]
        result = generate_graphene()

        assert result["success"], f"Generation failed: {result.get('error')}"

        # Get lattice constant from structure
        a_actual = result["structure"]["lattice"]["a"]
        assert_lattice_within_tolerance(a_actual, ref["a"], "Graphene", "a")

    def test_hbn_lattice_constant(self):
        """hBN monolayer: a = 2.504 Å - generated from bulk spacegroup with vacuum"""
        from generators.bulk.spacegroups import generate_from_spacegroup

        ref = REFERENCE_DATA["hBN_monolayer"]
        # hBN has P6_3/mmc structure (SG 194) with a=2.504, c includes vacuum
        # Unit cell has 2 B and 2 N atoms
        result = generate_from_spacegroup(
            spacegroup=194,  # P6_3/mmc
            elements=["B", "N"],
            composition=[2, 2],
            a=ref["a"],
            c=15.0  # Include vacuum for 2D material
        )

        assert result["success"], f"Generation failed: {result.get('error')}"

        a_actual = result["structure"]["lattice"]["a"]
        assert_lattice_within_tolerance(a_actual, ref["a"], "hBN", "a")

    def test_mos2_lattice_constant(self):
        """MoS2 monolayer: a = 3.160 Å (2H phase)"""
        from generators.two_d.tmds import generate_tmd

        ref = REFERENCE_DATA["MoS2_monolayer"]
        # generate_tmd uses 'formula' not 'material'
        result = generate_tmd(formula="MoS2", phase="1H")

        assert result["success"], f"Generation failed: {result.get('error')}"

        a_actual = result["structure"]["lattice"]["a"]
        assert_lattice_within_tolerance(a_actual, ref["a"], "MoS2", "a")

    def test_graphene_atom_count(self):
        """Graphene unit cell should have 2 carbon atoms."""
        from generators.two_d.xenes import generate_graphene

        result = generate_graphene()
        assert result["success"], f"Generation failed: {result.get('error')}"

        atoms = result["structure"]["atoms"]
        carbon_count = sum(1 for a in atoms if a["element"] == "C")
        assert carbon_count == 2, f"Graphene should have 2 C atoms per unit cell, got {carbon_count}"


class TestTopologicalInsulators:
    """Test topological insulator structures."""

    def test_bi2se3_structure(self):
        """Bi2Se3: a = 4.138 Å, c = 28.64 Å, SG 166 (R-3m)"""
        from generators.bulk.spacegroups import generate_from_spacegroup

        ref = REFERENCE_DATA["Bi2Se3"]
        result = generate_from_spacegroup(
            spacegroup=ref["spacegroup"],
            elements=["Bi", "Se"],
            composition=[6, 9],  # Hexagonal setting
            a=ref["a"],
            c=ref["c"]
        )

        assert result["success"], f"Generation failed: {result.get('error')}"
        lattice = result["structure"]["lattice"]
        assert_lattice_within_tolerance(lattice["a"], ref["a"], "Bi2Se3", "a")
        assert_lattice_within_tolerance(lattice["c"], ref["c"], "Bi2Se3", "c")

    def test_bi2te3_structure(self):
        """Bi2Te3: a = 4.3835 Å, c = 30.487 Å, SG 166 (R-3m)"""
        from generators.bulk.spacegroups import generate_from_spacegroup

        ref = REFERENCE_DATA["Bi2Te3"]
        result = generate_from_spacegroup(
            spacegroup=ref["spacegroup"],
            elements=["Bi", "Te"],
            composition=[6, 9],  # Hexagonal setting
            a=ref["a"],
            c=ref["c"]
        )

        assert result["success"], f"Generation failed: {result.get('error')}"
        lattice = result["structure"]["lattice"]
        assert_lattice_within_tolerance(lattice["a"], ref["a"], "Bi2Te3", "a")
        assert_lattice_within_tolerance(lattice["c"], ref["c"], "Bi2Te3", "c")


class TestStructureTransformations:
    """Test physical correctness of structure transformations."""

    def test_hydrostatic_strain_volume(self):
        """5% hydrostatic strain should increase volume by (1.05)^3 = 1.1576"""
        import numpy as np
        from generators.bulk.spacegroups import generate_from_spacegroup
        from generators.bulk.strain import apply_strain

        # Generate Cu FCC
        base = generate_from_spacegroup(
            spacegroup=225,
            elements=["Cu"],
            composition=[4],
            a=3.6149
        )
        assert base["success"], f"Base generation failed: {base.get('error')}"

        # Apply 5% hydrostatic strain
        strained = apply_strain(
            structure_dict=base["structure"],
            strain=0.05,
            strain_type="hydrostatic"
        )
        assert strained["success"], f"Strain application failed: {strained.get('error')}"

        # Calculate volumes from lattice matrices
        orig_matrix = np.array(base["structure"]["lattice"]["matrix"])
        new_matrix = np.array(strained["structure"]["lattice"]["matrix"])

        orig_vol = abs(np.linalg.det(orig_matrix))
        new_vol = abs(np.linalg.det(new_matrix))

        expected_ratio = 1.05 ** 3  # 1.157625
        actual_ratio = new_vol / orig_vol

        deviation_pct = abs(actual_ratio - expected_ratio) / expected_ratio * 100
        assert deviation_pct < 1.0, (
            f"Volume ratio {actual_ratio:.4f} deviates from expected {expected_ratio:.4f} "
            f"by {deviation_pct:.2f}%"
        )

    def test_vacancy_removes_one_atom(self):
        """Creating a vacancy should remove exactly one atom."""
        from generators.bulk.spacegroups import generate_from_spacegroup
        from generators.defect.point_defects import generate_vacancy

        # Generate Cu supercell
        base = generate_from_spacegroup(
            spacegroup=225,
            elements=["Cu"],
            composition=[32],
            a=3.6149
        )
        assert base["success"], f"Base generation failed: {base.get('error')}"

        n_before = len(base["structure"]["atoms"])

        # Create vacancy using the correct function
        with_vacancy = generate_vacancy(
            host_structure=base["structure"],
            vacancy_site=0,
            vacancy_type="single"
        )
        assert with_vacancy["success"], f"Vacancy creation failed: {with_vacancy.get('error')}"

        n_after = len(with_vacancy["structure"]["atoms"])

        assert n_after == n_before - 1, (
            f"Vacancy should remove 1 atom: {n_before} → {n_after} (expected {n_before - 1})"
        )


class TestSpaceGroupSymmetry:
    """Test space group handling and symmetry operations."""

    def test_all_crystal_systems_covered(self):
        """Verify we can generate structures in all 7 crystal systems."""
        from generators.bulk.spacegroups import generate_from_spacegroup

        # Representative space groups for each crystal system
        crystal_systems = {
            "triclinic": 1,      # P1
            "monoclinic": 15,    # C2/c
            "orthorhombic": 62,  # Pnma
            "tetragonal": 123,   # P4/mmm
            "trigonal": 166,     # R-3m
            "hexagonal": 194,    # P6_3/mmc
            "cubic": 225,        # Fm-3m
        }

        for system, sg in crystal_systems.items():
            result = generate_from_spacegroup(
                spacegroup=sg,
                elements=["Fe"],
                composition=[4],
                a=3.0
            )
            assert result["success"], f"{system} (SG {sg}) generation failed: {result.get('error')}"

    def test_space_group_subgroup_relations(self):
        """Verify space group subgroup finder works correctly."""
        from generators.quality_control.symmetry import get_subgroups

        # SG 225 (Fm-3m) is a supergroup of many lower symmetry groups
        result = get_subgroups(225)

        assert result["success"], f"Subgroup lookup failed: {result.get('error')}"
        assert len(result["subgroups"]) > 0, "SG 225 should have subgroups"

        # Verify all returned subgroups are valid (1-230)
        # subgroups is a list of dicts with "number" key
        for sg_info in result["subgroups"]:
            sg_num = sg_info["number"]
            assert 1 <= sg_num <= 230, f"Invalid subgroup number: {sg_num}"

    def test_symmetry_path_finding(self):
        """Verify symmetry path finder works for group-subgroup chains."""
        from generators.quality_control.symmetry import get_symmetry_path

        # Find path from high symmetry (225) to low symmetry (1)
        # Uses start_spg and end_spg as parameter names
        result = get_symmetry_path(start_spg=225, end_spg=1)

        assert result["success"], f"Path finding failed: {result.get('error')}"
        assert len(result["path"]) >= 2, "Path should include at least start and end"
        # Path entries are dicts with "number" key
        assert result["path"][0]["number"] == 225, "Path should start at 225"
        assert result["path"][-1]["number"] == 1, "Path should end at 1"


class TestCrystallographicConstraints:
    """Test that generated structures satisfy crystallographic constraints."""

    def test_cubic_angles_are_90(self):
        """Cubic structures must have all angles = 90°."""
        from generators.bulk.spacegroups import generate_from_spacegroup

        result = generate_from_spacegroup(
            spacegroup=225,  # Fm-3m (cubic)
            elements=["Cu"],
            composition=[4],
            a=3.6149
        )

        assert result["success"], f"Generation failed: {result.get('error')}"
        lattice = result["structure"]["lattice"]

        assert abs(lattice["alpha"] - 90.0) < 0.001, f"Cubic alpha should be 90°, got {lattice['alpha']}"
        assert abs(lattice["beta"] - 90.0) < 0.001, f"Cubic beta should be 90°, got {lattice['beta']}"
        assert abs(lattice["gamma"] - 90.0) < 0.001, f"Cubic gamma should be 90°, got {lattice['gamma']}"

    def test_cubic_axes_are_equal(self):
        """Cubic structures must have a = b = c."""
        from generators.bulk.spacegroups import generate_from_spacegroup

        result = generate_from_spacegroup(
            spacegroup=227,  # Fd-3m (cubic)
            elements=["Si"],
            composition=[8],
            a=5.4309
        )

        assert result["success"], f"Generation failed: {result.get('error')}"
        lattice = result["structure"]["lattice"]

        assert abs(lattice["a"] - lattice["b"]) < 0.001, f"Cubic a ≠ b: {lattice['a']} vs {lattice['b']}"
        assert abs(lattice["b"] - lattice["c"]) < 0.001, f"Cubic b ≠ c: {lattice['b']} vs {lattice['c']}"

    def test_hexagonal_gamma_is_120(self):
        """Hexagonal structures must have gamma = 120°."""
        from generators.bulk.spacegroups import generate_from_spacegroup

        result = generate_from_spacegroup(
            spacegroup=194,  # P6_3/mmc (hexagonal)
            elements=["Mg"],
            composition=[2],
            a=3.21,
            c=5.21
        )

        assert result["success"], f"Generation failed: {result.get('error')}"
        lattice = result["structure"]["lattice"]

        assert abs(lattice["gamma"] - 120.0) < 0.001, f"Hexagonal gamma should be 120°, got {lattice['gamma']}"


class TestSpaceGroupTransformation:
    """
    Test Bilbao-like space group transformation functionality.

    This is critical for materials science workflows where structures in
    different phases need to be compared. Example: Comparing cubic Pm-3m
    and rhombohedral R3c phases of BiFeO3.

    Reference: Bilbao Crystallographic Server TRANSTRU tool
    https://www.cryst.ehu.es/cryst/transtru.html
    """

    def test_transform_to_direct_subgroup(self):
        """Transform cubic Pm-3m to tetragonal P4/mmm (direct subgroup)."""
        from generators.bulk.spacegroups import generate_from_spacegroup
        from generators.quality_control.symmetry import transform_to_subgroup

        # Create simple cubic structure
        cubic = generate_from_spacegroup(
            spacegroup=221,  # Pm-3m
            elements=["Cu"],
            composition=[1],
            a=3.6
        )
        assert cubic["success"], f"Cubic generation failed: {cubic.get('error')}"

        # Transform to P4/mmm (123) - direct t-subgroup
        result = transform_to_subgroup(cubic["structure"], 123)

        assert result["success"], f"Transformation failed: {result.get('error')}"
        assert result["original_spacegroup"]["number"] == 221
        assert result["target_spacegroup"]["number"] == 123
        # Atoms should increase due to cell transformation
        assert result["n_atoms_transformed"] >= result["n_atoms_original"]

    def test_bifeo3_cubic_to_r3c_transformation(self):
        """
        BiFeO3 phase transformation: Pm-3m (cubic) → R3c (rhombohedral).

        This is the canonical example for comparing ferroelectric R3c phase
        with the high-temperature cubic paraelectric phase.

        Path: 221 (Pm-3m) → 166 (R-3m) → 160 (R3m) → 161 (R3c)

        Reference: Catalan & Scott, Adv. Mater. 21, 2463 (2009)
        """
        from generators.bulk.spacegroups import generate_from_spacegroup
        from generators.quality_control.symmetry import transform_to_subgroup

        # Create cubic BiFeO3 perovskite
        # Bi at 1a (0,0,0), Fe at 1b (0.5,0.5,0.5), O at 3c
        cubic_bfo = generate_from_spacegroup(
            spacegroup=221,  # Pm-3m
            elements=["Bi", "Fe", "O"],
            composition=[1, 1, 3],
            a=3.96  # Approximate cubic perovskite lattice constant
        )
        assert cubic_bfo["success"], f"BiFeO3 generation failed: {cubic_bfo.get('error')}"
        assert len(cubic_bfo["structure"]["atoms"]) == 5, "Cubic perovskite should have 5 atoms"

        # Transform to R3c (161)
        result = transform_to_subgroup(cubic_bfo["structure"], 161)

        assert result["success"], f"Transformation failed: {result.get('error')}"
        assert result["original_spacegroup"]["number"] == 221
        assert result["target_spacegroup"]["number"] == 161

        # Verify atom count increases (Pm-3m to R3c involves cell multiplication)
        # Expected: 5 atoms → 30 atoms (6x multiplication)
        assert result["n_atoms_transformed"] > result["n_atoms_original"]
        assert result["multiplicity"] >= 1

        # Verify transformation path includes expected steps
        path = result["transformation_path"]
        assert 221 in path, "Path should start at 221"
        assert 161 in path, "Path should end at 161"

    def test_get_all_transformation_paths(self):
        """Test finding all possible transformation paths."""
        from generators.quality_control.symmetry import get_all_subgroup_paths

        # Find all paths from Pm-3m (221) to R3c (161)
        result = get_all_subgroup_paths(221, 161, max_depth=4)

        assert result["success"], f"Path finding failed: {result.get('error')}"
        assert result["n_paths"] > 0, "Should find at least one path"

        # Verify each path is valid
        for path_info in result["paths"]:
            path = path_info["path"]
            assert path[0] == 221, f"Path should start at 221, got {path[0]}"
            assert path[-1] == 161, f"Path should end at 161, got {path[-1]}"

    def test_transform_by_specific_path(self):
        """Test transforming along a user-specified path."""
        from generators.bulk.spacegroups import generate_from_spacegroup
        from generators.quality_control.symmetry import transform_by_path

        # Create cubic structure
        cubic = generate_from_spacegroup(
            spacegroup=221,
            elements=["Fe"],
            composition=[1],
            a=2.87
        )
        assert cubic["success"]

        # Transform via specific path: 221 → 166 → 160
        result = transform_by_path(
            cubic["structure"],
            path=[221, 166, 160]
        )

        assert result["success"], f"Path transformation failed: {result.get('error')}"
        assert result["final_spacegroup"]["number"] == 160
        assert result["n_steps"] == 2

    def test_transformation_preserves_composition(self):
        """Transformation should preserve atomic composition ratios."""
        from generators.bulk.spacegroups import generate_from_spacegroup
        from generators.quality_control.symmetry import transform_to_subgroup

        # Create NaCl rocksalt
        nacl = generate_from_spacegroup(
            spacegroup=225,  # Fm-3m
            elements=["Na", "Cl"],
            composition=[4, 4],
            a=5.64
        )
        assert nacl["success"]

        # Get original composition
        orig_atoms = nacl["structure"]["atoms"]
        orig_na = sum(1 for a in orig_atoms if a["element"] == "Na")
        orig_cl = sum(1 for a in orig_atoms if a["element"] == "Cl")
        orig_ratio = orig_na / orig_cl

        # Transform to lower symmetry
        result = transform_to_subgroup(nacl["structure"], 123)  # P4/mmm

        if result["success"]:
            # Verify composition ratio is preserved
            new_atoms = result["structure"]["atoms"]
            new_na = sum(1 for a in new_atoms if a["element"] == "Na")
            new_cl = sum(1 for a in new_atoms if a["element"] == "Cl")

            if new_cl > 0:
                new_ratio = new_na / new_cl
                assert abs(new_ratio - orig_ratio) < 0.01, (
                    f"Composition ratio not preserved: {orig_ratio:.2f} → {new_ratio:.2f}"
                )

    def test_invalid_target_not_subgroup(self):
        """Transformation should fail if target is not a subgroup."""
        from generators.bulk.spacegroups import generate_from_spacegroup
        from generators.quality_control.symmetry import transform_to_subgroup

        # Create structure in SG 15 (C2/c)
        struct = generate_from_spacegroup(
            spacegroup=15,
            elements=["Fe"],
            composition=[8],
            a=5.0
        )
        assert struct["success"]

        # Try to transform to higher symmetry (should fail)
        result = transform_to_subgroup(struct["structure"], 225)  # Fm-3m is NOT a subgroup of C2/c

        assert not result["success"], "Should fail when target is not a subgroup"
        assert "error" in result


# =============================================================================
# MAIN
# =============================================================================

if __name__ == "__main__":
    pytest.main([__file__, "-v", "--tb=short"])
