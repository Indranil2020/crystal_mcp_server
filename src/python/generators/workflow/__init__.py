"""
Workflow Module for Crystal Structure Generation

Provides high-level workflow capabilities for LLM-driven structure generation:
- Parametric scans (z-scan, distance scan, rotation scan)
- Iterative structure editing with undo/redo
- Trajectory/animation generation for visualization
- Multi-structure operations and sweeps

All functions are designed to be called by an LLM through MCP with clear
parameter specifications and comprehensive return data.

Example usage:
    # Z-scan of molecule on surface
    generate_z_scan(
        surface_material="NaCl",
        surface_miller=(1, 0, 0),
        molecule="PTCDA",
        z_start=3.0,
        z_end=20.0,
        z_step=0.1
    )

    # Iterative editing - chain operations on structure
    edit_structure(
        structure=base_structure,
        operations=[
            {"operation": "make_supercell", "scaling": [2, 2, 2]},
            {"operation": "create_vacancy", "site_index": 0},
            {"operation": "apply_strain", "strain": 0.02, "strain_type": "biaxial"}
        ]
    )
"""

# Parametric Scans
from generators.workflow.parametric_scans import (
    # Basic scans
    generate_z_scan,
    generate_distance_scan,
    generate_rotation_scan,
    # Directional and spatial scans
    generate_directional_scan,
    generate_planar_scan,
    generate_radial_scan,
    # Lattice and strain scans
    generate_lattice_scan,
    generate_strain_tensor_scan,
    # Surface/adsorption scans
    generate_coverage_scan,
    # Reaction path scans
    generate_neb_path,
    generate_bond_scan,
    generate_dihedral_scan,
    # General parameter sweeps
    generate_parameter_sweep,
    generate_multi_parameter_scan,
    # Animation and utilities
    generate_trajectory_animation,
    generate_ptcda_nacl_zscan,
    get_available_scans,
)

# Iterative Editing
from generators.workflow.iterative_editing import (
    StructureEditor,
    edit_structure,
    compare_structures,
    get_available_edit_operations,
)

__all__ = [
    # Basic scans
    "generate_z_scan",
    "generate_distance_scan",
    "generate_rotation_scan",
    # Directional and spatial scans
    "generate_directional_scan",
    "generate_planar_scan",
    "generate_radial_scan",
    # Lattice and strain scans
    "generate_lattice_scan",
    "generate_strain_tensor_scan",
    # Surface/adsorption scans
    "generate_coverage_scan",
    # Reaction path scans
    "generate_neb_path",
    "generate_bond_scan",
    "generate_dihedral_scan",
    # General parameter sweeps
    "generate_parameter_sweep",
    "generate_multi_parameter_scan",
    # Animation and utilities
    "generate_trajectory_animation",
    "generate_ptcda_nacl_zscan",
    "get_available_scans",
    # Iterative Editing
    "StructureEditor",
    "edit_structure",
    "compare_structures",
    "get_available_edit_operations",
]
