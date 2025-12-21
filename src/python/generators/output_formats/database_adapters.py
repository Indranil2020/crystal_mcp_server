"""
output_formats/database_adapters.py - Database Integration Adapters

Provides export functionality for major computational materials databases:
- ASE (Atomic Simulation Environment)
- AiiDA (Automated Interactive Infrastructure and Database for Computational Science)
- AFLOW (Automatic Flow for Materials Discovery)
- Materials Project JSON format
"""

from typing import Dict, Any, List, Optional, Union
import json
import numpy as np
from datetime import datetime


def _structure_to_dict(structure) -> Dict[str, Any]:
    """Convert pymatgen Structure to serializable dict."""
    if hasattr(structure, 'as_dict'):
        return structure.as_dict()
    return structure


def export_ase_atoms(
    structure,
    include_constraints: bool = False,
    include_calculator: bool = False,
    info: Optional[Dict[str, Any]] = None
) -> Dict[str, Any]:
    """
    Export structure to ASE Atoms JSON format.

    ASE Atoms format is widely used in the atomistic simulation community
    and can be directly loaded into ASE for further calculations.

    Args:
        structure: Crystal structure (pymatgen Structure or dict)
        include_constraints: Include constraint information
        include_calculator: Include calculator settings placeholder
        info: Additional metadata to include in atoms.info

    Returns:
        ASE-compatible JSON structure
    """
    from pymatgen.core import Structure

    # Handle input structure
    if isinstance(structure, dict):
        if "lattice" in structure and "sites" in structure:
            struct = Structure.from_dict(structure)
        elif "lattice" in structure and "atoms" in structure:
            from pymatgen.core import Lattice
            lat = Lattice(structure["lattice"]["matrix"])
            species = [s["element"] for s in structure["atoms"]]
            coords = [s["coords"] for s in structure["atoms"]]
            struct = Structure(lat, species, coords)
        else:
            return {"success": False, "error": {"code": "INVALID_STRUCTURE",
                    "message": "Cannot parse structure dict"}}
    else:
        struct = structure

    # Build ASE Atoms compatible dict
    cell = struct.lattice.matrix.tolist()
    positions = [site.coords.tolist() for site in struct.sites]
    symbols = [str(site.specie) for site in struct.sites]
    numbers = [site.specie.Z for site in struct.sites]

    # Calculate scaled positions (fractional coordinates)
    scaled_positions = [site.frac_coords.tolist() for site in struct.sites]

    # Build the ASE-compatible JSON
    ase_dict = {
        "cell": cell,
        "pbc": [True, True, True],  # Periodic boundary conditions
        "positions": positions,
        "numbers": numbers,
        "symbols": symbols,
        "scaled_positions": scaled_positions,
    }

    # Add optional info dict
    if info:
        ase_dict["info"] = info
    else:
        ase_dict["info"] = {
            "source": "crystal-mcp-server",
            "export_time": datetime.now().isoformat(),
            "formula": struct.composition.reduced_formula
        }

    # Add constraints placeholder
    if include_constraints:
        ase_dict["constraints"] = []

    # Add calculator placeholder
    if include_calculator:
        ase_dict["calculator"] = {
            "name": None,
            "parameters": {}
        }

    # Generate ASE-loadable Python code
    ase_code = f'''# ASE Atoms loading code
from ase import Atoms
import numpy as np

cell = {cell}
positions = {positions}
symbols = {symbols}
pbc = [True, True, True]

atoms = Atoms(
    symbols=symbols,
    positions=positions,
    cell=cell,
    pbc=pbc
)
'''

    return {
        "success": True,
        "format": "ase_atoms",
        "ase_dict": ase_dict,
        "ase_code": ase_code,
        "n_atoms": len(symbols),
        "formula": struct.composition.reduced_formula,
        "volume": struct.lattice.volume
    }


def export_aiida_structuredata(
    structure,
    node_label: Optional[str] = None,
    node_description: Optional[str] = None,
    extras: Optional[Dict[str, Any]] = None
) -> Dict[str, Any]:
    """
    Export structure to AiiDA StructureData format.

    AiiDA is a workflow management system for computational materials science.
    This export creates a dictionary that can be used to create a StructureData node.

    Args:
        structure: Crystal structure (pymatgen Structure or dict)
        node_label: Label for the AiiDA node
        node_description: Description for the AiiDA node
        extras: Extra metadata to attach to node

    Returns:
        AiiDA StructureData compatible format
    """
    from pymatgen.core import Structure

    # Handle input structure
    if isinstance(structure, dict):
        if "lattice" in structure and "sites" in structure:
            struct = Structure.from_dict(structure)
        elif "lattice" in structure and "atoms" in structure:
            from pymatgen.core import Lattice
            lat = Lattice(structure["lattice"]["matrix"])
            species = [s["element"] for s in structure["atoms"]]
            coords = [s["coords"] for s in structure["atoms"]]
            struct = Structure(lat, species, coords)
        else:
            return {"success": False, "error": {"code": "INVALID_STRUCTURE",
                    "message": "Cannot parse structure dict"}}
    else:
        struct = structure

    # AiiDA StructureData uses cell and sites format
    cell = struct.lattice.matrix.tolist()

    # Build sites list in AiiDA format
    sites = []
    kinds = {}  # Track unique kinds

    for i, site in enumerate(struct.sites):
        symbol = str(site.specie)

        # Track kinds (AiiDA's way of handling species)
        if symbol not in kinds:
            kinds[symbol] = {
                "name": symbol,
                "symbols": [symbol],
                "weights": [1.0],
                "mass": site.specie.atomic_mass
            }

        sites.append({
            "position": site.coords.tolist(),
            "kind_name": symbol
        })

    # Build AiiDA-compatible structure
    aiida_dict = {
        "cell": cell,
        "pbc1": True,
        "pbc2": True,
        "pbc3": True,
        "kinds": list(kinds.values()),
        "sites": sites
    }

    # Build node attributes
    node_attrs = {
        "label": node_label or struct.composition.reduced_formula,
        "description": node_description or f"Structure exported from crystal-mcp-server",
    }

    if extras:
        node_attrs["extras"] = extras
    else:
        node_attrs["extras"] = {
            "source": "crystal-mcp-server",
            "export_time": datetime.now().isoformat(),
            "spacegroup": struct.get_space_group_info()[0] if hasattr(struct, 'get_space_group_info') else None
        }

    # Generate AiiDA loading code
    aiida_code = f'''# AiiDA StructureData creation code
from aiida.orm import StructureData
from aiida.plugins import DataFactory
import numpy as np

StructureData = DataFactory('core.structure')

structure = StructureData(cell={cell})

# Add kinds and sites
'''
    for kind in kinds.values():
        aiida_code += f"structure.append_kind(Kind(name='{kind['name']}', symbols={kind['symbols']}))\n"

    for i, site in enumerate(sites):
        aiida_code += f"structure.append_site(Site(kind_name='{site['kind_name']}', position={site['position']}))\n"

    aiida_code += f"\nstructure.label = '{node_attrs['label']}'\n"
    aiida_code += f"structure.description = '{node_attrs['description']}'\n"

    return {
        "success": True,
        "format": "aiida_structuredata",
        "aiida_dict": aiida_dict,
        "node_attributes": node_attrs,
        "aiida_code": aiida_code,
        "n_atoms": len(sites),
        "n_kinds": len(kinds),
        "formula": struct.composition.reduced_formula
    }


def export_aflow(
    structure,
    auid: Optional[str] = None,
    compound: Optional[str] = None,
    keywords: Optional[List[str]] = None
) -> Dict[str, Any]:
    """
    Export structure to AFLOW AFLUX format.

    AFLOW is a high-throughput framework for materials discovery.
    This creates an AFLOW-compatible JSON entry.

    Args:
        structure: Crystal structure (pymatgen Structure or dict)
        auid: AFLOW Unique Identifier (auto-generated if not provided)
        compound: Compound name
        keywords: Search keywords

    Returns:
        AFLOW-compatible JSON structure
    """
    from pymatgen.core import Structure
    from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

    # Handle input structure
    if isinstance(structure, dict):
        if "lattice" in structure and "sites" in structure:
            struct = Structure.from_dict(structure)
        elif "lattice" in structure and "atoms" in structure:
            from pymatgen.core import Lattice
            lat = Lattice(structure["lattice"]["matrix"])
            species = [s["element"] for s in structure["atoms"]]
            coords = [s["coords"] for s in structure["atoms"]]
            struct = Structure(lat, species, coords)
        else:
            return {"success": False, "error": {"code": "INVALID_STRUCTURE",
                    "message": "Cannot parse structure dict"}}
    else:
        struct = structure

    # Get symmetry information
    sga = SpacegroupAnalyzer(struct, symprec=0.1)
    spg_symbol = sga.get_space_group_symbol()
    spg_number = sga.get_space_group_number()
    crystal_system = sga.get_crystal_system()

    # Get lattice parameters
    a, b, c = struct.lattice.abc
    alpha, beta, gamma = struct.lattice.angles

    # Build AFLOW entry
    aflow_entry = {
        # Identifiers
        "auid": auid or f"aflow:generated:{datetime.now().strftime('%Y%m%d%H%M%S')}",
        "aurl": f"aflowlib.org/generated",
        "compound": compound or struct.composition.reduced_formula,

        # Composition
        "species": list(set(str(s.specie) for s in struct.sites)),
        "species_pp": list(set(str(s.specie) for s in struct.sites)),  # Pseudopotential species
        "composition": [float(x) for x in struct.composition.get_el_amt_dict().values()],
        "natoms": len(struct.sites),
        "stoichiometry": struct.composition.reduced_formula,

        # Lattice
        "geometry": [a, b, c, alpha, beta, gamma],
        "lattice_system_relax": crystal_system,
        "volume_cell": struct.lattice.volume,
        "volume_atom": struct.lattice.volume / len(struct.sites),
        "density": struct.density,

        # Symmetry
        "spacegroup_relax": spg_symbol,
        "sg": spg_number,
        "sg2": spg_number,
        "Bravais_lattice_relax": sga.get_lattice_type(),
        "crystal_class": crystal_system,
        "point_group_relax": sga.get_point_group_symbol(),

        # Positions
        "positions_cartesian": [site.coords.tolist() for site in struct.sites],
        "positions_fractional": [site.frac_coords.tolist() for site in struct.sites],

        # Keywords
        "keywords": keywords or ["generated", "crystal-mcp-server"],

        # Metadata
        "catalog": "generated",
        "data_api": "1.0",
        "data_source": "crystal-mcp-server"
    }

    # Generate AFLUX query URL (for reference)
    aflux_url = f"http://aflowlib.org/API/aflux/?species({','.join(aflow_entry['species'])})"

    return {
        "success": True,
        "format": "aflow",
        "aflow_entry": aflow_entry,
        "aflux_url": aflux_url,
        "formula": struct.composition.reduced_formula,
        "spacegroup": f"{spg_symbol} (#{spg_number})"
    }


def export_materials_project(
    structure,
    material_id: Optional[str] = None,
    task_id: Optional[str] = None,
    deprecated: bool = False
) -> Dict[str, Any]:
    """
    Export structure to Materials Project JSON format.

    Compatible with Materials Project API v2 format.

    Args:
        structure: Crystal structure (pymatgen Structure or dict)
        material_id: MP material ID (auto-generated if not provided)
        task_id: Task ID
        deprecated: Mark as deprecated

    Returns:
        Materials Project compatible JSON
    """
    from pymatgen.core import Structure
    from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

    # Handle input structure
    if isinstance(structure, dict):
        if "lattice" in structure and "sites" in structure:
            struct = Structure.from_dict(structure)
        elif "lattice" in structure and "atoms" in structure:
            from pymatgen.core import Lattice
            lat = Lattice(structure["lattice"]["matrix"])
            species = [s["element"] for s in structure["atoms"]]
            coords = [s["coords"] for s in structure["atoms"]]
            struct = Structure(lat, species, coords)
        else:
            return {"success": False, "error": {"code": "INVALID_STRUCTURE",
                    "message": "Cannot parse structure dict"}}
    else:
        struct = structure

    # Get symmetry
    sga = SpacegroupAnalyzer(struct, symprec=0.1)

    # Build MP format
    mp_entry = {
        "material_id": material_id or f"mp-generated-{datetime.now().strftime('%Y%m%d%H%M%S')}",
        "task_id": task_id,
        "deprecated": deprecated,

        # Formula
        "formula_pretty": struct.composition.reduced_formula,
        "formula_anonymous": struct.composition.anonymized_formula,
        "chemsys": "-".join(sorted(set(str(s.specie) for s in struct.sites))),
        "nelements": len(set(str(s.specie) for s in struct.sites)),
        "elements": sorted(set(str(s.specie) for s in struct.sites)),
        "nsites": len(struct.sites),

        # Composition
        "composition": struct.composition.as_dict(),
        "composition_reduced": struct.composition.to_reduced_dict,

        # Structure (full pymatgen format)
        "structure": struct.as_dict(),

        # Symmetry
        "symmetry": {
            "crystal_system": sga.get_crystal_system(),
            "symbol": sga.get_space_group_symbol(),
            "number": sga.get_space_group_number(),
            "point_group": sga.get_point_group_symbol(),
            "hall": sga.get_hall()
        },

        # Properties
        "volume": struct.lattice.volume,
        "density": struct.density,
        "density_atomic": len(struct.sites) / struct.lattice.volume,

        # Metadata
        "last_updated": datetime.now().isoformat(),
        "origins": [{"name": "crystal-mcp-server", "task_id": task_id}],
        "database_IDs": {}
    }

    return {
        "success": True,
        "format": "materials_project",
        "mp_entry": mp_entry,
        "material_id": mp_entry["material_id"],
        "formula": struct.composition.reduced_formula
    }


def export_optimade(
    structure,
    entry_id: Optional[str] = None,
    provider: str = "crystal-mcp-server"
) -> Dict[str, Any]:
    """
    Export structure to OPTIMADE JSON format.

    OPTIMADE (Open Databases Integration for Materials Design) is a
    standardized API for materials databases.

    Args:
        structure: Crystal structure (pymatgen Structure or dict)
        entry_id: Entry identifier
        provider: Database provider name

    Returns:
        OPTIMADE-compatible JSON structure
    """
    from pymatgen.core import Structure
    from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

    # Handle input structure
    if isinstance(structure, dict):
        if "lattice" in structure and "sites" in structure:
            struct = Structure.from_dict(structure)
        elif "lattice" in structure and "atoms" in structure:
            from pymatgen.core import Lattice
            lat = Lattice(structure["lattice"]["matrix"])
            species = [s["element"] for s in structure["atoms"]]
            coords = [s["coords"] for s in structure["atoms"]]
            struct = Structure(lat, species, coords)
        else:
            return {"success": False, "error": {"code": "INVALID_STRUCTURE",
                    "message": "Cannot parse structure dict"}}
    else:
        struct = structure

    # Get elements
    elements = sorted(set(str(s.specie) for s in struct.sites))
    element_ratios = []
    comp_dict = struct.composition.get_el_amt_dict()
    total = sum(comp_dict.values())
    for el in elements:
        element_ratios.append(comp_dict.get(el, 0) / total)

    # Get symmetry
    sga = SpacegroupAnalyzer(struct, symprec=0.1)

    # Build OPTIMADE entry
    optimade_entry = {
        "id": entry_id or f"{provider}:{datetime.now().strftime('%Y%m%d%H%M%S')}",
        "type": "structures",
        "attributes": {
            # Required properties
            "elements": elements,
            "nelements": len(elements),
            "elements_ratios": element_ratios,
            "chemical_formula_descriptive": struct.composition.reduced_formula,
            "chemical_formula_reduced": struct.composition.reduced_formula,
            "chemical_formula_anonymous": struct.composition.anonymized_formula,
            "dimension_types": [1, 1, 1],  # 3D periodic
            "nperiodic_dimensions": 3,
            "lattice_vectors": struct.lattice.matrix.tolist(),
            "cartesian_site_positions": [site.coords.tolist() for site in struct.sites],
            "nsites": len(struct.sites),
            "species_at_sites": [str(site.specie) for site in struct.sites],
            "species": [
                {"name": el, "chemical_symbols": [el], "concentration": [1.0]}
                for el in elements
            ],

            # Optional properties
            "space_group_symbol": sga.get_space_group_symbol(),
            "space_group_it_number": sga.get_space_group_number(),

            # Provider-specific
            "_crystal_mcp_server_source": "generated",
            "_crystal_mcp_server_timestamp": datetime.now().isoformat()
        }
    }

    return {
        "success": True,
        "format": "optimade",
        "optimade_entry": optimade_entry,
        "entry_id": optimade_entry["id"],
        "formula": struct.composition.reduced_formula
    }


def export_jarvis(
    structure,
    jid: Optional[str] = None
) -> Dict[str, Any]:
    """
    Export structure to JARVIS-DFT format.

    JARVIS (Joint Automated Repository for Various Integrated Simulations)
    is a database for materials science.

    Args:
        structure: Crystal structure
        jid: JARVIS ID

    Returns:
        JARVIS-compatible JSON structure
    """
    from pymatgen.core import Structure
    from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

    # Handle input structure
    if isinstance(structure, dict):
        if "lattice" in structure and "sites" in structure:
            struct = Structure.from_dict(structure)
        elif "lattice" in structure and "atoms" in structure:
            from pymatgen.core import Lattice
            lat = Lattice(structure["lattice"]["matrix"])
            species = [s["element"] for s in structure["atoms"]]
            coords = [s["coords"] for s in structure["atoms"]]
            struct = Structure(lat, species, coords)
        else:
            return {"success": False, "error": {"code": "INVALID_STRUCTURE",
                    "message": "Cannot parse structure dict"}}
    else:
        struct = structure

    # Get symmetry
    sga = SpacegroupAnalyzer(struct, symprec=0.1)

    # Build JARVIS format (Atoms-like)
    jarvis_entry = {
        "jid": jid or f"JVASP-generated-{datetime.now().strftime('%Y%m%d%H%M%S')}",
        "atoms": {
            "lattice_mat": struct.lattice.matrix.tolist(),
            "coords": [site.frac_coords.tolist() for site in struct.sites],
            "elements": [str(site.specie) for site in struct.sites],
            "cartesian": False,
            "props": ["" for _ in struct.sites]
        },
        "formula": struct.composition.reduced_formula,
        "spg_number": sga.get_space_group_number(),
        "spg_symbol": sga.get_space_group_symbol(),
        "crys": sga.get_crystal_system(),
        "natoms": len(struct.sites),
        "source": "crystal-mcp-server"
    }

    return {
        "success": True,
        "format": "jarvis",
        "jarvis_entry": jarvis_entry,
        "jid": jarvis_entry["jid"],
        "formula": struct.composition.reduced_formula
    }
