"""
molecule/universal_molecule.py - Universal Molecule Generation

Generate 3D molecular structures from any identifier:
- Common names (e.g., "aspirin", "PTCDA")
- SMILES strings (e.g., "CC(=O)OC1=CC=CC=C1C(=O)O")
- IUPAC names (e.g., "perylene-3,4,9,10-tetracarboxylic dianhydride")
- PubChem CIDs (e.g., 2244)
- InChI strings

Uses RDKit for 3D conformer generation and PubChem API for name resolution.
"""

from typing import Dict, Any, Optional, Tuple, List, Union
import logging
import re
import importlib.util

# Check module availability using importlib (no try/except)
RDKIT_AVAILABLE = importlib.util.find_spec("rdkit") is not None
REQUESTS_AVAILABLE = importlib.util.find_spec("requests") is not None
PYMATGEN_AVAILABLE = importlib.util.find_spec("pymatgen") is not None

# Initialize module references
Chem = None
AllChem = None
Descriptors = None
requests = None
Molecule = None

# Import modules only if available
if RDKIT_AVAILABLE:
    from rdkit import Chem
    from rdkit.Chem import AllChem, Descriptors

if REQUESTS_AVAILABLE:
    import requests

if PYMATGEN_AVAILABLE:
    from pymatgen.core import Molecule

# Import local database for fallback
from .small_molecules import MOLECULE_DATABASE

# Import unified molecule database (SQLite-based)
from .molecule_database import (
    MoleculeDatabase, 
    get_molecule_database, 
    lookup_molecule as db_lookup_molecule
)

logger = logging.getLogger(__name__)


# DEPRECATED: Legacy aliases kept for backward compatibility
# All new molecules should be added to the SQLite database via db_manager
MOLECULE_ALIASES = {
    # Populated by user configuration if needed.
    # Default is empty to rely on robust database lookups.
}


def detect_input_type(identifier: str) -> str:
    """
    Auto-detect the type of molecular identifier.
    
    Args:
        identifier: The molecular identifier string
        
    Returns:
        One of: "smiles", "inchi", "cid", "iupac", "name"
    """
    identifier = identifier.strip()
    
    # Check for InChI
    if identifier.startswith("InChI="):
        return "inchi"
    
    # Check for PubChem CID (pure numeric)
    if identifier.isdigit():
        return "cid"
    
    # Check for SMILES patterns (contains typical SMILES characters)
    smiles_indicators = ['=', '#', '@', '/', '\\', '[', ']', '(', ')']
    if any(char in identifier for char in smiles_indicators):
        # Validate it's actually a parseable SMILES
        if RDKIT_AVAILABLE:
            mol = Chem.MolFromSmiles(identifier)
            if mol is not None:
                return "smiles"
    
    # Check for simple element combinations that might be SMILES
    if re.match(r'^[A-Za-z0-9\(\)\[\]=#@\\/+-]+$', identifier) and len(identifier) > 2:
        if RDKIT_AVAILABLE:
            mol = Chem.MolFromSmiles(identifier)
            if mol is not None:
                return "smiles"
    
    # Check if it looks like IUPAC (contains numbers, dashes, commas in specific patterns)
    iupac_patterns = [
        r'\d+-',           # Numbers followed by dash
        r'-\d+,\d+-',      # Position indicators like -3,4-
        r'yl$',            # Ends with -yl
        r'oic acid$',      # Carboxylic acid ending
        r'one$',           # Ketone ending
        r'ol$',            # Alcohol ending
        r'ene$',           # Alkene ending
        r'ane$',           # Alkane ending
    ]
    if any(re.search(pattern, identifier.lower()) for pattern in iupac_patterns):
        return "iupac"
    
    # Default to common name
    return "name"


def smiles_to_3d_structure(
    smiles: str, 
    optimize: bool = True,
    num_conformers: int = 1,
    max_attempts: int = 5
) -> Optional[Dict[str, Any]]:
    """
    Generate 3D molecular structure from SMILES string using RDKit.
    
    Args:
        smiles: SMILES string
        optimize: Whether to optimize geometry with force field
        num_conformers: Number of conformers to generate (returns lowest energy)
        max_attempts: Maximum embedding attempts
        
    Returns:
        Dict with atoms, coords, and molecule data, or None if failed
    """
    if not RDKIT_AVAILABLE:
        return None
    
    try:
        # Parse SMILES
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            logger.error(f"Failed to parse SMILES: {smiles}")
            return None
        
        # Add hydrogens explicitly
        mol = Chem.AddHs(mol)
        
        # Generate 3D conformers using ETKDG
        # Use try/except for version compatibility
        try:
            params = AllChem.ETKDGv3()
            params.randomSeed = 42
            params.useRandomCoords = True
        except:
            # Fallback for older RDKit versions
            params = AllChem.ETKDG()
            params.randomSeed = 42
        
        # Embed molecule
        conf_ids = AllChem.EmbedMultipleConfs(mol, numConfs=num_conformers, params=params)
        
        if len(conf_ids) == 0:
            # Fallback: try with random coordinates and single conf
            result_code = AllChem.EmbedMolecule(mol, randomSeed=42, useRandomCoords=True, maxAttempts=max_attempts)
            if result_code == 0:
                conf_ids = [0]
        
        if len(conf_ids) == 0:
            logger.warning(f"Failed to generate 3D conformer for: {smiles}. Falling back to 2D coordinates.")
            # Fallback to 2D coordinates
            # Compute2DCoords returns 0 on success
            if AllChem.Compute2DCoords(mol) == 0:
                # Force Z coordinates to 0.0 for 2D
                # Note: Compute2DCoords puts them on Z=0 anyway
                conf_ids = [0]
            else:
                logger.error(f"Failed to generate 2D coordinates for: {smiles}")
                return None
        
        # Optimize geometry
        if optimize:
            # Try MMFF94 first, fall back to UFF
            try:
                results = AllChem.MMFFOptimizeMoleculeConfs(mol, maxIters=500)
                if results and all(r[0] == 0 for r in results):  # All converged
                    pass  # MMFF succeeded
                else:
                    AllChem.UFFOptimizeMoleculeConfs(mol, maxIters=500)
            except:
                try:
                    AllChem.UFFOptimizeMoleculeConfs(mol, maxIters=500)
                except:
                    pass  # Use unoptimized coords
        
        # Get lowest energy conformer
        best_conf_id = conf_ids[0]
        if len(conf_ids) > 1 and optimize:
            try:
                energies = []
                for conf_id in conf_ids:
                    ff = AllChem.MMFFGetMoleculeForceField(mol, AllChem.MMFFGetMoleculeProperties(mol), confId=conf_id)
                    if ff:
                        energies.append((conf_id, ff.CalcEnergy()))
                    else:
                        ff = AllChem.UFFGetMoleculeForceField(mol, confId=conf_id)
                        if ff:
                            energies.append((conf_id, ff.CalcEnergy()))
                if energies:
                    best_conf_id = min(energies, key=lambda x: x[1])[0]
            except:
                pass
        
        # Extract coordinates
        conformer = mol.GetConformer(best_conf_id)
        atoms = []
        coords = []
        
        for i, atom in enumerate(mol.GetAtoms()):
            atoms.append(atom.GetSymbol())
            pos = conformer.GetAtomPosition(i)
            coords.append([pos.x, pos.y, pos.z])
        
        # Calculate molecular properties
        mol_no_h = Chem.RemoveHs(mol)
        mol_weight = Descriptors.MolWt(mol)
        formula = Chem.rdMolDescriptors.CalcMolFormula(mol)
        
        return {
            "success": True,
            "source": "rdkit",
            "smiles": smiles,
            "canonical_smiles": Chem.MolToSmiles(mol_no_h),
            "formula": formula,
            "molecular_weight": mol_weight,
            "n_atoms": len(atoms),
            "atoms": atoms,
            "coords": coords,
            "optimized": optimize,
        }
        
    except Exception as e:
        logger.error(f"Error generating 3D structure from SMILES '{smiles}': {e}")
        return None


def fetch_from_pubchem(
    identifier: str, 
    by: str = "name"
) -> Optional[Dict[str, Any]]:
    """
    Fetch molecule data from PubChem API.
    
    Args:
        identifier: Molecule identifier (name, CID, or SMILES)
        by: Search type - "name", "cid", or "smiles"
        
    Returns:
        Dict with SMILES and other data, or None if not found
    """
    if not REQUESTS_AVAILABLE:
        return None
    
    base_url = "https://pubchem.ncbi.nlm.nih.gov/rest/pug"
    
    try:
        # First, get the CID
        if by == "cid":
            cid = identifier
        else:
            # Search by name or SMILES
            if by == "smiles":
                # URL encode the SMILES
                import urllib.parse
                encoded = urllib.parse.quote(identifier, safe='')
                search_url = f"{base_url}/compound/smiles/{encoded}/cids/JSON"
            else:
                search_url = f"{base_url}/compound/name/{identifier}/cids/JSON"
            
            response = requests.get(search_url, timeout=10)
            if response.status_code != 200:
                return None
            
            data = response.json()
            cids = data.get("IdentifierList", {}).get("CID", [])
            if not cids:
                return None
            cid = cids[0]
        
        # Get compound properties including SMILES
        # Request both Canonical and Isomeric SMILES for robustness
        props_url = f"{base_url}/compound/cid/{cid}/property/CanonicalSMILES,IsomericSMILES,MolecularFormula,MolecularWeight,IUPACName/JSON"
        response = requests.get(props_url, timeout=10)
        if response.status_code != 200:
            return None
        
        data = response.json()
        props = data.get("PropertyTable", {}).get("Properties", [])
        if not props:
            return None
        
        prop = props[0]
        
        # Try to get 3D coordinates directly from PubChem
        coords_3d = None
        try:
            sdf_url = f"{base_url}/compound/cid/{cid}/record/SDF?record_type=3d"
            sdf_response = requests.get(sdf_url, timeout=15)
            if sdf_response.status_code == 200:
                # Parse SDF to extract coordinates
                coords_3d = parse_sdf_coordinates(sdf_response.text)
        except:
            pass
        
        return {
            "success": True,
            "source": "pubchem",
            "cid": cid,
            "smiles": prop.get("IsomericSMILES") or prop.get("CanonicalSMILES") or prop.get("ConnectivitySMILES"),
            "formula": prop.get("MolecularFormula"),
            "molecular_weight": prop.get("MolecularWeight"),
            "iupac_name": prop.get("IUPACName"),
            "coords_3d": coords_3d,
        }
        
    except requests.exceptions.Timeout:
        logger.warning(f"PubChem API timeout for: {identifier}")
        return None
    except Exception as e:
        logger.error(f"Error fetching from PubChem: {e}")
        return None


def parse_sdf_coordinates(sdf_text: str) -> Optional[Dict[str, Any]]:
    """
    Parse atom coordinates from SDF file content.
    
    Args:
        sdf_text: SDF file content as string
        
    Returns:
        Dict with atoms and coords, or None if parsing failed
    """
    try:
        lines = sdf_text.strip().split('\n')
        
        # Find the counts line (4th line typically)
        counts_line = None
        for i, line in enumerate(lines):
            parts = line.split()
            if len(parts) >= 2:
                try:
                    n_atoms = int(parts[0])
                    n_bonds = int(parts[1])
                    if n_atoms > 0 and n_atoms < 1000:
                        counts_line = i
                        break
                except ValueError:
                    continue
        
        if counts_line is None:
            return None
        
        n_atoms = int(lines[counts_line].split()[0])
        
        atoms = []
        coords = []
        
        # Read atom block
        for i in range(counts_line + 1, counts_line + 1 + n_atoms):
            if i >= len(lines):
                break
            parts = lines[i].split()
            if len(parts) >= 4:
                x, y, z = float(parts[0]), float(parts[1]), float(parts[2])
                element = parts[3]
                atoms.append(element)
                coords.append([x, y, z])
        
        if len(atoms) == n_atoms:
            return {"atoms": atoms, "coords": coords}
        
        return None
        
    except Exception as e:
        logger.error(f"Error parsing SDF: {e}")
        return None


def resolve_to_smiles(identifier: str, input_type: str = "auto") -> Optional[str]:
    """
    Resolve any molecular identifier to a SMILES string.
    
    Args:
        identifier: The molecular identifier
        input_type: Type of identifier ("auto", "name", "smiles", "iupac", "cid", "inchi")
        
    Returns:
        SMILES string or None if resolution failed
    """
    if input_type == "auto":
        input_type = detect_input_type(identifier)
    
    # If already SMILES, validate and return
    if input_type == "smiles":
        if RDKIT_AVAILABLE:
            mol = Chem.MolFromSmiles(identifier)
            if mol:
                return identifier
        return identifier
    
    # If InChI, convert to SMILES
    if input_type == "inchi" and RDKIT_AVAILABLE:
        mol = Chem.MolFromInchi(identifier)
        if mol:
            return Chem.MolToSmiles(mol)
    
    # Check local aliases first
    identifier_lower = identifier.lower() if isinstance(identifier, str) else str(identifier).lower()
    alias_target = None
    
    if identifier in MOLECULE_ALIASES:
        alias_target = MOLECULE_ALIASES[identifier]
    elif identifier_lower in MOLECULE_ALIASES:
        alias_target = MOLECULE_ALIASES[identifier_lower]
        
    if alias_target:
        # If alias leads to another name/SMILES, check what it is
        # If it's a SMILES, return it
        if RDKIT_AVAILABLE and Chem.MolFromSmiles(alias_target):
             return alias_target
        # Otherwise, treat it as a refined identifier and continue (recurse effectively)
        identifier = alias_target
        # Re-check input type for the new identifier
        if input_type == "auto":
             new_type = detect_input_type(identifier)
             if new_type != "name": # If it resolved to SMILES/CID
                 if new_type == "smiles": return identifier
                 # For CID/InChI, flow through
                 input_type = new_type
    
    # Check local database
    
    # Check local database
    if identifier.upper() in MOLECULE_DATABASE:
        # Local DB has coords but not SMILES, so we can't easily get SMILES
        # Return None to indicate we should use local DB directly
        return None
    
    # Try PubChem
    if input_type == "cid":
        result = fetch_from_pubchem(identifier, by="cid")
    else:
        result = fetch_from_pubchem(identifier, by="name")
    
    if result and result.get("smiles"):
        return result["smiles"]
    
    # Try OPSIN for IUPAC names (via web service)
    if input_type == "iupac":
        try:
            import urllib.parse
            encoded_name = urllib.parse.quote(identifier)
            opsin_url = f"https://opsin.ch.cam.ac.uk/opsin/{encoded_name}.smi"
            response = requests.get(opsin_url, timeout=10)
            if response.status_code == 200:
                smiles = response.text.strip()
                if smiles and not smiles.startswith("<!"):
                    return smiles
        except:
            pass
    
    return None


def generate_molecule_universal(
    identifier: str,
    input_type: str = "auto",
    optimize: bool = True,
    allow_external: bool = True
) -> Dict[str, Any]:
    """
    Generate 3D molecular structure from any identifier.
    
    This is the main entry point for universal molecule generation.
    It handles:
    - Common names (aspirin, PTCDA, benzene, etc.)
    - SMILES strings
    - IUPAC names
    - PubChem CIDs
    - InChI strings
    
    Args:
        identifier: Molecule identifier (name, SMILES, IUPAC, CID, or InChI)
        input_type: Type of identifier ("auto", "name", "smiles", "iupac", "cid", "inchi")
        optimize: Whether to optimize geometry with force field
        allow_external: Whether to query external APIs (PubChem, OPSIN)
        
    Returns:
        Dict with molecule structure data and metadata
    """
    identifier = identifier.strip()
    
    # Detect input type if auto
    if input_type == "auto":
        input_type = detect_input_type(identifier)
    
    # Check local database first (for molecules we have accurate hand-crafted coords)
    if identifier.upper() in MOLECULE_DATABASE or identifier in MOLECULE_DATABASE:
        key = identifier.upper() if identifier.upper() in MOLECULE_DATABASE else identifier
        info = MOLECULE_DATABASE[key]
        
        # Build pymatgen molecule if available
        molecule_data = None
        if PYMATGEN_AVAILABLE:
            mol = Molecule(info["atoms"], info["coords"])
            molecule_data = {
                "species": [str(s) for s in mol.species],
                "coords": [list(c) for c in mol.cart_coords]
            }
        
        return {
            "success": True,
            "source": "local_database",
            "identifier": identifier,
            "formula": key,
            "n_atoms": len(info["atoms"]),
            "atoms": info["atoms"],
            "coords": info["coords"],
            "spin_multiplicity": 2 * info.get("spin", 0) + 1,
            "is_radical": info.get("radical", False),
            "symmetry": info.get("symmetry", ""),
            "is_linear": info.get("linear", False),
            "molecule": molecule_data,
        }
    
    # Check unified SQLite molecule database (ChEMBL, PubChem, ZINC, etc.)
    db_result = db_lookup_molecule(identifier)
    if db_result:
        smiles = db_result.get("smiles")
        if smiles:
            result = smiles_to_3d_structure(smiles, optimize=optimize)
            if result:
                result["identifier"] = identifier
                result["source"] = f"database_{db_result.get('source', 'local')}"
                result["formula"] = db_result.get("formula") or result.get("formula")
                result["category"] = db_result.get("category")
                result["pubchem_cid"] = db_result.get("pubchem_cid")
                result["chembl_id"] = db_result.get("chembl_id")
                
                # Build pymatgen molecule
                if PYMATGEN_AVAILABLE and Molecule is not None:
                    mol = Molecule(result["atoms"], result["coords"])
                    result["molecule"] = {
                        "species": [str(s) for s in mol.species],
                        "coords": [list(c) for c in mol.cart_coords]
                    }
                
                return result
    
    # Alias / Redirect handling
    # This handles "C60" -> "Buckminsterfullerene", "Vitamin B12" -> "Cyanocobalamin"
    if identifier in MOLECULE_ALIASES or identifier.lower() in MOLECULE_ALIASES:
        alias_key = identifier if identifier in MOLECULE_ALIASES else identifier.lower()
        target = MOLECULE_ALIASES[alias_key]
        
        # If target looks like SMILES (rudimentary check), try generating
        is_smiles = False
        if RDKIT_AVAILABLE and Chem.MolFromSmiles(target):
             is_smiles = True
        
        if is_smiles:
            result = smiles_to_3d_structure(target, optimize=optimize)
            if result:
                result["identifier"] = identifier
                result["alias_matched"] = alias_key
                # Build pymatgen molecule
                if PYMATGEN_AVAILABLE and Molecule is not None:
                    mol = Molecule(result["atoms"], result["coords"])
                    result["molecule"] = {
                        "species": [str(s) for s in mol.species],
                        "coords": [list(c) for c in mol.cart_coords]
                    }
                return result
        else:
            # Target is a Name (Redirect)
            identifier = target
            # Re-detect type for the new name
            input_type = detect_input_type(identifier)
    
    # If input is SMILES, generate directly
    if input_type == "smiles":
        result = smiles_to_3d_structure(identifier, optimize=optimize)
        if result:
            result["identifier"] = identifier
            
            # Build pymatgen molecule
            if PYMATGEN_AVAILABLE:
                mol = Molecule(result["atoms"], result["coords"])
                result["molecule"] = {
                    "species": [str(s) for s in mol.species],
                    "coords": [list(c) for c in mol.cart_coords]
                }
            
            return result
        return {
            "success": False,
            "error": {
                "code": "SMILES_PARSE_ERROR",
                "message": f"Failed to parse SMILES: {identifier}",
                "suggestion": "Check SMILES syntax or try a different representation"
            }
        }
    
    # If InChI, convert and generate
    if input_type == "inchi" and RDKIT_AVAILABLE:
        mol = Chem.MolFromInchi(identifier)
        if mol:
            smiles = Chem.MolToSmiles(mol)
            result = smiles_to_3d_structure(smiles, optimize=optimize)
            if result:
                result["identifier"] = identifier
                result["inchi"] = identifier
                
                # Build pymatgen molecule
                if PYMATGEN_AVAILABLE:
                    mol = Molecule(result["atoms"], result["coords"])
                    result["molecule"] = {
                        "species": [str(s) for s in mol.species],
                        "coords": [list(c) for c in mol.cart_coords]
                    }
                
                return result
    
    # External resolution (if allowed)
    if not allow_external:
        return {
            "success": False,
            "error": {
                "code": "NOT_FOUND_OFFLINE",
                "message": f"Molecule '{identifier}' not found in local database",
                "suggestion": "Enable external APIs or provide SMILES directly"
            }
        }
    
    # Try PubChem
    pubchem_result = None
    if input_type == "cid":
        pubchem_result = fetch_from_pubchem(identifier, by="cid")
    else:
        pubchem_result = fetch_from_pubchem(identifier, by="name")
    
    if pubchem_result and pubchem_result.get("success"):
        # Check if PubChem returned 3D coords
        if pubchem_result.get("coords_3d"):
            coord_data = pubchem_result["coords_3d"]
            
            # Build pymatgen molecule
            molecule_data = None
            if PYMATGEN_AVAILABLE:
                mol = Molecule(coord_data["atoms"], coord_data["coords"])
                molecule_data = {
                    "species": [str(s) for s in mol.species],
                    "coords": [list(c) for c in mol.cart_coords]
                }
            
            return {
                "success": True,
                "source": "pubchem_3d",
                "identifier": identifier,
                "pubchem_cid": pubchem_result.get("cid"),
                "smiles": pubchem_result.get("smiles"),
                "formula": pubchem_result.get("formula"),
                "molecular_weight": pubchem_result.get("molecular_weight"),
                "iupac_name": pubchem_result.get("iupac_name"),
                "n_atoms": len(coord_data["atoms"]),
                "atoms": coord_data["atoms"],
                "coords": coord_data["coords"],
                "molecule": molecule_data,
            }
        
        
        # Generate 3D from SMILES
        smiles = pubchem_result.get("smiles")
        if smiles:
            result = smiles_to_3d_structure(smiles, optimize=optimize)
            if result:
                result["identifier"] = identifier
                result["pubchem_cid"] = pubchem_result.get("cid")
                result["iupac_name"] = pubchem_result.get("iupac_name")
                result["source"] = "pubchem_smiles_rdkit"
                
                # Cache to local database for future offline use
                db = get_molecule_database()
                db.add_molecule(
                    smiles=smiles,
                    name=identifier,
                    source="pubchem_cache",
                    formula=pubchem_result.get("formula"),
                    pubchem_cid=pubchem_result.get("cid"),
                )
                
                # Build pymatgen molecule
                if PYMATGEN_AVAILABLE:
                    mol = Molecule(result["atoms"], result["coords"])
                    result["molecule"] = {
                        "species": [str(s) for s in mol.species],
                        "coords": [list(c) for c in mol.cart_coords]
                    }
                
                return result
    
    # Try OPSIN for IUPAC-like names
    if input_type in ["iupac", "name"] and REQUESTS_AVAILABLE and requests is not None:
        import urllib.parse
        encoded_name = urllib.parse.quote(identifier)
        opsin_url = f"https://opsin.ch.cam.ac.uk/opsin/{encoded_name}.smi"
        response = requests.get(opsin_url, timeout=10)
        
        if response.status_code == 200:
            smiles = response.text.strip()
            # Validate response is SMILES, not HTML error
            if smiles and not smiles.startswith("<!") and not smiles.startswith("<"):
                result = smiles_to_3d_structure(smiles, optimize=optimize)
                if result:
                    result["identifier"] = identifier
                    result["source"] = "opsin_rdkit"
                    
                    # Cache to local database
                    db = get_molecule_database()
                    db.add_molecule(
                        smiles=smiles,
                        name=identifier,
                        source="opsin_cache",
                    )
                    
                    # Build pymatgen molecule
                    if PYMATGEN_AVAILABLE:
                        mol = Molecule(result["atoms"], result["coords"])
                        result["molecule"] = {
                            "species": [str(s) for s in mol.species],
                            "coords": [list(c) for c in mol.cart_coords]
                        }
                    
                    return result
    
    # All methods failed
    return {
        "success": False,
        "error": {
            "code": "MOLECULE_NOT_FOUND",
            "message": f"Could not resolve molecule: {identifier}",
            "detected_type": input_type,
            "suggestions": [
                "Check spelling of the molecule name",
                "Try the IUPAC systematic name",
                "Provide the SMILES string directly",
                f"Search PubChem for '{identifier}' to find the correct identifier"
            ]
        }
    }


def get_supported_molecules() -> Dict[str, Any]:
    """
    Get information about supported molecule generation capabilities.
    
    Returns:
        Dict with local database contents and supported input types
    """
    return {
        "local_database": list(MOLECULE_DATABASE.keys()),
        "local_aliases": list(MOLECULE_ALIASES.keys()),
        "supported_input_types": ["name", "smiles", "iupac", "cid", "inchi"],
        "external_sources": ["pubchem", "opsin"],
        "rdkit_available": RDKIT_AVAILABLE,
        "requests_available": REQUESTS_AVAILABLE,
        "capabilities": {
            "arbitrary_smiles": RDKIT_AVAILABLE,
            "name_resolution": REQUESTS_AVAILABLE,
            "iupac_resolution": REQUESTS_AVAILABLE,
            "3d_optimization": RDKIT_AVAILABLE,
        }
    }
