"""
molecule/conformers.py - Molecular Conformer Generation

Provides conformer generation and analysis:
- Systematic conformer search
- MMFF94/UFF force field optimization
- Rotamer enumeration
- Tautomer/isomer generation
- Protonation state prediction
"""

from typing import Dict, Any, List, Optional, Tuple, Union
import numpy as np


def generate_conformers(
    smiles: str,
    n_conformers: int = 10,
    force_field: str = "MMFF94",
    energy_window_kcal: float = 10.0,
    rms_threshold: float = 0.5,
    optimize: bool = True,
    random_seed: Optional[int] = None
) -> Dict[str, Any]:
    """
    Generate molecular conformers using RDKit.

    Uses systematic torsion driving and force field optimization
    to generate diverse low-energy conformers.

    Args:
        smiles: SMILES string of molecule
        n_conformers: Maximum number of conformers to generate
        force_field: "MMFF94", "MMFF94s", "UFF"
        energy_window_kcal: Keep conformers within this energy of minimum
        rms_threshold: RMS threshold for conformer clustering (Angstroms)
        optimize: Whether to optimize geometries
        random_seed: Random seed for reproducibility

    Returns:
        List of conformers with energies and coordinates
    """
    from rdkit import Chem
    from rdkit.Chem import AllChem, Descriptors
    from rdkit.Chem.rdMolDescriptors import CalcNumRotatableBonds

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return {
            "success": False,
            "error": {
                "code": "INVALID_SMILES",
                "message": f"Could not parse SMILES: {smiles}"
            }
        }

    # Add hydrogens
    mol = Chem.AddHs(mol)

    # Get molecular properties
    n_rotatable = CalcNumRotatableBonds(mol)
    n_atoms = mol.GetNumAtoms()
    molecular_weight = Descriptors.ExactMolWt(mol)

    # Determine number of conformers to generate (more for flexible molecules)
    n_to_generate = max(n_conformers * 3, n_rotatable * 5, 50)
    n_to_generate = min(n_to_generate, 500)  # Cap at reasonable number

    # Set random seed
    params = AllChem.ETKDGv3()
    if random_seed is not None:
        params.randomSeed = random_seed
    params.useSmallRingTorsions = True
    params.useMacrocycleTorsions = True
    params.enforceChirality = True

    # Generate conformers
    conf_ids = AllChem.EmbedMultipleConfs(mol, numConfs=n_to_generate, params=params)

    if len(conf_ids) == 0:
        # Try without strict distance geometry
        params.useRandomCoords = True
        conf_ids = AllChem.EmbedMultipleConfs(mol, numConfs=n_to_generate, params=params)

    if len(conf_ids) == 0:
        return {
            "success": False,
            "error": {
                "code": "EMBEDDING_FAILED",
                "message": "Could not generate 3D conformers"
            }
        }

    # Set up force field
    ff_class = None
    if force_field.upper() in ["MMFF94", "MMFF94S"]:
        ff_getter = lambda mol, conf_id: AllChem.MMFFGetMoleculeForceField(
            mol, AllChem.MMFFGetMoleculeProperties(mol, mmffVariant=force_field),
            confId=conf_id
        )
    else:  # UFF
        ff_getter = lambda mol, conf_id: AllChem.UFFGetMoleculeForceField(mol, confId=conf_id)

    # Optimize and calculate energies
    conformer_data = []

    for conf_id in conf_ids:
        ff = ff_getter(mol, conf_id)
        if ff is None:
            continue

        if optimize:
            # Optimize geometry
            converged = ff.Minimize(maxIts=500)

        energy = ff.CalcEnergy()

        # Get coordinates
        conf = mol.GetConformer(conf_id)
        coords = []
        symbols = []
        for i in range(mol.GetNumAtoms()):
            pos = conf.GetAtomPosition(i)
            coords.append([pos.x, pos.y, pos.z])
            symbols.append(mol.GetAtomWithIdx(i).GetSymbol())

        conformer_data.append({
            "conf_id": int(conf_id),
            "energy_kcal_mol": float(energy),
            "coords": coords,
            "symbols": symbols,
            "converged": converged if optimize else None
        })

    # Sort by energy
    conformer_data.sort(key=lambda x: x["energy_kcal_mol"])

    # Filter by energy window
    if conformer_data:
        min_energy = conformer_data[0]["energy_kcal_mol"]
        conformer_data = [c for c in conformer_data
                        if c["energy_kcal_mol"] - min_energy <= energy_window_kcal]

    # Cluster by RMS to remove duplicates
    unique_conformers = []
    for conf in conformer_data:
        is_unique = True
        for unique in unique_conformers:
            # Calculate RMS between conformers
            coords1 = np.array(conf["coords"])
            coords2 = np.array(unique["coords"])
            rms = np.sqrt(np.mean((coords1 - coords2) ** 2))
            if rms < rms_threshold:
                is_unique = False
                break
        if is_unique:
            unique_conformers.append(conf)
        if len(unique_conformers) >= n_conformers:
            break

    # Calculate relative energies
    if unique_conformers:
        min_e = unique_conformers[0]["energy_kcal_mol"]
        for conf in unique_conformers:
            conf["relative_energy_kcal_mol"] = conf["energy_kcal_mol"] - min_e
            # Boltzmann weight at 298 K
            conf["boltzmann_weight"] = np.exp(-conf["relative_energy_kcal_mol"] / (0.001987 * 298))

        # Normalize Boltzmann weights
        total_weight = sum(c["boltzmann_weight"] for c in unique_conformers)
        for conf in unique_conformers:
            conf["boltzmann_population"] = conf["boltzmann_weight"] / total_weight

    return {
        "success": True,
        "smiles": smiles,
        "n_atoms": n_atoms,
        "n_rotatable_bonds": n_rotatable,
        "molecular_weight": molecular_weight,
        "force_field": force_field,
        "n_conformers_generated": len(conf_ids),
        "n_conformers_unique": len(unique_conformers),
        "energy_window_kcal": energy_window_kcal,
        "rms_threshold_A": rms_threshold,
        "conformers": unique_conformers
    }


def generate_tautomers(
    smiles: str,
    max_tautomers: int = 20,
    include_scores: bool = True
) -> Dict[str, Any]:
    """
    Generate tautomers of a molecule.

    Uses RDKit's tautomer enumeration to find all reasonable tautomeric forms.

    Args:
        smiles: SMILES string of molecule
        max_tautomers: Maximum number of tautomers to generate
        include_scores: Include tautomer stability scores

    Returns:
        List of tautomers with SMILES and scores
    """
    from rdkit import Chem
    from rdkit.Chem.MolStandardize import rdMolStandardize

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return {
            "success": False,
            "error": {"code": "INVALID_SMILES", "message": f"Could not parse: {smiles}"}
        }

    # Get tautomer enumerator
    enumerator = rdMolStandardize.TautomerEnumerator()

    # Enumerate tautomers
    tautomers = enumerator.Enumerate(mol)
    tautomer_list = list(tautomers)[:max_tautomers]

    # Get canonical tautomer
    canonical = enumerator.Canonicalize(mol)
    canonical_smiles = Chem.MolToSmiles(canonical)

    results = []
    for i, taut in enumerate(tautomer_list):
        taut_smiles = Chem.MolToSmiles(taut)
        result = {
            "id": i,
            "smiles": taut_smiles,
            "is_canonical": taut_smiles == canonical_smiles
        }

        if include_scores:
            # RDKit tautomer score (higher = more stable)
            score = enumerator.ScoreTautomer(taut)
            result["score"] = float(score)

        results.append(result)

    # Sort by score if available
    if include_scores:
        results.sort(key=lambda x: x.get("score", 0), reverse=True)

    return {
        "success": True,
        "input_smiles": smiles,
        "canonical_tautomer": canonical_smiles,
        "n_tautomers": len(results),
        "tautomers": results
    }


def enumerate_stereoisomers(
    smiles: str,
    include_unassigned: bool = True,
    max_isomers: int = 100
) -> Dict[str, Any]:
    """
    Enumerate all stereoisomers of a molecule.

    Generates all possible E/Z and R/S stereoisomers.

    Args:
        smiles: SMILES string
        include_unassigned: Include stereoisomers with unassigned centers
        max_isomers: Maximum number to generate

    Returns:
        List of stereoisomers
    """
    from rdkit import Chem
    from rdkit.Chem import AllChem
    from rdkit.Chem.EnumerateStereoisomers import EnumerateStereoisomers, StereoEnumerationOptions

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return {
            "success": False,
            "error": {"code": "INVALID_SMILES", "message": f"Could not parse: {smiles}"}
        }

    # Count stereocenters
    chiral_centers = Chem.FindMolChiralCenters(mol, includeUnassigned=True)
    n_chiral = len(chiral_centers)

    # Count double bonds that can be E/Z
    n_ez = 0
    for bond in mol.GetBonds():
        if bond.GetBondType() == Chem.BondType.DOUBLE:
            if bond.GetStereo() != Chem.BondStereo.STEREONONE:
                n_ez += 1

    # Set enumeration options
    opts = StereoEnumerationOptions(
        tryEmbedding=True,
        onlyUnassigned=not include_unassigned,
        maxIsomers=max_isomers
    )

    # Enumerate
    isomers = list(EnumerateStereoisomers(mol, options=opts))

    results = []
    for i, iso in enumerate(isomers):
        iso_smiles = Chem.MolToSmiles(iso, isomericSmiles=True)
        centers = Chem.FindMolChiralCenters(iso, includeUnassigned=False)

        results.append({
            "id": i,
            "smiles": iso_smiles,
            "n_defined_centers": len(centers),
            "stereocenters": [{"atom_idx": c[0], "config": c[1]} for c in centers]
        })

    return {
        "success": True,
        "input_smiles": smiles,
        "n_chiral_centers": n_chiral,
        "n_e_z_bonds": n_ez,
        "max_possible_isomers": 2 ** (n_chiral + n_ez),
        "n_isomers_found": len(results),
        "stereoisomers": results
    }


def predict_protonation_states(
    smiles: str,
    pH: float = 7.4,
    pH_range: Optional[Tuple[float, float]] = None
) -> Dict[str, Any]:
    """
    Predict protonation states at given pH.

    Uses pKa predictions to determine likely protonation states.

    Args:
        smiles: SMILES string
        pH: Target pH for protonation state
        pH_range: If provided, generate states for pH range

    Returns:
        Protonation states with populations
    """
    from rdkit import Chem
    from rdkit.Chem import Descriptors, AllChem

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return {
            "success": False,
            "error": {"code": "INVALID_SMILES", "message": f"Could not parse: {smiles}"}
        }

    # Identify ionizable groups
    ionizable_groups = []

    # Common ionizable patterns with approximate pKa values
    patterns = {
        "carboxylic_acid": ("[CX3](=O)[OX1H1]", 4.0, "acid"),
        "sulfonic_acid": ("[SX4](=O)(=O)[OX1H1]", -1.0, "acid"),
        "phosphoric_acid": ("[PX4](=O)([OX1H1])([OX1H1])[OX1H1]", 2.0, "acid"),
        "phenol": ("[OX1H1][c]", 10.0, "acid"),
        "primary_amine": ("[NX3H2;!$(NC=O)]", 10.0, "base"),
        "secondary_amine": ("[NX3H1;!$(NC=O)]", 10.0, "base"),
        "tertiary_amine": ("[NX3H0;!$(NC=O)]", 9.0, "base"),
        "guanidine": ("[NX3][CX3](=[NX2+])[NX3]", 13.0, "base"),
        "imidazole": ("[nX2]1[cX3][nX3H1][cX3][cX3]1", 6.0, "base"),
        "pyridine": ("[nX2]1[cX3][cX3][cX3][cX3][cX3]1", 5.2, "base"),
    }

    for name, (pattern, pka, type_) in patterns.items():
        patt = Chem.MolFromSmarts(pattern)
        if patt and mol.HasSubstructMatch(patt):
            matches = mol.GetSubstructMatches(patt)
            for match in matches:
                ionizable_groups.append({
                    "name": name,
                    "atom_indices": list(match),
                    "pKa": pka,
                    "type": type_
                })

    # Calculate protonation probabilities at target pH
    for group in ionizable_groups:
        pka = group["pKa"]
        if group["type"] == "acid":
            # For acids: fraction deprotonated
            group["fraction_ionized"] = 1 / (1 + 10 ** (pka - pH))
            group["state_at_pH"] = "deprotonated" if group["fraction_ionized"] > 0.5 else "protonated"
        else:
            # For bases: fraction protonated
            group["fraction_ionized"] = 1 / (1 + 10 ** (pH - pka))
            group["state_at_pH"] = "protonated" if group["fraction_ionized"] > 0.5 else "neutral"

    # Calculate overall charge at pH
    total_charge = 0
    for group in ionizable_groups:
        if group["type"] == "acid":
            total_charge -= group["fraction_ionized"]
        else:
            total_charge += group["fraction_ionized"]

    # Generate major protonation state SMILES (simplified)
    # In reality, would modify the molecule

    return {
        "success": True,
        "input_smiles": smiles,
        "pH": pH,
        "n_ionizable_groups": len(ionizable_groups),
        "ionizable_groups": ionizable_groups,
        "predicted_charge": round(total_charge, 2),
        "isoelectric_point": None,  # Would require more complex calculation
        "note": "pKa values are approximate; use specialized software for accuracy"
    }


def generate_rotamers(
    smiles: str,
    rotatable_bond_idx: Optional[int] = None,
    n_rotamers: int = 12,
    optimize_each: bool = True
) -> Dict[str, Any]:
    """
    Generate rotamers by systematic torsion scanning.

    Rotates around specified or all rotatable bonds.

    Args:
        smiles: SMILES string
        rotatable_bond_idx: Specific bond index to rotate (None = all)
        n_rotamers: Number of rotamers per bond (360/n spacing)
        optimize_each: Optimize geometry after rotation

    Returns:
        Rotamer structures with energies
    """
    from rdkit import Chem
    from rdkit.Chem import AllChem
    from rdkit.Chem.Lipinski import RotatableBondSmarts

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return {
            "success": False,
            "error": {"code": "INVALID_SMILES", "message": f"Could not parse: {smiles}"}
        }

    mol = Chem.AddHs(mol)

    # Generate initial 3D structure
    AllChem.EmbedMolecule(mol, AllChem.ETKDGv3())
    AllChem.MMFFOptimizeMolecule(mol)

    # Find rotatable bonds
    rot_bond_pattern = Chem.MolFromSmarts(RotatableBondSmarts)
    rotatable_bonds = mol.GetSubstructMatches(rot_bond_pattern)

    if not rotatable_bonds:
        return {
            "success": True,
            "smiles": smiles,
            "n_rotatable_bonds": 0,
            "message": "No rotatable bonds found",
            "rotamers": []
        }

    # Select bonds to rotate
    if rotatable_bond_idx is not None:
        if rotatable_bond_idx >= len(rotatable_bonds):
            return {
                "success": False,
                "error": {
                    "code": "INVALID_BOND_INDEX",
                    "message": f"Bond index {rotatable_bond_idx} out of range (0-{len(rotatable_bonds)-1})"
                }
            }
        bonds_to_scan = [rotatable_bonds[rotatable_bond_idx]]
    else:
        bonds_to_scan = rotatable_bonds

    angle_step = 360.0 / n_rotamers
    rotamers = []

    for bond_atoms in bonds_to_scan:
        atom1, atom2 = bond_atoms

        # Get atoms defining the torsion (need 4 atoms)
        atom0 = None
        atom3 = None

        for neighbor in mol.GetAtomWithIdx(atom1).GetNeighbors():
            if neighbor.GetIdx() != atom2:
                atom0 = neighbor.GetIdx()
                break

        for neighbor in mol.GetAtomWithIdx(atom2).GetNeighbors():
            if neighbor.GetIdx() != atom1:
                atom3 = neighbor.GetIdx()
                break

        if atom0 is None or atom3 is None:
            continue

        torsion = (atom0, atom1, atom2, atom3)

        # Scan torsion
        for i in range(n_rotamers):
            angle = i * angle_step

            # Set torsion angle
            mol_copy = Chem.Mol(mol)
            AllChem.EmbedMolecule(mol_copy, AllChem.ETKDGv3())
            conf = mol_copy.GetConformer()

            Chem.rdMolTransforms.SetDihedralDeg(conf, *torsion, angle)

            if optimize_each:
                ff = AllChem.MMFFGetMoleculeForceField(
                    mol_copy,
                    AllChem.MMFFGetMoleculeProperties(mol_copy)
                )
                if ff:
                    ff.Minimize(maxIts=100)
                    energy = ff.CalcEnergy()
                else:
                    energy = None
            else:
                energy = None

            # Get coordinates
            coords = []
            symbols = []
            for j in range(mol_copy.GetNumAtoms()):
                pos = conf.GetAtomPosition(j)
                coords.append([pos.x, pos.y, pos.z])
                symbols.append(mol_copy.GetAtomWithIdx(j).GetSymbol())

            rotamers.append({
                "bond": [atom1, atom2],
                "torsion_atoms": list(torsion),
                "angle_deg": angle,
                "energy_kcal_mol": float(energy) if energy else None,
                "coords": coords,
                "symbols": symbols
            })

    # Sort by energy if available
    rotamers_with_energy = [r for r in rotamers if r["energy_kcal_mol"] is not None]
    if rotamers_with_energy:
        min_energy = min(r["energy_kcal_mol"] for r in rotamers_with_energy)
        for r in rotamers:
            if r["energy_kcal_mol"] is not None:
                r["relative_energy_kcal_mol"] = r["energy_kcal_mol"] - min_energy

    return {
        "success": True,
        "smiles": smiles,
        "n_rotatable_bonds": len(rotatable_bonds),
        "n_rotamers": len(rotamers),
        "angle_step_deg": angle_step,
        "rotamers": rotamers
    }


def generate_ligand_conformers(
    smiles: str,
    binding_site_shape: str = "cavity",
    n_conformers: int = 50,
    diversity_threshold: float = 0.5
) -> Dict[str, Any]:
    """
    Generate conformers suitable for docking studies.

    Focuses on extended conformers for cavity binding or
    compact conformers for surface binding.

    Args:
        smiles: SMILES string
        binding_site_shape: "cavity", "surface", "channel"
        n_conformers: Number of conformers to generate
        diversity_threshold: RMS threshold for diversity

    Returns:
        Docking-ready conformers
    """
    from rdkit import Chem
    from rdkit.Chem import AllChem, Descriptors

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return {
            "success": False,
            "error": {"code": "INVALID_SMILES", "message": f"Could not parse: {smiles}"}
        }

    mol = Chem.AddHs(mol)

    # Generate many conformers
    params = AllChem.ETKDGv3()
    params.useSmallRingTorsions = True
    params.numThreads = 0  # Use all cores

    n_to_generate = n_conformers * 5
    AllChem.EmbedMultipleConfs(mol, numConfs=n_to_generate, params=params)

    # Optimize and calculate properties
    conformers = []

    for conf_id in range(mol.GetNumConformers()):
        # Optimize
        ff = AllChem.MMFFGetMoleculeForceField(
            mol, AllChem.MMFFGetMoleculeProperties(mol), confId=conf_id
        )
        if ff:
            ff.Minimize(maxIts=200)
            energy = ff.CalcEnergy()
        else:
            continue

        # Calculate shape descriptors
        conf = mol.GetConformer(conf_id)
        coords = np.array([list(conf.GetAtomPosition(i)) for i in range(mol.GetNumAtoms())])

        # Calculate radius of gyration
        center = coords.mean(axis=0)
        rg = np.sqrt(np.mean(np.sum((coords - center) ** 2, axis=1)))

        # Calculate principal moments (simplified)
        distances = np.linalg.norm(coords - center, axis=1)
        max_extent = distances.max()

        # Asphericity (0 = spherical, 1 = linear)
        asphericity = 1 - (rg / max_extent) if max_extent > 0 else 0

        conformers.append({
            "conf_id": conf_id,
            "energy_kcal_mol": float(energy),
            "radius_of_gyration_A": float(rg),
            "max_extent_A": float(max_extent),
            "asphericity": float(asphericity),
            "coords": coords.tolist(),
            "symbols": [mol.GetAtomWithIdx(i).GetSymbol() for i in range(mol.GetNumAtoms())]
        })

    # Select based on binding site shape
    if binding_site_shape == "cavity":
        # Prefer compact conformers
        conformers.sort(key=lambda x: x["radius_of_gyration_A"])
    elif binding_site_shape == "channel":
        # Prefer extended conformers
        conformers.sort(key=lambda x: -x["asphericity"])
    else:  # surface
        # Prefer flat conformers (low asphericity, medium extension)
        conformers.sort(key=lambda x: x["asphericity"])

    # Select diverse conformers
    selected = []
    for conf in conformers:
        is_diverse = True
        for sel in selected:
            coords1 = np.array(conf["coords"])
            coords2 = np.array(sel["coords"])
            rms = np.sqrt(np.mean((coords1 - coords2) ** 2))
            if rms < diversity_threshold:
                is_diverse = False
                break
        if is_diverse:
            selected.append(conf)
        if len(selected) >= n_conformers:
            break

    return {
        "success": True,
        "smiles": smiles,
        "binding_site_shape": binding_site_shape,
        "n_conformers_requested": n_conformers,
        "n_conformers_generated": len(selected),
        "conformers": selected
    }
