"""
Annotated Molecule Data Structures

Provides rich molecular representation with:
- Stereochemistry information (R/S, E/Z, axial/equatorial)
- Ring system analysis
- Functional group detection
- Atom labeling
- Edit history tracking
"""

from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple, Any, Union
from enum import Enum
import importlib.util
import logging

# Check module availability
RDKIT_AVAILABLE = importlib.util.find_spec("rdkit") is not None

# Conditional imports
if RDKIT_AVAILABLE:
    from rdkit import Chem
    from rdkit.Chem import AllChem, Descriptors, rdMolDescriptors
    from rdkit.Chem.rdchem import ChiralType, BondStereo

logger = logging.getLogger(__name__)


class StereoConfiguration(Enum):
    """Stereochemical configuration types."""
    R = "R"
    S = "S"
    UNKNOWN = "unknown"
    NONE = "none"


class DoubleBondStereo(Enum):
    """Double bond stereochemistry."""
    E = "E"
    Z = "Z"
    UNKNOWN = "unknown"
    NONE = "none"


class RingType(Enum):
    """Common ring system types."""
    BENZENE = "benzene"
    NAPHTHALENE = "naphthalene"
    CYCLOHEXANE = "cyclohexane"
    CYCLOPENTANE = "cyclopentane"
    FURAN = "furan"
    THIOPHENE = "thiophene"
    PYRROLE = "pyrrole"
    PYRIDINE = "pyridine"
    IMIDAZOLE = "imidazole"
    OXETANE = "oxetane"
    EPOXIDE = "epoxide"
    LACTONE = "lactone"
    LACTAM = "lactam"
    MACROCYCLE = "macrocycle"
    UNKNOWN = "unknown"


class FunctionalGroupType(Enum):
    """Common functional group types."""
    HYDROXYL = "hydroxyl"
    CARBOXYL = "carboxyl"
    AMINO = "amino"
    CARBONYL = "carbonyl"
    ESTER = "ester"
    ETHER = "ether"
    AMIDE = "amide"
    NITRO = "nitro"
    NITRILE = "nitrile"
    SULFHYDRYL = "sulfhydryl"
    PHOSPHATE = "phosphate"
    HALIDE = "halide"
    METHYL = "methyl"
    ETHYL = "ethyl"
    PHENYL = "phenyl"
    ACETYL = "acetyl"
    UNKNOWN = "unknown"


@dataclass
class AtomInfo:
    """Information about a single atom."""
    idx: int
    element: str
    formal_charge: int = 0
    hybridization: str = "SP3"
    is_aromatic: bool = False
    num_hydrogens: int = 0
    in_ring: bool = False
    ring_sizes: List[int] = field(default_factory=list)

    def to_dict(self) -> Dict[str, Any]:
        return {
            "idx": self.idx,
            "element": self.element,
            "formal_charge": self.formal_charge,
            "hybridization": self.hybridization,
            "is_aromatic": self.is_aromatic,
            "num_hydrogens": self.num_hydrogens,
            "in_ring": self.in_ring,
            "ring_sizes": self.ring_sizes,
        }

    @classmethod
    def from_dict(cls, data: Dict[str, Any]) -> "AtomInfo":
        return cls(
            idx=data["idx"],
            element=data["element"],
            formal_charge=data.get("formal_charge", 0),
            hybridization=data.get("hybridization", "SP3"),
            is_aromatic=data.get("is_aromatic", False),
            num_hydrogens=data.get("num_hydrogens", 0),
            in_ring=data.get("in_ring", False),
            ring_sizes=data.get("ring_sizes", []),
        )


@dataclass
class StereocenterInfo:
    """Information about a stereocenter."""
    atom_idx: int
    configuration: StereoConfiguration
    neighbors: List[int] = field(default_factory=list)
    is_chiral: bool = True
    cip_code: Optional[str] = None

    def to_dict(self) -> Dict[str, Any]:
        return {
            "atom_idx": self.atom_idx,
            "configuration": self.configuration.value,
            "neighbors": self.neighbors,
            "is_chiral": self.is_chiral,
            "cip_code": self.cip_code,
        }

    @classmethod
    def from_dict(cls, data: Dict[str, Any]) -> "StereocenterInfo":
        config_value = data.get("configuration", "unknown")
        config = StereoConfiguration(config_value) if config_value in [e.value for e in StereoConfiguration] else StereoConfiguration.UNKNOWN
        return cls(
            atom_idx=data["atom_idx"],
            configuration=config,
            neighbors=data.get("neighbors", []),
            is_chiral=data.get("is_chiral", True),
            cip_code=data.get("cip_code"),
        )


@dataclass
class DoubleBondStereoInfo:
    """Information about double bond stereochemistry."""
    atom1_idx: int
    atom2_idx: int
    configuration: DoubleBondStereo
    substituents: Tuple[int, int, int, int] = (0, 0, 0, 0)  # 4 substituent atoms

    def to_dict(self) -> Dict[str, Any]:
        return {
            "atom1_idx": self.atom1_idx,
            "atom2_idx": self.atom2_idx,
            "configuration": self.configuration.value,
            "substituents": list(self.substituents),
        }

    @classmethod
    def from_dict(cls, data: Dict[str, Any]) -> "DoubleBondStereoInfo":
        config_value = data.get("configuration", "unknown")
        config = DoubleBondStereo(config_value) if config_value in [e.value for e in DoubleBondStereo] else DoubleBondStereo.UNKNOWN
        return cls(
            atom1_idx=data["atom1_idx"],
            atom2_idx=data["atom2_idx"],
            configuration=config,
            substituents=tuple(data.get("substituents", [0, 0, 0, 0])),
        )


@dataclass
class RingInfo:
    """Information about a ring system."""
    ring_id: int
    atom_indices: List[int]
    size: int
    ring_type: RingType
    is_aromatic: bool
    is_fused: bool = False
    fused_with: List[int] = field(default_factory=list)
    bridgehead_atoms: List[int] = field(default_factory=list)

    def to_dict(self) -> Dict[str, Any]:
        return {
            "ring_id": self.ring_id,
            "atom_indices": self.atom_indices,
            "size": self.size,
            "ring_type": self.ring_type.value,
            "is_aromatic": self.is_aromatic,
            "is_fused": self.is_fused,
            "fused_with": self.fused_with,
            "bridgehead_atoms": self.bridgehead_atoms,
        }

    @classmethod
    def from_dict(cls, data: Dict[str, Any]) -> "RingInfo":
        ring_type_value = data.get("ring_type", "unknown")
        ring_type = RingType(ring_type_value) if ring_type_value in [e.value for e in RingType] else RingType.UNKNOWN
        return cls(
            ring_id=data["ring_id"],
            atom_indices=data["atom_indices"],
            size=data["size"],
            ring_type=ring_type,
            is_aromatic=data.get("is_aromatic", False),
            is_fused=data.get("is_fused", False),
            fused_with=data.get("fused_with", []),
            bridgehead_atoms=data.get("bridgehead_atoms", []),
        )


@dataclass
class FunctionalGroupInfo:
    """Information about a functional group."""
    group_id: int
    group_type: FunctionalGroupType
    atom_indices: List[int]
    attachment_atom: int
    smarts_pattern: str = ""

    def to_dict(self) -> Dict[str, Any]:
        return {
            "group_id": self.group_id,
            "group_type": self.group_type.value,
            "atom_indices": self.atom_indices,
            "attachment_atom": self.attachment_atom,
            "smarts_pattern": self.smarts_pattern,
        }

    @classmethod
    def from_dict(cls, data: Dict[str, Any]) -> "FunctionalGroupInfo":
        group_type_value = data.get("group_type", "unknown")
        group_type = FunctionalGroupType(group_type_value) if group_type_value in [e.value for e in FunctionalGroupType] else FunctionalGroupType.UNKNOWN
        return cls(
            group_id=data["group_id"],
            group_type=group_type,
            atom_indices=data["atom_indices"],
            attachment_atom=data["attachment_atom"],
            smarts_pattern=data.get("smarts_pattern", ""),
        )


# SMARTS patterns for functional group detection
FUNCTIONAL_GROUP_PATTERNS = {
    FunctionalGroupType.HYDROXYL: "[OX2H]",
    FunctionalGroupType.CARBOXYL: "[CX3](=O)[OX2H1]",
    FunctionalGroupType.AMINO: "[NX3;H2,H1;!$(NC=O)]",
    FunctionalGroupType.CARBONYL: "[CX3]=[OX1]",
    FunctionalGroupType.ESTER: "[CX3](=O)[OX2][#6]",
    FunctionalGroupType.ETHER: "[OD2]([#6])[#6]",
    FunctionalGroupType.AMIDE: "[NX3][CX3](=[OX1])[#6]",
    FunctionalGroupType.NITRO: "[$([NX3](=O)=O),$([NX3+](=O)[O-])]",
    FunctionalGroupType.NITRILE: "[NX1]#[CX2]",
    FunctionalGroupType.SULFHYDRYL: "[SX2H]",
    FunctionalGroupType.HALIDE: "[F,Cl,Br,I]",
    FunctionalGroupType.METHYL: "[CH3]",
    FunctionalGroupType.PHENYL: "c1ccccc1",
    FunctionalGroupType.ACETYL: "[CX3](=O)[CH3]",
}

# SMARTS patterns for ring type detection
RING_TYPE_PATTERNS = {
    RingType.BENZENE: "c1ccccc1",
    RingType.NAPHTHALENE: "c1ccc2ccccc2c1",
    RingType.FURAN: "c1ccoc1",
    RingType.THIOPHENE: "c1ccsc1",
    RingType.PYRROLE: "c1cc[nH]c1",
    RingType.PYRIDINE: "c1ccncc1",
    RingType.IMIDAZOLE: "c1cnc[nH]1",
    RingType.OXETANE: "C1COC1",
    RingType.EPOXIDE: "C1OC1",
}


@dataclass
class AnnotatedMolecule:
    """
    Rich molecular representation with comprehensive annotations.

    This is the core data structure for molecular editing, containing:
    - SMILES representation
    - RDKit molecule object
    - Stereochemistry information
    - Ring system analysis
    - Functional group detection
    - Atom labeling
    - Edit history
    """
    smiles: str
    canonical_smiles: str = ""
    inchi: str = ""
    formula: str = ""
    molecular_weight: float = 0.0

    # Atom information
    atoms: List[AtomInfo] = field(default_factory=list)
    atom_labels: Dict[int, str] = field(default_factory=dict)

    # 3D coordinates (if available)
    coordinates: Optional[List[Tuple[float, float, float]]] = None
    has_3d: bool = False

    # Stereochemistry
    stereocenters: List[StereocenterInfo] = field(default_factory=list)
    double_bond_stereo: List[DoubleBondStereoInfo] = field(default_factory=list)

    # Ring systems
    rings: List[RingInfo] = field(default_factory=list)

    # Functional groups
    functional_groups: List[FunctionalGroupInfo] = field(default_factory=list)

    # Edit history
    edit_history: List[Dict[str, Any]] = field(default_factory=list)

    # Internal: RDKit mol object (not serialized)
    _mol: Any = field(default=None, repr=False)

    def get_mol(self) -> Optional[Any]:
        """Get RDKit molecule object, creating from SMILES if needed."""
        if self._mol is not None:
            return self._mol
        if not RDKIT_AVAILABLE:
            return None
        mol = Chem.MolFromSmiles(self.smiles)
        if mol is not None:
            self._mol = mol
        return self._mol

    def get_rwmol(self) -> Optional[Any]:
        """Get editable RWMol object."""
        mol = self.get_mol()
        if mol is None:
            return None
        return Chem.RWMol(mol)

    def to_dict(self) -> Dict[str, Any]:
        """Serialize to dictionary for JSON transport."""
        return {
            "smiles": self.smiles,
            "canonical_smiles": self.canonical_smiles,
            "inchi": self.inchi,
            "formula": self.formula,
            "molecular_weight": self.molecular_weight,
            "atoms": [a.to_dict() for a in self.atoms],
            "atom_labels": self.atom_labels,
            "coordinates": self.coordinates,
            "has_3d": self.has_3d,
            "stereocenters": [s.to_dict() for s in self.stereocenters],
            "double_bond_stereo": [d.to_dict() for d in self.double_bond_stereo],
            "rings": [r.to_dict() for r in self.rings],
            "functional_groups": [f.to_dict() for f in self.functional_groups],
            "edit_history": self.edit_history,
        }

    @classmethod
    def from_dict(cls, data: Dict[str, Any]) -> "AnnotatedMolecule":
        """Deserialize from dictionary."""
        return cls(
            smiles=data["smiles"],
            canonical_smiles=data.get("canonical_smiles", ""),
            inchi=data.get("inchi", ""),
            formula=data.get("formula", ""),
            molecular_weight=data.get("molecular_weight", 0.0),
            atoms=[AtomInfo.from_dict(a) for a in data.get("atoms", [])],
            atom_labels=data.get("atom_labels", {}),
            coordinates=data.get("coordinates"),
            has_3d=data.get("has_3d", False),
            stereocenters=[StereocenterInfo.from_dict(s) for s in data.get("stereocenters", [])],
            double_bond_stereo=[DoubleBondStereoInfo.from_dict(d) for d in data.get("double_bond_stereo", [])],
            rings=[RingInfo.from_dict(r) for r in data.get("rings", [])],
            functional_groups=[FunctionalGroupInfo.from_dict(f) for f in data.get("functional_groups", [])],
            edit_history=data.get("edit_history", []),
        )

    @classmethod
    def from_smiles(cls, smiles: str, auto_annotate: bool = True) -> Optional["AnnotatedMolecule"]:
        """
        Create AnnotatedMolecule from SMILES string with automatic annotation.

        Args:
            smiles: SMILES string
            auto_annotate: Whether to automatically detect stereocenters, rings, etc.

        Returns:
            AnnotatedMolecule or None if SMILES is invalid
        """
        if not RDKIT_AVAILABLE:
            logger.error("RDKit not available for SMILES parsing")
            return None

        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            logger.error(f"Invalid SMILES: {smiles}")
            return None

        # Get canonical SMILES
        canonical = Chem.MolToSmiles(mol, canonical=True)

        # Calculate molecular properties
        formula = rdMolDescriptors.CalcMolFormula(mol)
        mw = Descriptors.MolWt(mol)

        # Create molecule with basic info
        annotated = cls(
            smiles=smiles,
            canonical_smiles=canonical,
            formula=formula,
            molecular_weight=mw,
            _mol=mol,
        )

        if auto_annotate:
            annotated._annotate_atoms()
            annotated._annotate_stereocenters()
            annotated._annotate_double_bond_stereo()
            annotated._annotate_rings()
            annotated._annotate_functional_groups()
            annotated._generate_atom_labels()

        return annotated

    def _annotate_atoms(self) -> None:
        """Extract atom information from RDKit molecule."""
        mol = self.get_mol()
        if mol is None:
            return

        ring_info = mol.GetRingInfo()
        atom_rings = ring_info.AtomRings()

        self.atoms = []
        for atom in mol.GetAtoms():
            idx = atom.GetIdx()

            # Get ring sizes for this atom
            ring_sizes = [len(ring) for ring in atom_rings if idx in ring]

            info = AtomInfo(
                idx=idx,
                element=atom.GetSymbol(),
                formal_charge=atom.GetFormalCharge(),
                hybridization=str(atom.GetHybridization()),
                is_aromatic=atom.GetIsAromatic(),
                num_hydrogens=atom.GetTotalNumHs(),
                in_ring=atom.IsInRing(),
                ring_sizes=ring_sizes,
            )
            self.atoms.append(info)

    def _annotate_stereocenters(self) -> None:
        """Detect and annotate stereocenters."""
        mol = self.get_mol()
        if mol is None:
            return

        # Assign stereochemistry
        Chem.AssignStereochemistry(mol, cleanIt=True, force=True)

        # Find chiral centers
        chiral_centers = Chem.FindMolChiralCenters(mol, includeUnassigned=True)

        self.stereocenters = []
        for atom_idx, cip_code in chiral_centers:
            atom = mol.GetAtomWithIdx(atom_idx)
            neighbors = [n.GetIdx() for n in atom.GetNeighbors()]

            config = StereoConfiguration.UNKNOWN
            if cip_code == "R":
                config = StereoConfiguration.R
            elif cip_code == "S":
                config = StereoConfiguration.S
            elif cip_code == "?":
                config = StereoConfiguration.UNKNOWN

            self.stereocenters.append(StereocenterInfo(
                atom_idx=atom_idx,
                configuration=config,
                neighbors=neighbors,
                is_chiral=True,
                cip_code=cip_code,
            ))

    def _annotate_double_bond_stereo(self) -> None:
        """Detect and annotate double bond stereochemistry."""
        mol = self.get_mol()
        if mol is None:
            return

        self.double_bond_stereo = []
        for bond in mol.GetBonds():
            if bond.GetBondType() != Chem.BondType.DOUBLE:
                continue

            stereo = bond.GetStereo()
            if stereo == BondStereo.STEREONONE:
                continue

            config = DoubleBondStereo.UNKNOWN
            if stereo == BondStereo.STEREOE:
                config = DoubleBondStereo.E
            elif stereo == BondStereo.STEREOZ:
                config = DoubleBondStereo.Z

            atom1_idx = bond.GetBeginAtomIdx()
            atom2_idx = bond.GetEndAtomIdx()

            # Get stereo atoms
            stereo_atoms = bond.GetStereoAtoms()
            substituents = tuple(stereo_atoms) if len(stereo_atoms) == 4 else (0, 0, 0, 0)

            self.double_bond_stereo.append(DoubleBondStereoInfo(
                atom1_idx=atom1_idx,
                atom2_idx=atom2_idx,
                configuration=config,
                substituents=substituents,
            ))

    def _annotate_rings(self) -> None:
        """Detect and annotate ring systems."""
        mol = self.get_mol()
        if mol is None:
            return

        ring_info = mol.GetRingInfo()
        atom_rings = ring_info.AtomRings()

        self.rings = []
        for ring_id, ring_atoms in enumerate(atom_rings):
            atom_list = list(ring_atoms)
            size = len(atom_list)

            # Check aromaticity
            is_aromatic = all(mol.GetAtomWithIdx(idx).GetIsAromatic() for idx in atom_list)

            # Determine ring type
            ring_type = self._classify_ring(mol, atom_list)

            # Check fusion
            fused_with = []
            for other_id, other_atoms in enumerate(atom_rings):
                if other_id != ring_id:
                    shared = set(atom_list) & set(other_atoms)
                    if len(shared) >= 2:
                        fused_with.append(other_id)

            self.rings.append(RingInfo(
                ring_id=ring_id,
                atom_indices=atom_list,
                size=size,
                ring_type=ring_type,
                is_aromatic=is_aromatic,
                is_fused=len(fused_with) > 0,
                fused_with=fused_with,
            ))

    def _classify_ring(self, mol, atom_indices: List[int]) -> RingType:
        """Classify a ring based on its structure."""
        size = len(atom_indices)

        # Get elements in ring
        elements = [mol.GetAtomWithIdx(idx).GetSymbol() for idx in atom_indices]
        is_aromatic = all(mol.GetAtomWithIdx(idx).GetIsAromatic() for idx in atom_indices)

        # 3-membered rings
        if size == 3:
            if "O" in elements and elements.count("C") == 2:
                return RingType.EPOXIDE
            return RingType.UNKNOWN

        # 4-membered rings
        if size == 4:
            if "O" in elements and elements.count("C") == 3:
                return RingType.OXETANE
            return RingType.CYCLOHEXANE  # Generic 4-membered

        # 5-membered rings
        if size == 5:
            if is_aromatic:
                if "O" in elements:
                    return RingType.FURAN
                if "S" in elements:
                    return RingType.THIOPHENE
                if "N" in elements:
                    if elements.count("N") == 1:
                        return RingType.PYRROLE
                    elif elements.count("N") == 2:
                        return RingType.IMIDAZOLE
            return RingType.CYCLOPENTANE

        # 6-membered rings
        if size == 6:
            if is_aromatic:
                if "N" in elements:
                    return RingType.PYRIDINE
                return RingType.BENZENE
            return RingType.CYCLOHEXANE

        # Larger rings
        if size >= 8:
            return RingType.MACROCYCLE

        return RingType.UNKNOWN

    def _annotate_functional_groups(self) -> None:
        """Detect functional groups using SMARTS patterns."""
        mol = self.get_mol()
        if mol is None:
            return

        self.functional_groups = []
        group_id = 0

        for group_type, smarts in FUNCTIONAL_GROUP_PATTERNS.items():
            pattern = Chem.MolFromSmarts(smarts)
            if pattern is None:
                continue

            matches = mol.GetSubstructMatches(pattern)
            for match in matches:
                atom_indices = list(match)

                # Determine attachment atom (first atom in pattern)
                attachment = atom_indices[0] if atom_indices else -1

                self.functional_groups.append(FunctionalGroupInfo(
                    group_id=group_id,
                    group_type=group_type,
                    atom_indices=atom_indices,
                    attachment_atom=attachment,
                    smarts_pattern=smarts,
                ))
                group_id += 1

    def _generate_atom_labels(self) -> None:
        """Generate standard atom labels (C-1, C-2, N-1, etc.)."""
        mol = self.get_mol()
        if mol is None:
            return

        element_counts: Dict[str, int] = {}
        self.atom_labels = {}

        for atom in mol.GetAtoms():
            symbol = atom.GetSymbol()
            count = element_counts.get(symbol, 0) + 1
            element_counts[symbol] = count
            self.atom_labels[atom.GetIdx()] = f"{symbol}-{count}"

    def get_atom_by_label(self, label: str) -> Optional[int]:
        """Get atom index by label (e.g., 'C-7')."""
        for idx, lbl in self.atom_labels.items():
            if lbl == label:
                return idx
        return None

    def get_atoms_by_smarts(self, smarts: str) -> List[Tuple[int, ...]]:
        """Get atom indices matching a SMARTS pattern."""
        mol = self.get_mol()
        if mol is None:
            return []

        pattern = Chem.MolFromSmarts(smarts)
        if pattern is None:
            return []

        return list(mol.GetSubstructMatches(pattern))

    def get_ring_by_type(self, ring_type: Union[str, RingType]) -> List[RingInfo]:
        """Get rings of a specific type."""
        if isinstance(ring_type, str):
            ring_type = RingType(ring_type) if ring_type in [e.value for e in RingType] else RingType.UNKNOWN
        return [r for r in self.rings if r.ring_type == ring_type]

    def get_functional_group_by_type(self, group_type: Union[str, FunctionalGroupType]) -> List[FunctionalGroupInfo]:
        """Get functional groups of a specific type."""
        if isinstance(group_type, str):
            group_type = FunctionalGroupType(group_type) if group_type in [e.value for e in FunctionalGroupType] else FunctionalGroupType.UNKNOWN
        return [g for g in self.functional_groups if g.group_type == group_type]

    def add_edit_to_history(self, operation_type: str, target: str, params: Dict[str, Any], success: bool) -> None:
        """Record an edit operation in history."""
        import datetime
        self.edit_history.append({
            "operation": operation_type,
            "target": target,
            "params": params,
            "success": success,
            "timestamp": datetime.datetime.now().isoformat(),
        })

    def update_from_mol(self, mol) -> None:
        """Update annotations after molecular modifications."""
        self._mol = mol
        self.smiles = Chem.MolToSmiles(mol)
        self.canonical_smiles = Chem.MolToSmiles(mol, canonical=True)
        self.formula = rdMolDescriptors.CalcMolFormula(mol)
        self.molecular_weight = Descriptors.MolWt(mol)

        # Re-annotate
        self._annotate_atoms()
        self._annotate_stereocenters()
        self._annotate_double_bond_stereo()
        self._annotate_rings()
        self._annotate_functional_groups()
        self._generate_atom_labels()
