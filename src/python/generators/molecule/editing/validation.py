"""
Molecular Validation Module

Provides post-edit validation pipeline for molecular structures:
- Valence checking
- Aromaticity consistency
- Ring strain analysis
- Stereochemistry validity
- Steric clash detection
- Chemical reasonability checks

All validation uses guard clauses and explicit checks - no try/except blocks.
"""

from dataclasses import dataclass, field
from typing import Dict, List, Optional, Set, Tuple, Any, Callable
from enum import Enum
import logging
import importlib.util

# Check module availability
RDKIT_AVAILABLE = importlib.util.find_spec("rdkit") is not None

if RDKIT_AVAILABLE:
    from rdkit import Chem
    from rdkit.Chem import Descriptors

logger = logging.getLogger(__name__)


class ValidationSeverity(Enum):
    """Severity levels for validation issues."""
    INFO = "info"
    WARNING = "warning"
    ERROR = "error"
    CRITICAL = "critical"


class ValidationType(Enum):
    """Types of validation checks."""
    VALENCE = "valence"
    AROMATICITY = "aromaticity"
    RING_STRAIN = "ring_strain"
    STEREOCHEMISTRY = "stereochemistry"
    STERIC_CLASH = "steric_clash"
    CONNECTIVITY = "connectivity"
    FORMAL_CHARGE = "formal_charge"
    HYBRIDIZATION = "hybridization"
    BOND_GEOMETRY = "bond_geometry"
    CHEMICAL_REASONABILITY = "chemical_reasonability"
    SANITIZATION = "sanitization"


@dataclass
class ValidationIssue:
    """Represents a single validation issue."""
    validation_type: ValidationType
    severity: ValidationSeverity
    message: str
    atom_indices: List[int] = field(default_factory=list)
    bond_indices: List[int] = field(default_factory=list)
    details: Dict[str, Any] = field(default_factory=dict)

    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary representation."""
        return {
            "type": self.validation_type.value,
            "severity": self.severity.value,
            "message": self.message,
            "atom_indices": self.atom_indices,
            "bond_indices": self.bond_indices,
            "details": self.details,
        }

    @classmethod
    def from_dict(cls, data: Dict[str, Any]) -> "ValidationIssue":
        """Create from dictionary."""
        return cls(
            validation_type=ValidationType(data.get("type", "sanitization")),
            severity=ValidationSeverity(data.get("severity", "warning")),
            message=data.get("message", "Unknown issue"),
            atom_indices=data.get("atom_indices", []),
            bond_indices=data.get("bond_indices", []),
            details=data.get("details", {}),
        )


@dataclass
class ValidationResult:
    """Result of validation pipeline execution."""
    valid: bool
    issues: List[ValidationIssue] = field(default_factory=list)
    warnings_count: int = 0
    errors_count: int = 0
    critical_count: int = 0
    checks_performed: List[str] = field(default_factory=list)
    molecule_info: Dict[str, Any] = field(default_factory=dict)

    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary representation."""
        return {
            "valid": self.valid,
            "issues": [issue.to_dict() for issue in self.issues],
            "warnings_count": self.warnings_count,
            "errors_count": self.errors_count,
            "critical_count": self.critical_count,
            "checks_performed": self.checks_performed,
            "molecule_info": self.molecule_info,
        }

    @classmethod
    def from_dict(cls, data: Dict[str, Any]) -> "ValidationResult":
        """Create from dictionary."""
        issues = [ValidationIssue.from_dict(i) for i in data.get("issues", [])]
        return cls(
            valid=data.get("valid", True),
            issues=issues,
            warnings_count=data.get("warnings_count", 0),
            errors_count=data.get("errors_count", 0),
            critical_count=data.get("critical_count", 0),
            checks_performed=data.get("checks_performed", []),
            molecule_info=data.get("molecule_info", {}),
        )

    def add_issue(self, issue: ValidationIssue) -> None:
        """Add an issue and update counts."""
        self.issues.append(issue)
        if issue.severity == ValidationSeverity.WARNING:
            self.warnings_count += 1
        elif issue.severity == ValidationSeverity.ERROR:
            self.errors_count += 1
            self.valid = False
        elif issue.severity == ValidationSeverity.CRITICAL:
            self.critical_count += 1
            self.valid = False

    def get_issues_by_type(self, validation_type: ValidationType) -> List[ValidationIssue]:
        """Get all issues of a specific type."""
        return [i for i in self.issues if i.validation_type == validation_type]

    def get_issues_by_severity(self, severity: ValidationSeverity) -> List[ValidationIssue]:
        """Get all issues of a specific severity."""
        return [i for i in self.issues if i.severity == severity]


# Standard valence expectations for common elements
STANDARD_VALENCES: Dict[str, List[int]] = {
    "H": [1],
    "C": [4],
    "N": [3, 5],  # Nitrogen can be trivalent or pentavalent
    "O": [2],
    "F": [1],
    "P": [3, 5],  # Phosphorus can be trivalent or pentavalent
    "S": [2, 4, 6],  # Sulfur has multiple valence states
    "Cl": [1, 3, 5, 7],
    "Br": [1, 3, 5, 7],
    "I": [1, 3, 5, 7],
    "B": [3],
    "Si": [4],
    "Se": [2, 4, 6],
    "As": [3, 5],
}

# Ring strain energies (kcal/mol) for small rings
RING_STRAIN_ENERGIES: Dict[int, float] = {
    3: 27.5,   # Cyclopropane - very strained
    4: 26.3,   # Cyclobutane - very strained
    5: 6.2,    # Cyclopentane - slightly strained
    6: 0.0,    # Cyclohexane - no strain
    7: 6.2,    # Cycloheptane - slight strain
    8: 9.7,    # Cyclooctane - transannular strain
}

# Strain thresholds
STRAIN_WARNING_THRESHOLD = 10.0  # kcal/mol
STRAIN_ERROR_THRESHOLD = 30.0   # kcal/mol

# Steric clash distance threshold (Angstroms)
STERIC_CLASH_THRESHOLD = 1.5  # Atoms closer than this are clashing

# Common unreasonable substructures (SMARTS patterns)
UNREASONABLE_PATTERNS: Dict[str, str] = {
    "pentavalent_carbon": "[#6X5]",
    "hexavalent_carbon": "[#6X6]",
    "peroxide_chain": "[O][O][O]",  # Polyperoxide
    "carbene": "[#6v2]",  # Carbene (unless stabilized)
    "nitrene": "[#7v1]",  # Nitrene
    "allene_stack": "C=C=C=C=C",  # Cumulated double bonds
}


def check_valence(mol) -> List[ValidationIssue]:
    """
    Check valence validity for all atoms.

    Args:
        mol: RDKit molecule object

    Returns:
        List of ValidationIssue objects for valence problems
    """
    if not RDKIT_AVAILABLE:
        return []
    if mol is None:
        return [ValidationIssue(
            validation_type=ValidationType.VALENCE,
            severity=ValidationSeverity.CRITICAL,
            message="No molecule provided for valence check",
        )]

    issues = []

    for atom in mol.GetAtoms():
        idx = atom.GetIdx()
        symbol = atom.GetSymbol()

        # Get total valence including implicit hydrogens
        total_valence = atom.GetTotalValence()
        formal_charge = atom.GetFormalCharge()

        # Adjust expected valence for formal charge
        expected_valences = STANDARD_VALENCES.get(symbol, [])

        if not expected_valences:
            # Unknown element - skip but note
            continue

        # Check if valence is reasonable
        # Positive charge reduces expected valence, negative increases
        adjusted_valences = []
        for v in expected_valences:
            adjusted = v - formal_charge
            if adjusted > 0:
                adjusted_valences.append(adjusted)

        if not adjusted_valences:
            adjusted_valences = expected_valences

        is_valid = total_valence in adjusted_valences

        # Allow some tolerance for aromatic systems
        if not is_valid and atom.GetIsAromatic():
            # Aromatic atoms can have different valence patterns
            is_valid = total_valence in [v - 1 for v in expected_valences] or \
                       total_valence in [v + 1 for v in expected_valences]

        if not is_valid:
            severity = ValidationSeverity.ERROR
            if total_valence > max(expected_valences) + 2:
                severity = ValidationSeverity.CRITICAL

            issues.append(ValidationIssue(
                validation_type=ValidationType.VALENCE,
                severity=severity,
                message=f"Atom {symbol} at index {idx} has valence {total_valence}, "
                       f"expected {expected_valences}",
                atom_indices=[idx],
                details={
                    "element": symbol,
                    "actual_valence": total_valence,
                    "expected_valences": expected_valences,
                    "formal_charge": formal_charge,
                },
            ))

    return issues


def check_aromaticity(mol) -> List[ValidationIssue]:
    """
    Check aromaticity consistency.

    Args:
        mol: RDKit molecule object

    Returns:
        List of ValidationIssue objects for aromaticity problems
    """
    if not RDKIT_AVAILABLE:
        return []
    if mol is None:
        return []

    issues = []

    # Get ring info
    ring_info = mol.GetRingInfo()
    atom_rings = ring_info.AtomRings()

    if not atom_rings:
        return []

    for ring_atoms in atom_rings:
        if len(ring_atoms) < 3:
            continue

        # Count aromatic atoms in ring
        aromatic_count = sum(1 for idx in ring_atoms if mol.GetAtomWithIdx(idx).GetIsAromatic())

        # Ring should be either fully aromatic or not aromatic at all
        if 0 < aromatic_count < len(ring_atoms):
            issues.append(ValidationIssue(
                validation_type=ValidationType.AROMATICITY,
                severity=ValidationSeverity.WARNING,
                message=f"Ring has partial aromaticity: {aromatic_count}/{len(ring_atoms)} atoms aromatic",
                atom_indices=list(ring_atoms),
                details={
                    "ring_size": len(ring_atoms),
                    "aromatic_atoms": aromatic_count,
                },
            ))

        # Check Huckel's rule for aromatic rings (4n+2 pi electrons)
        if aromatic_count == len(ring_atoms):
            pi_electrons = 0
            for idx in ring_atoms:
                atom = mol.GetAtomWithIdx(idx)
                symbol = atom.GetSymbol()

                # Simplified pi electron counting
                if symbol == "C":
                    pi_electrons += 1
                elif symbol == "N":
                    # Nitrogen contributes 1 or 2 depending on if it has lone pair in ring
                    if atom.GetTotalNumHs() == 0:
                        pi_electrons += 2
                    else:
                        pi_electrons += 1
                elif symbol == "O":
                    pi_electrons += 2
                elif symbol == "S":
                    pi_electrons += 2
                else:
                    pi_electrons += 1

            # Check 4n+2 rule
            remainder = (pi_electrons - 2) % 4
            if remainder != 0:
                issues.append(ValidationIssue(
                    validation_type=ValidationType.AROMATICITY,
                    severity=ValidationSeverity.INFO,
                    message=f"Aromatic ring may not satisfy Hückel's rule: {pi_electrons} π electrons",
                    atom_indices=list(ring_atoms),
                    details={
                        "pi_electrons": pi_electrons,
                        "ring_size": len(ring_atoms),
                    },
                ))

    return issues


def check_ring_strain(mol) -> List[ValidationIssue]:
    """
    Analyze ring strain in the molecule.

    Args:
        mol: RDKit molecule object

    Returns:
        List of ValidationIssue objects for ring strain problems
    """
    if not RDKIT_AVAILABLE:
        return []
    if mol is None:
        return []

    issues = []

    ring_info = mol.GetRingInfo()
    atom_rings = ring_info.AtomRings()

    if not atom_rings:
        return []

    for ring_atoms in atom_rings:
        ring_size = len(ring_atoms)

        # Get strain energy for this ring size
        strain_energy = RING_STRAIN_ENERGIES.get(ring_size, 0.0)

        if strain_energy >= STRAIN_ERROR_THRESHOLD:
            issues.append(ValidationIssue(
                validation_type=ValidationType.RING_STRAIN,
                severity=ValidationSeverity.ERROR,
                message=f"Ring of size {ring_size} has high strain energy: {strain_energy:.1f} kcal/mol",
                atom_indices=list(ring_atoms),
                details={
                    "ring_size": ring_size,
                    "strain_energy": strain_energy,
                },
            ))
        elif strain_energy >= STRAIN_WARNING_THRESHOLD:
            issues.append(ValidationIssue(
                validation_type=ValidationType.RING_STRAIN,
                severity=ValidationSeverity.WARNING,
                message=f"Ring of size {ring_size} has notable strain: {strain_energy:.1f} kcal/mol",
                atom_indices=list(ring_atoms),
                details={
                    "ring_size": ring_size,
                    "strain_energy": strain_energy,
                },
            ))

        # Check for bridgehead double bonds (Bredt's rule)
        if ring_size <= 7:
            for atom_idx in ring_atoms:
                atom = mol.GetAtomWithIdx(atom_idx)
                # Check if bridgehead (in multiple rings)
                num_rings_containing = sum(1 for r in atom_rings if atom_idx in r)

                if num_rings_containing >= 2:
                    # Check for double bonds at bridgehead
                    for bond in atom.GetBonds():
                        if bond.GetBondType() == Chem.BondType.DOUBLE:
                            other_idx = bond.GetOtherAtomIdx(atom_idx)
                            # Check if the double bond is part of a small ring
                            for r in atom_rings:
                                if atom_idx in r and other_idx in r and len(r) <= 7:
                                    issues.append(ValidationIssue(
                                        validation_type=ValidationType.RING_STRAIN,
                                        severity=ValidationSeverity.WARNING,
                                        message=f"Potential Bredt's rule violation: "
                                               f"double bond at bridgehead in {len(r)}-membered ring",
                                        atom_indices=[atom_idx, other_idx],
                                        details={
                                            "ring_size": len(r),
                                            "bridgehead_atom": atom_idx,
                                        },
                                    ))
                                    break

    return issues


def check_stereochemistry(mol) -> List[ValidationIssue]:
    """
    Validate stereochemistry assignments.

    Args:
        mol: RDKit molecule object

    Returns:
        List of ValidationIssue objects for stereochemistry problems
    """
    if not RDKIT_AVAILABLE:
        return []
    if mol is None:
        return []

    issues = []

    # Find stereocenters
    chiral_centers = Chem.FindMolChiralCenters(mol, includeUnassigned=True)

    for atom_idx, stereo_type in chiral_centers:
        atom = mol.GetAtomWithIdx(atom_idx)

        # Check for unassigned stereocenters
        if stereo_type == "?":
            issues.append(ValidationIssue(
                validation_type=ValidationType.STEREOCHEMISTRY,
                severity=ValidationSeverity.WARNING,
                message=f"Unassigned stereocenter at atom {atom_idx} ({atom.GetSymbol()})",
                atom_indices=[atom_idx],
                details={
                    "element": atom.GetSymbol(),
                    "chirality": "unassigned",
                },
            ))

        # Check for proper tetrahedral geometry
        neighbors = list(atom.GetNeighbors())
        if len(neighbors) != 4:
            # Stereocenter should have 4 substituents (including implicit H)
            total_neighbors = len(neighbors) + atom.GetTotalNumHs()
            if total_neighbors != 4:
                issues.append(ValidationIssue(
                    validation_type=ValidationType.STEREOCHEMISTRY,
                    severity=ValidationSeverity.ERROR,
                    message=f"Stereocenter at atom {atom_idx} has {total_neighbors} substituents, "
                           f"expected 4 for tetrahedral geometry",
                    atom_indices=[atom_idx],
                    details={
                        "substituent_count": total_neighbors,
                        "expected": 4,
                    },
                ))

    # Check double bond stereochemistry
    for bond in mol.GetBonds():
        if bond.GetBondType() == Chem.BondType.DOUBLE:
            stereo = bond.GetStereo()

            # Get atoms
            begin_idx = bond.GetBeginAtomIdx()
            end_idx = bond.GetEndAtomIdx()
            begin_atom = mol.GetAtomWithIdx(begin_idx)
            end_atom = mol.GetAtomWithIdx(end_idx)

            # Check if E/Z stereochemistry is possible
            begin_neighbors = [n.GetIdx() for n in begin_atom.GetNeighbors() if n.GetIdx() != end_idx]
            end_neighbors = [n.GetIdx() for n in end_atom.GetNeighbors() if n.GetIdx() != begin_idx]

            # E/Z is only relevant if both atoms have 2 different substituents
            if len(begin_neighbors) >= 1 and len(end_neighbors) >= 1:
                if stereo == Chem.BondStereo.STEREONONE:
                    # Check if substituents are different
                    if len(set(begin_neighbors)) > 1 or len(set(end_neighbors)) > 1:
                        issues.append(ValidationIssue(
                            validation_type=ValidationType.STEREOCHEMISTRY,
                            severity=ValidationSeverity.INFO,
                            message=f"Double bond {begin_idx}-{end_idx} may have unspecified E/Z stereochemistry",
                            bond_indices=[bond.GetIdx()],
                            atom_indices=[begin_idx, end_idx],
                            details={
                                "stereo": "unspecified",
                            },
                        ))

    return issues


def check_steric_clashes(mol) -> List[ValidationIssue]:
    """
    Detect steric clashes in 3D structure.

    Args:
        mol: RDKit molecule object with 3D coordinates

    Returns:
        List of ValidationIssue objects for steric clash problems
    """
    if not RDKIT_AVAILABLE:
        return []
    if mol is None:
        return []

    # Check if molecule has 3D coordinates
    if mol.GetNumConformers() == 0:
        return []

    issues = []
    conf = mol.GetConformer()
    num_atoms = mol.GetNumAtoms()

    # Get bonded pairs to exclude from clash detection
    bonded_pairs: Set[Tuple[int, int]] = set()
    for bond in mol.GetBonds():
        i, j = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
        bonded_pairs.add((min(i, j), max(i, j)))

    # Get 1-3 pairs (atoms separated by one bond) to exclude
    one_three_pairs: Set[Tuple[int, int]] = set()
    for atom in mol.GetAtoms():
        neighbors = [n.GetIdx() for n in atom.GetNeighbors()]
        for i in range(len(neighbors)):
            for j in range(i + 1, len(neighbors)):
                pair = (min(neighbors[i], neighbors[j]), max(neighbors[i], neighbors[j]))
                one_three_pairs.add(pair)

    # Check distances between all non-bonded, non-1-3 atom pairs
    clashes: List[Tuple[int, int, float]] = []

    for i in range(num_atoms):
        pos_i = conf.GetAtomPosition(i)
        for j in range(i + 1, num_atoms):
            pair = (i, j)

            # Skip bonded pairs and 1-3 pairs
            if pair in bonded_pairs or pair in one_three_pairs:
                continue

            pos_j = conf.GetAtomPosition(j)
            distance = pos_i.Distance(pos_j)

            # Get van der Waals radii sum for these atoms
            atom_i = mol.GetAtomWithIdx(i)
            atom_j = mol.GetAtomWithIdx(j)

            # Approximate vdW radii
            vdw_i = Chem.GetPeriodicTable().GetRvdw(atom_i.GetAtomicNum())
            vdw_j = Chem.GetPeriodicTable().GetRvdw(atom_j.GetAtomicNum())

            # Clash if distance is less than sum of radii minus tolerance
            min_distance = (vdw_i + vdw_j) * 0.7  # 70% of sum as threshold

            if distance < min_distance:
                clashes.append((i, j, distance))

    # Report clashes
    for i, j, distance in clashes:
        atom_i = mol.GetAtomWithIdx(i)
        atom_j = mol.GetAtomWithIdx(j)

        severity = ValidationSeverity.WARNING
        if distance < STERIC_CLASH_THRESHOLD:
            severity = ValidationSeverity.ERROR

        issues.append(ValidationIssue(
            validation_type=ValidationType.STERIC_CLASH,
            severity=severity,
            message=f"Steric clash between {atom_i.GetSymbol()}{i} and {atom_j.GetSymbol()}{j}: "
                   f"{distance:.2f} Å",
            atom_indices=[i, j],
            details={
                "distance": distance,
                "atom1": {"idx": i, "element": atom_i.GetSymbol()},
                "atom2": {"idx": j, "element": atom_j.GetSymbol()},
            },
        ))

    return issues


def check_connectivity(mol) -> List[ValidationIssue]:
    """
    Check molecular connectivity for fragments and isolated atoms.

    Args:
        mol: RDKit molecule object

    Returns:
        List of ValidationIssue objects for connectivity problems
    """
    if not RDKIT_AVAILABLE:
        return []
    if mol is None:
        return []

    issues = []

    # Get molecular fragments
    fragments = Chem.GetMolFrags(mol)

    if len(fragments) > 1:
        issues.append(ValidationIssue(
            validation_type=ValidationType.CONNECTIVITY,
            severity=ValidationSeverity.WARNING,
            message=f"Molecule has {len(fragments)} disconnected fragments",
            details={
                "fragment_count": len(fragments),
                "fragment_sizes": [len(f) for f in fragments],
            },
        ))

        # Identify small fragments (might be leaving groups or artifacts)
        for i, frag in enumerate(fragments):
            if len(frag) == 1:
                atom_idx = frag[0]
                atom = mol.GetAtomWithIdx(atom_idx)
                issues.append(ValidationIssue(
                    validation_type=ValidationType.CONNECTIVITY,
                    severity=ValidationSeverity.WARNING,
                    message=f"Isolated atom detected: {atom.GetSymbol()} at index {atom_idx}",
                    atom_indices=[atom_idx],
                    details={
                        "fragment_index": i,
                        "element": atom.GetSymbol(),
                    },
                ))

    return issues


def check_formal_charges(mol) -> List[ValidationIssue]:
    """
    Check formal charge distribution.

    Args:
        mol: RDKit molecule object

    Returns:
        List of ValidationIssue objects for charge problems
    """
    if not RDKIT_AVAILABLE:
        return []
    if mol is None:
        return []

    issues = []

    total_charge = 0
    charged_atoms = []

    for atom in mol.GetAtoms():
        charge = atom.GetFormalCharge()
        if charge != 0:
            total_charge += charge
            charged_atoms.append((atom.GetIdx(), atom.GetSymbol(), charge))

    # Report high net charge
    if abs(total_charge) > 2:
        issues.append(ValidationIssue(
            validation_type=ValidationType.FORMAL_CHARGE,
            severity=ValidationSeverity.WARNING,
            message=f"High net formal charge: {total_charge:+d}",
            details={
                "total_charge": total_charge,
                "charged_atoms": charged_atoms,
            },
        ))

    # Check for unusual charges on specific elements
    unusual_charges = {
        "C": {-2, -3, 2, 3, 4},  # Unusual carbon charges
        "N": {-3, 3, 4},  # Unusual nitrogen charges
        "O": {2, -3},  # Unusual oxygen charges
    }

    for idx, symbol, charge in charged_atoms:
        if symbol in unusual_charges and charge in unusual_charges[symbol]:
            issues.append(ValidationIssue(
                validation_type=ValidationType.FORMAL_CHARGE,
                severity=ValidationSeverity.WARNING,
                message=f"Unusual formal charge {charge:+d} on {symbol} at index {idx}",
                atom_indices=[idx],
                details={
                    "element": symbol,
                    "charge": charge,
                },
            ))

    return issues


def check_chemical_reasonability(mol) -> List[ValidationIssue]:
    """
    Check for chemically unreasonable structures.

    Args:
        mol: RDKit molecule object

    Returns:
        List of ValidationIssue objects for unreasonable structures
    """
    if not RDKIT_AVAILABLE:
        return []
    if mol is None:
        return []

    issues = []

    # Check for unreasonable substructures
    for pattern_name, smarts in UNREASONABLE_PATTERNS.items():
        pattern = Chem.MolFromSmarts(smarts)
        if pattern is None:
            continue

        matches = mol.GetSubstructMatches(pattern)
        if matches:
            for match in matches:
                issues.append(ValidationIssue(
                    validation_type=ValidationType.CHEMICAL_REASONABILITY,
                    severity=ValidationSeverity.ERROR,
                    message=f"Unreasonable structure detected: {pattern_name}",
                    atom_indices=list(match),
                    details={
                        "pattern_name": pattern_name,
                        "smarts": smarts,
                    },
                ))

    # Check molecular weight reasonability (warn if very large)
    mw = Descriptors.MolWt(mol)
    if mw > 2000:
        issues.append(ValidationIssue(
            validation_type=ValidationType.CHEMICAL_REASONABILITY,
            severity=ValidationSeverity.INFO,
            message=f"Large molecular weight: {mw:.1f} Da",
            details={"molecular_weight": mw},
        ))

    # Check for too many radicals
    radical_count = 0
    for atom in mol.GetAtoms():
        radical_count += atom.GetNumRadicalElectrons()

    if radical_count > 2:
        issues.append(ValidationIssue(
            validation_type=ValidationType.CHEMICAL_REASONABILITY,
            severity=ValidationSeverity.WARNING,
            message=f"Multiple radical centers: {radical_count} unpaired electrons",
            details={"radical_count": radical_count},
        ))

    return issues


def check_sanitization(mol) -> List[ValidationIssue]:
    """
    Attempt to sanitize molecule and report issues.

    Args:
        mol: RDKit molecule object

    Returns:
        List of ValidationIssue objects for sanitization problems
    """
    if not RDKIT_AVAILABLE:
        return []
    if mol is None:
        return [ValidationIssue(
            validation_type=ValidationType.SANITIZATION,
            severity=ValidationSeverity.CRITICAL,
            message="No molecule provided for sanitization",
        )]

    issues = []

    # Make a copy for sanitization testing
    mol_copy = Chem.RWMol(mol)

    # Try sanitization with individual flags to identify specific problems
    sanitize_flags = [
        (Chem.SanitizeFlags.SANITIZE_FINDRADICALS, "radical detection"),
        (Chem.SanitizeFlags.SANITIZE_KEKULIZE, "kekulization"),
        (Chem.SanitizeFlags.SANITIZE_SETAROMATICITY, "aromaticity"),
        (Chem.SanitizeFlags.SANITIZE_SETCONJUGATION, "conjugation"),
        (Chem.SanitizeFlags.SANITIZE_SETHYBRIDIZATION, "hybridization"),
        (Chem.SanitizeFlags.SANITIZE_CLEANUP, "cleanup"),
        (Chem.SanitizeFlags.SANITIZE_PROPERTIES, "properties"),
    ]

    for flag, flag_name in sanitize_flags:
        mol_test = Chem.RWMol(mol)
        result = Chem.SanitizeMol(mol_test, sanitizeOps=flag, catchErrors=True)

        if result != Chem.SanitizeFlags.SANITIZE_NONE:
            issues.append(ValidationIssue(
                validation_type=ValidationType.SANITIZATION,
                severity=ValidationSeverity.ERROR,
                message=f"Sanitization failed at step: {flag_name}",
                details={
                    "sanitize_flag": flag_name,
                    "error_code": int(result),
                },
            ))

    # Try full sanitization
    full_result = Chem.SanitizeMol(mol_copy, catchErrors=True)

    if full_result != Chem.SanitizeFlags.SANITIZE_NONE:
        issues.append(ValidationIssue(
            validation_type=ValidationType.SANITIZATION,
            severity=ValidationSeverity.CRITICAL,
            message=f"Full sanitization failed with code {int(full_result)}",
            details={"error_code": int(full_result)},
        ))

    return issues


# Validator registry
VALIDATORS: Dict[ValidationType, Callable] = {
    ValidationType.VALENCE: check_valence,
    ValidationType.AROMATICITY: check_aromaticity,
    ValidationType.RING_STRAIN: check_ring_strain,
    ValidationType.STEREOCHEMISTRY: check_stereochemistry,
    ValidationType.STERIC_CLASH: check_steric_clashes,
    ValidationType.CONNECTIVITY: check_connectivity,
    ValidationType.FORMAL_CHARGE: check_formal_charges,
    ValidationType.CHEMICAL_REASONABILITY: check_chemical_reasonability,
    ValidationType.SANITIZATION: check_sanitization,
}


class ValidationPipeline:
    """
    Configurable validation pipeline for molecular structures.

    The pipeline runs multiple validators in sequence and aggregates results.
    Validators can be enabled/disabled and configured individually.

    Usage:
        pipeline = ValidationPipeline()
        pipeline.disable(ValidationType.STERIC_CLASH)  # Skip expensive check
        result = pipeline.validate(mol)
    """

    def __init__(self):
        """Initialize with all validators enabled."""
        self.enabled_validators: Set[ValidationType] = set(VALIDATORS.keys())
        self.validator_config: Dict[ValidationType, Dict[str, Any]] = {}

    def enable(self, validation_type: ValidationType) -> "ValidationPipeline":
        """Enable a validator."""
        self.enabled_validators.add(validation_type)
        return self

    def disable(self, validation_type: ValidationType) -> "ValidationPipeline":
        """Disable a validator."""
        self.enabled_validators.discard(validation_type)
        return self

    def enable_only(self, *validation_types: ValidationType) -> "ValidationPipeline":
        """Enable only the specified validators."""
        self.enabled_validators = set(validation_types)
        return self

    def configure(self, validation_type: ValidationType, **kwargs) -> "ValidationPipeline":
        """Configure a specific validator."""
        self.validator_config[validation_type] = kwargs
        return self

    def validate(self, mol) -> ValidationResult:
        """
        Run validation pipeline on molecule.

        Args:
            mol: RDKit molecule object

        Returns:
            ValidationResult with all issues found
        """
        result = ValidationResult(valid=True)

        # Add molecule info
        if mol is not None and RDKIT_AVAILABLE:
            result.molecule_info = {
                "num_atoms": mol.GetNumAtoms(),
                "num_bonds": mol.GetNumBonds(),
                "num_rings": mol.GetRingInfo().NumRings(),
                "smiles": Chem.MolToSmiles(mol),
            }

        # Run each enabled validator
        for validation_type in self.enabled_validators:
            validator = VALIDATORS.get(validation_type)
            if validator is None:
                continue

            result.checks_performed.append(validation_type.value)

            # Run validator
            issues = validator(mol)

            # Add issues to result
            for issue in issues:
                result.add_issue(issue)

        return result

    def quick_validate(self, mol) -> ValidationResult:
        """
        Run quick validation (sanitization and valence only).

        Args:
            mol: RDKit molecule object

        Returns:
            ValidationResult with critical issues
        """
        saved_validators = self.enabled_validators.copy()
        self.enable_only(ValidationType.SANITIZATION, ValidationType.VALENCE)
        result = self.validate(mol)
        self.enabled_validators = saved_validators
        return result

    def full_validate(self, mol) -> ValidationResult:
        """
        Run full validation with all checks.

        Args:
            mol: RDKit molecule object

        Returns:
            ValidationResult with all issues
        """
        saved_validators = self.enabled_validators.copy()
        self.enabled_validators = set(VALIDATORS.keys())
        result = self.validate(mol)
        self.enabled_validators = saved_validators
        return result


def validate_molecule(mol, checks: Optional[List[ValidationType]] = None) -> ValidationResult:
    """
    Convenience function to validate a molecule.

    Args:
        mol: RDKit molecule object or SMILES string
        checks: Optional list of specific checks to run (default: all)

    Returns:
        ValidationResult with all issues found
    """
    if not RDKIT_AVAILABLE:
        return ValidationResult(
            valid=False,
            issues=[ValidationIssue(
                validation_type=ValidationType.SANITIZATION,
                severity=ValidationSeverity.CRITICAL,
                message="RDKit not available for validation",
            )],
        )

    # Handle SMILES input
    if isinstance(mol, str):
        mol = Chem.MolFromSmiles(mol)
        if mol is None:
            return ValidationResult(
                valid=False,
                issues=[ValidationIssue(
                    validation_type=ValidationType.SANITIZATION,
                    severity=ValidationSeverity.CRITICAL,
                    message="Invalid SMILES string",
                )],
            )

    pipeline = ValidationPipeline()

    if checks is not None:
        pipeline.enable_only(*checks)

    return pipeline.validate(mol)


def validate_smiles(smiles: str) -> ValidationResult:
    """
    Validate a molecule from SMILES string.

    Args:
        smiles: SMILES string

    Returns:
        ValidationResult with all issues found
    """
    return validate_molecule(smiles)
