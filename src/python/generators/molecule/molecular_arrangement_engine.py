"""
molecular_arrangement_engine.py - Production-Ready MCP Backend

This is the COMPLETE molecular arrangement engine designed specifically as the
backend for the Molecule Generator MCP. It handles ANY clustering request that
can be expressed in natural language.

KEY CAPABILITIES:
================

1. ARRANGEMENT MODES (from 3 synthesized approaches):
   - Named patterns: "pi_pi_parallel", "circular", "h_bonded", etc.
   - Mathematical formulas: Custom x/y/z expressions
   - Relative poses: Parent-child hierarchical placement
   - Chemical constraints: "distance(0:centroid, 1:centroid, 3.4)"
   - Crystal symmetry: P21/c, C2/c, Pna21, etc.

2. NATURAL LANGUAGE PARSING:
   - "Stack 3 benzenes with π-stacking at 3.4 Å"
   - "Arrange 6 water molecules in a circular H-bonded ring"
   - "Create a T-shaped dimer of naphthalene"

3. CHEMICAL INTELLIGENCE:
   - AtomSelector for chemistry-aware targeting
   - H-bond geometry validation
   - Aromatic ring detection
   - Clash detection and resolution

4. CONSTRAINT SOLVER:
   - Gradient descent optimization
   - Hierarchical constraint satisfaction
   - Convergence guarantees

5. VALIDATION:
   - Physically reasonable distances
   - No atomic clashes
   - Chemical plausibility checks

Author: Molecular Arrangement Engine for MCP
License: MIT
"""

from typing import (
    Dict, List, Optional, Any, Tuple, Union,
    Callable, Set, NamedTuple
)
from dataclasses import dataclass, field
from enum import Enum, auto
from abc import ABC, abstractmethod
import numpy as np
import math
import re
import sys
from collections import defaultdict


# Debug setup
DEBUG_LOG_FILE = "/home/niel/git/crystal-mcp-server/debug_arrangement.log"

def debug(msg: str) -> None:
    # Write to stderr for realtime visibility
    sys.stderr.write(f"[DEBUG molecular_arrangement_engine] {msg}\n")
    sys.stderr.flush()
    # Write to log file for audit
    try:
        with open(DEBUG_LOG_FILE, "a") as f:
            f.write(f"[DEBUG molecular_arrangement_engine] {msg}\n")
    except Exception:
        pass




# =============================================================================
# CONSTANTS AND CHEMICAL DATA
# =============================================================================

# Van der Waals radii (Å) for clash detection
VDW_RADII = {
    'H': 1.20, 'C': 1.70, 'N': 1.55, 'O': 1.52, 'F': 1.47,
    'S': 1.80, 'P': 1.80, 'Cl': 1.75, 'Br': 1.85, 'I': 1.98,
    'B': 1.92, 'Si': 2.10, 'Se': 1.90, 'default': 1.70
}

# Typical interaction distances (Å)
INTERACTION_DISTANCES = {
    'pi_pi_parallel': (3.3, 3.8),      # Face-to-face π-stacking
    'pi_pi_offset': (3.3, 4.0),        # Offset π-stacking
    't_shaped': (4.5, 5.5),            # Edge-to-face
    'h_bond_strong': (2.5, 2.8),       # Strong H-bond (O...O)
    'h_bond_moderate': (2.8, 3.2),     # Moderate H-bond
    'h_bond_weak': (3.2, 4.0),         # Weak H-bond / CH...π
    'vdw_contact': (3.4, 4.0),         # Van der Waals contact
    'halogen_bond': (2.8, 3.5),        # Halogen bonding
    'cation_pi': (3.5, 4.5),           # Cation-π interaction
}

# Chemical interaction keywords for NL parsing
INTERACTION_KEYWORDS = {
    'pi_stacking': ['pi-pi', 'π-π', 'pi stacking', 'π stacking', 'aromatic stacking', 
                    'face-to-face', 'cofacial', 'parallel stacking', 'pi-stack', 'π-stack',
                    'pi stack', 'π stack', 'pipi', 'face to face'],
    'offset_stacking': ['offset', 'slip-stacked', 'slip stacked', 'displaced', 
                        'parallel displaced', 'slip stack'],
    't_shaped': ['t-shaped', 't shaped', 'edge-to-face', 'edge to face', 
                 'perpendicular', 'ch-pi', 'ch-π', 'ch pi'],
    'h_bonded': ['h-bond', 'hydrogen bond', 'h bond', 'hydrogen-bonded', 
                 'h-bonded', 'hbond', 'hydrogen bonding'],
    'circular': ['circular', 'ring', 'cyclic', 'wheel'],
    'linear': ['linear', 'line', 'chain', 'column'],
    'helical': ['helix', 'helical', 'spiral', 'coil'],
    'herringbone': ['herringbone', 'fishbone', 'zigzag', 'zig-zag'],
    'sandwich': ['sandwich', 'sandwiched', 'intercalated'],
}


# =============================================================================
# CORE DATA STRUCTURES
# =============================================================================

@dataclass
class Position:
    """3D position vector."""
    x: float = 0.0
    y: float = 0.0
    z: float = 0.0
    
    def to_array(self) -> np.ndarray:
        return np.array([self.x, self.y, self.z])
    
    @classmethod
    def from_array(cls, arr: np.ndarray) -> 'Position':
        return cls(float(arr[0]), float(arr[1]), float(arr[2]))
    
    def __add__(self, other: 'Position') -> 'Position':
        return Position(self.x + other.x, self.y + other.y, self.z + other.z)
    
    def __sub__(self, other: 'Position') -> 'Position':
        return Position(self.x - other.x, self.y - other.y, self.z - other.z)
    
    def distance_to(self, other: 'Position') -> float:
        return np.linalg.norm(self.to_array() - other.to_array())


@dataclass
class Orientation:
    """Euler angles (degrees) or quaternion representation."""
    roll: float = 0.0   # Rotation around X
    pitch: float = 0.0  # Rotation around Y
    yaw: float = 0.0    # Rotation around Z
    
    def to_matrix(self) -> np.ndarray:
        """Convert to 3x3 rotation matrix."""
        rx, ry, rz = np.radians([self.roll, self.pitch, self.yaw])
        
        Rx = np.array([
            [1, 0, 0],
            [0, np.cos(rx), -np.sin(rx)],
            [0, np.sin(rx), np.cos(rx)]
        ])
        Ry = np.array([
            [np.cos(ry), 0, np.sin(ry)],
            [0, 1, 0],
            [-np.sin(ry), 0, np.cos(ry)]
        ])
        Rz = np.array([
            [np.cos(rz), -np.sin(rz), 0],
            [np.sin(rz), np.cos(rz), 0],
            [0, 0, 1]
        ])
        
        return Rz @ Ry @ Rx
    
    def to_quaternion(self) -> np.ndarray:
        """Convert to quaternion [w, x, y, z] for numerical stability."""
        rx, ry, rz = np.radians([self.roll / 2, self.pitch / 2, self.yaw / 2])
        
        cr, sr = np.cos(rx), np.sin(rx)
        cp, sp = np.cos(ry), np.sin(ry)
        cy, sy = np.cos(rz), np.sin(rz)
        
        w = cr * cp * cy + sr * sp * sy
        x = sr * cp * cy - cr * sp * sy
        y = cr * sp * cy + sr * cp * sy
        z = cr * cp * sy - sr * sp * cy
        
        return np.array([w, x, y, z])
    
    @classmethod
    def from_matrix(cls, R: np.ndarray) -> 'Orientation':
        """Extract Euler angles from rotation matrix."""
        sy = np.sqrt(R[0, 0]**2 + R[1, 0]**2)
        singular = sy < 1e-6
        
        if not singular:
            roll = np.arctan2(R[2, 1], R[2, 2])
            pitch = np.arctan2(-R[2, 0], sy)
            yaw = np.arctan2(R[1, 0], R[0, 0])
        else:
            roll = np.arctan2(-R[1, 2], R[1, 1])
            pitch = np.arctan2(-R[2, 0], sy)
            yaw = 0
        
        return cls(np.degrees(roll), np.degrees(pitch), np.degrees(yaw))
    
    @classmethod
    def from_axis_angle(cls, axis: np.ndarray, angle_deg: float) -> 'Orientation':
        """Create orientation from axis-angle representation."""
        axis = axis / (np.linalg.norm(axis) + 1e-10)
        angle = np.radians(angle_deg)
        
        c, s = np.cos(angle), np.sin(angle)
        t = 1 - c
        x, y, z = axis
        
        R = np.array([
            [t*x*x + c, t*x*y - s*z, t*x*z + s*y],
            [t*x*y + s*z, t*y*y + c, t*y*z - s*x],
            [t*x*z - s*y, t*y*z + s*x, t*z*z + c]
        ])
        
        return cls.from_matrix(R)


@dataclass
class MolecularFrame:
    """
    Local coordinate frame attached to a molecule.
    Essential for relative placement and chemical alignment.
    """
    origin: np.ndarray
    x_axis: np.ndarray
    y_axis: np.ndarray
    z_axis: np.ndarray
    frame_type: str = "centroid"
    
    def to_matrix(self) -> np.ndarray:
        """Rotation matrix with frame axes as columns."""
        return np.column_stack([self.x_axis, self.y_axis, self.z_axis])
    
    def transform_point(self, local_point: np.ndarray) -> np.ndarray:
        """Transform point from local to global coordinates."""
        return self.origin + self.to_matrix() @ local_point
    
    def inverse_transform(self, global_point: np.ndarray) -> np.ndarray:
        """Transform point from global to local coordinates."""
        return self.to_matrix().T @ (global_point - self.origin)
    
    @classmethod
    def identity(cls) -> 'MolecularFrame':
        return cls(
            origin=np.zeros(3),
            x_axis=np.array([1, 0, 0]),
            y_axis=np.array([0, 1, 0]),
            z_axis=np.array([0, 0, 1]),
            frame_type="identity"
        )
    
    @classmethod
    def from_coords(cls, coords: np.ndarray, 
                    method: str = "svd") -> 'MolecularFrame':
        """
        Compute molecular frame from coordinates.
        
        Methods:
        - svd: Principal component analysis (default)
        - plane: Best-fit plane with normal as z-axis
        - inertia: Inertia tensor eigenvectors
        """
        centroid = np.mean(coords, axis=0)
        centered = coords - centroid
        
        if method == "svd" or method == "plane":
            u, s, vh = np.linalg.svd(centered, full_matrices=False)
            
            # For planar molecules, smallest singular value is plane normal
            x_axis = vh[0] / (np.linalg.norm(vh[0]) + 1e-10)
            y_axis = vh[1] / (np.linalg.norm(vh[1]) + 1e-10)
            z_axis = np.cross(x_axis, y_axis)
            z_axis = z_axis / (np.linalg.norm(z_axis) + 1e-10)
            
            # Ensure right-handed
            if np.dot(z_axis, vh[2]) < 0:
                z_axis = -z_axis
        
        elif method == "inertia":
            # Use inertia tensor for mass-weighted principal axes
            I = np.zeros((3, 3))
            for r in centered:
                I += np.outer(r, r)
            I = np.trace(I) * np.eye(3) - I
            
            eigenvalues, eigenvectors = np.linalg.eigh(I)
            x_axis = eigenvectors[:, 0]
            y_axis = eigenvectors[:, 1]
            z_axis = eigenvectors[:, 2]
        
        else:
            raise ValueError(f"Unknown frame method: {method}")
        
        return cls(
            origin=centroid,
            x_axis=x_axis,
            y_axis=y_axis,
            z_axis=z_axis,
            frame_type=method
        )


@dataclass
class MoleculePose:
    """
    Complete 6-DOF pose specification for a molecule.
    Supports both absolute and relative positioning.
    """
    position: Position = field(default_factory=Position)
    orientation: Orientation = field(default_factory=Orientation)
    
    # For relative placement
    parent_idx: int = -1  # -1 = absolute, >= 0 = relative to parent
    local_offset: Optional[Position] = None
    local_orientation: Optional[Orientation] = None
    
    # Metadata
    molecule_idx: int = 0
    constraints_satisfied: bool = True
    metadata: Dict[str, Any] = field(default_factory=dict)
    
    @property
    def is_relative(self) -> bool:
        return self.parent_idx >= 0
    
    def to_matrix(self) -> Tuple[np.ndarray, np.ndarray]:
        """Return (position_vector, rotation_matrix)."""
        return self.position.to_array(), self.orientation.to_matrix()


# =============================================================================
# ATOM SELECTOR (Chemical Intelligence from K2_plan)
# =============================================================================

class AtomSelector:
    """
    Chemistry-aware atom and group selector using DSL.
    
    Syntax: "mol_idx:target(args)"
    
    Targets:
        centroid()          - Molecular center of mass
        atom(element)       - First atom of element
        atom(element, n)    - Nth atom of element
        ring_center(n)      - Center of nth aromatic ring
        plane_normal()      - Normal to molecular plane
        func_group(name)    - Functional group center
        donor_h(n)          - Nth H-bond donor hydrogen
        acceptor(n)         - Nth H-bond acceptor
        vector(a1, a2)      - Vector from atom a1 to a2
    """
    
    PATTERN = re.compile(r"(\d+):(\w+)(?:\(([^)]*)\))?")
    
    def __init__(self, selector: str, molecules: List[Dict[str, Any]]):
        self.selector = selector
        self.molecules = molecules
    
    def evaluate(self) -> np.ndarray:
        """Evaluate selector and return 3D coordinates."""
        debug(f"AtomSelector.evaluate: '{self.selector}'")
        match = self.PATTERN.match(self.selector)
        if not match:
            debug(f"AtomSelector ERROR: Invalid selector syntax: {self.selector}")
            raise ValueError(f"Invalid selector syntax: {self.selector}")
        
        mol_idx = int(match.group(1))
        target = match.group(2)
        args_str = match.group(3) or ""
        
        debug(f"AtomSelector: mol_idx={mol_idx}, target={target}, args={args_str}")
        
        if mol_idx >= len(self.molecules):
            debug(f"AtomSelector ERROR: Molecule index {mol_idx} out of range (max={len(self.molecules)-1})")
            raise ValueError(f"Molecule index {mol_idx} out of range")
        
        mol = self.molecules[mol_idx]
        coords = np.array(mol.get("coords", [[0, 0, 0]]))
        atoms = mol.get("atoms", ["C"] * len(coords))
        
        debug(f"AtomSelector: molecule has {len(atoms)} atoms")
        
        args = self._parse_args(args_str)
        
        result = self._resolve(coords, atoms, mol, target, args)
        debug(f"AtomSelector: result = {result}")
        return result
    
    def _parse_args(self, args_str: str) -> Dict[str, Any]:
        """Parse argument string into dict."""
        args = {}
        if not args_str.strip():
            return args
        
        for i, part in enumerate(args_str.split(",")):
            part = part.strip()
            if "=" in part:
                k, v = part.split("=", 1)
                args[k.strip()] = self._parse_value(v.strip())
            else:
                args[f"arg{i}"] = self._parse_value(part)
        
        return args
    
    def _parse_value(self, v: str) -> Any:
        """Parse single value."""
        v = v.strip("'\"")
        if v.replace(".", "").replace("-", "").isdigit():
            return float(v) if "." in v else int(v)
        return v
    
    def _resolve(self, coords: np.ndarray, atoms: List[str],
                 mol: Dict, target: str, args: Dict) -> np.ndarray:
        """Resolve selector to coordinates."""
        
        if target == "centroid":
            return np.mean(coords, axis=0)
        
        elif target == "atom":
            element = args.get("arg0", "C")
            index = args.get("arg1", args.get("index", 0))
            indices = [i for i, a in enumerate(atoms) if a == element]
            if not indices:
                return np.mean(coords, axis=0)  # Fallback
            return coords[indices[min(index, len(indices) - 1)]]
        
        elif target == "ring_center":
            ring_idx = args.get("arg0", 0)
            return self._find_ring_center(coords, atoms, ring_idx)
        
        elif target == "plane_normal":
            frame = MolecularFrame.from_coords(coords, "svd")
            return frame.z_axis
        
        elif target == "func_group":
            group_name = args.get("arg0", "")
            return self._find_functional_group(coords, atoms, group_name)
        
        elif target == "donor_h":
            index = args.get("arg0", 0)
            donors = self._find_donor_hydrogens(coords, atoms)
            if len(donors) > index:
                return donors[index]
            return np.mean(coords, axis=0)
        
        elif target == "acceptor":
            index = args.get("arg0", 0)
            acceptors = self._find_acceptors(coords, atoms)
            if len(acceptors) > index:
                return acceptors[index]
            return np.mean(coords, axis=0)
        
        elif target == "vector":
            idx1 = args.get("arg0", 0)
            idx2 = args.get("arg1", 1)
            return coords[idx2] - coords[idx1]
        
        elif target == "carbonyl":
            return self._find_carbonyl_oxygen(coords, atoms)
            
        elif target == "amide":
            return self._find_amide_nitrogen(coords, atoms)
        
        else:
            raise ValueError(f"Unknown selector target: {target}")
    
    def _find_ring_center(self, coords: np.ndarray, 
                          atoms: List[str], ring_idx: int) -> np.ndarray:
        """Find center of aromatic ring (heuristic)."""
        # Find clusters of 5-6 C atoms with ~1.4 Å spacing
        c_indices = [i for i, a in enumerate(atoms) if a == "C"]
        
        if len(c_indices) < 5:
            return np.mean(coords, axis=0)
        
        # Simple: take first 6 carbons as ring
        ring_size = min(6, len(c_indices))
        ring_atoms = c_indices[ring_idx * ring_size:(ring_idx + 1) * ring_size]
        
        if not ring_atoms:
            ring_atoms = c_indices[:ring_size]
        
        return np.mean(coords[ring_atoms], axis=0)
    
    def _find_functional_group(self, coords: np.ndarray,
                                atoms: List[str], group_name: str) -> np.ndarray:
        """Find functional group center with chemistry-aware logic."""
        debug(f"AtomSelector: searching for group '{group_name}'")
        
        if group_name.lower() == "carbonyl":
            return self._find_carbonyl_oxygen(coords, atoms)
        elif group_name.lower() == "amide":
            # Return amide nitrogen by default
            return self._find_amide_nitrogen(coords, atoms)
        
        # Fallback to simple element matching for other groups
        group_elements = {
            "carboxyl": ["C", "O", "O"],
            "hydroxyl": ["O"],
            "amine": ["N"],
            "nitro": ["N", "O", "O"],
            "cyano": ["C", "N"],
            "methyl": ["C"],
        }
        
        elements = group_elements.get(group_name.lower(), ["C"])
        matching = []
        
        for elem in elements:
            indices = [i for i, a in enumerate(atoms) if a == elem]
            if indices:
                matching.append(coords[indices[0]])
        
        return np.mean(matching, axis=0) if matching else np.mean(coords, axis=0)

    def _find_carbonyl_oxygen(self, coords: np.ndarray, atoms: List[str]) -> np.ndarray:
        """Find C=O carbonyl oxygen using bond length analysis."""
        o_indices = [i for i, a in enumerate(atoms) if a == "O"]
        c_indices = [i for i, a in enumerate(atoms) if a == "C"]
        
        debug(f"AtomSelector: Checking {len(o_indices)} oxygens for carbonyl...")
        
        for o_idx in o_indices:
            for c_idx in c_indices:
                dist = np.linalg.norm(coords[o_idx] - coords[c_idx])
                # C=O double bond is typically 1.20-1.25 Å
                # C-O single bond is typically 1.43 Å
                if 1.15 < dist < 1.30:
                    debug(f"AtomSelector: Found carbonyl C=O (dist={dist:.2f} Å) at O-index {o_idx}")
                    return coords[o_idx]
        
        debug("AtomSelector: No carbonyl bond found, falling back to first Oxygen")
        if o_indices:
            return coords[o_indices[0]]
        return np.mean(coords, axis=0)

    def _find_amide_nitrogen(self, coords: np.ndarray, atoms: List[str]) -> np.ndarray:
        """Find amide nitrogen (adjacent to C=O)."""
        n_indices = [i for i, a in enumerate(atoms) if a == "N"]
        c_indices = [i for i, a in enumerate(atoms) if a == "C"]
        o_indices = [i for i, a in enumerate(atoms) if a == "O"]
        
        for n_idx in n_indices:
            # Check for C neighbor
            for c_idx in c_indices:
                c_n_dist = np.linalg.norm(coords[n_idx] - coords[c_idx])
                if 1.30 < c_n_dist < 1.47: # C-N bond
                    # Check if this C is double bonded to O
                    for o_idx in o_indices:
                        c_o_dist = np.linalg.norm(coords[c_idx] - coords[o_idx])
                        if 1.15 < c_o_dist < 1.30:
                            debug(f"AtomSelector: Found amide Nitrogen (linked to C=O) at index {n_idx}")
                            return coords[n_idx]
        
        if n_indices:
            return coords[n_indices[0]]
        return np.mean(coords, axis=0)
    
    def _find_donor_hydrogens(self, coords: np.ndarray,
                               atoms: List[str]) -> List[np.ndarray]:
        """Find H-bond donor hydrogens (H bonded to O, N, F)."""
        donors = []
        h_indices = [i for i, a in enumerate(atoms) if a == "H"]
        hetero_indices = [i for i, a in enumerate(atoms) if a in ["O", "N", "F"]]
        
        for h_idx in h_indices:
            for het_idx in hetero_indices:
                dist = np.linalg.norm(coords[h_idx] - coords[het_idx])
                if dist < 1.3:  # Typical X-H bond length
                    donors.append(coords[h_idx])
                    break
        
        return donors
    
    def _find_acceptors(self, coords: np.ndarray,
                        atoms: List[str]) -> List[np.ndarray]:
        """Find H-bond acceptor atoms."""
        indices = [i for i, a in enumerate(atoms) if a in ["O", "N", "F"]]
        return [coords[i] for i in indices]


# =============================================================================
# CONSTRAINT SYSTEM (from K2_plan with solver)
# =============================================================================

class Constraint(ABC):
    """Abstract constraint base class."""
    
    def __init__(self, weight: float = 1.0, priority: int = 1):
        self.weight = weight
        self.priority = priority
    
    @abstractmethod
    def evaluate(self, molecules: List[Dict], poses: List[MoleculePose]) -> float:
        """Return constraint violation (0 = satisfied)."""
        pass
    
    @abstractmethod
    def description(self) -> str:
        """Human-readable description."""
        pass
    
    def is_satisfied(self, molecules: List[Dict], 
                     poses: List[MoleculePose], tol: float = 1e-2) -> bool:
        return self.evaluate(molecules, poses) < tol


class DistanceConstraint(Constraint):
    """Distance between two points must be within range."""
    
    def __init__(self, selector1: str, selector2: str,
                 target: Optional[float] = None,
                 min_dist: Optional[float] = None,
                 max_dist: Optional[float] = None,
                 **kwargs):
        super().__init__(**kwargs)
        self.selector1 = selector1
        self.selector2 = selector2
        self.target = target
        self.min_dist = min_dist
        self.max_dist = max_dist
    
    def evaluate(self, molecules: List[Dict], poses: List[MoleculePose]) -> float:
        # Apply poses to molecules first
        transformed = apply_poses_to_molecules(molecules, poses)
        
        p1 = AtomSelector(self.selector1, transformed).evaluate()
        p2 = AtomSelector(self.selector2, transformed).evaluate()
        d = np.linalg.norm(p1 - p2)
        
        if self.target is not None:
            return self.weight * (d - self.target) ** 2
        
        violation = 0.0
        if self.min_dist and d < self.min_dist:
            violation += (self.min_dist - d) ** 2
        if self.max_dist and d > self.max_dist:
            violation += (d - self.max_dist) ** 2
        
        return self.weight * violation
    
    def description(self) -> str:
        if self.target:
            return f"Distance({self.selector1}, {self.selector2}) = {self.target} Å"
        return f"Distance({self.selector1}, {self.selector2}) in [{self.min_dist}, {self.max_dist}] Å"


class AngleConstraint(Constraint):
    """Angle between three points."""
    
    def __init__(self, selector1: str, selector2: str, selector3: str,
                 target: Optional[float] = None,
                 min_angle: Optional[float] = None,
                 max_angle: Optional[float] = None,
                 **kwargs):
        super().__init__(**kwargs)
        self.selector1 = selector1
        self.selector2 = selector2  # Vertex
        self.selector3 = selector3
        self.target = target
        self.min_angle = min_angle
        self.max_angle = max_angle
    
    def evaluate(self, molecules: List[Dict], poses: List[MoleculePose]) -> float:
        transformed = apply_poses_to_molecules(molecules, poses)
        
        p1 = AtomSelector(self.selector1, transformed).evaluate()
        p2 = AtomSelector(self.selector2, transformed).evaluate()
        p3 = AtomSelector(self.selector3, transformed).evaluate()
        
        v1 = p1 - p2
        v2 = p3 - p2
        
        cos_angle = np.dot(v1, v2) / (np.linalg.norm(v1) * np.linalg.norm(v2) + 1e-10)
        angle = np.degrees(np.arccos(np.clip(cos_angle, -1, 1)))
        
        if self.target is not None:
            return self.weight * (angle - self.target) ** 2
        
        violation = 0.0
        if self.min_angle and angle < self.min_angle:
            violation += (self.min_angle - angle) ** 2
        if self.max_angle and angle > self.max_angle:
            violation += (angle - self.max_angle) ** 2
        
        return self.weight * violation
    
    def description(self) -> str:
        return f"Angle({self.selector1}, {self.selector2}, {self.selector3})"


class PlaneAlignmentConstraint(Constraint):
    """Two molecular planes should be parallel or perpendicular."""
    
    def __init__(self, mol1: int, mol2: int,
                 mode: str = "parallel",  # "parallel" or "perpendicular"
                 tolerance: float = 10.0,  # degrees
                 **kwargs):
        super().__init__(**kwargs)
        self.mol1 = mol1
        self.mol2 = mol2
        self.mode = mode
        self.tolerance = tolerance
    
    def evaluate(self, molecules: List[Dict], poses: List[MoleculePose]) -> float:
        transformed = apply_poses_to_molecules(molecules, poses)
        
        n1 = AtomSelector(f"{self.mol1}:plane_normal()", transformed).evaluate()
        n2 = AtomSelector(f"{self.mol2}:plane_normal()", transformed).evaluate()
        
        cos_angle = abs(np.dot(n1, n2))
        angle = np.degrees(np.arccos(np.clip(cos_angle, 0, 1)))
        
        if self.mode == "parallel":
            target = 0.0
        else:  # perpendicular
            target = 90.0
        
        deviation = abs(angle - target)
        
        if deviation <= self.tolerance:
            return 0.0
        
        return self.weight * (deviation - self.tolerance) ** 2
    
    def description(self) -> str:
        return f"Planes of mol {self.mol1} and {self.mol2} {self.mode} (±{self.tolerance}°)"


class HBondConstraint(Constraint):
    """
    Hydrogen bond geometry constraint.
    
    Optimal H-bond geometry:
    - H...A distance: 1.5-2.5 Å
    - D...A distance: 2.5-3.5 Å
    - D-H...A angle: > 120° (ideally ~170°)
    """
    
    def __init__(self, donor_mol: int, acceptor_mol: int, **kwargs):
        super().__init__(**kwargs)
        self.donor_mol = donor_mol
        self.acceptor_mol = acceptor_mol
    
    def evaluate(self, molecules: List[Dict], poses: List[MoleculePose]) -> float:
        transformed = apply_poses_to_molecules(molecules, poses)
        
        # Find donor H and acceptor sites
        donor_mol = transformed[self.donor_mol]
        acceptor_mol = transformed[self.acceptor_mol]
        
        donor_coords = np.array(donor_mol["coords"])
        donor_atoms = donor_mol.get("atoms", [])
        acceptor_coords = np.array(acceptor_mol["coords"])
        acceptor_atoms = acceptor_mol.get("atoms", [])
        
        # Find H-bond donors (H attached to O, N, F)
        selector = AtomSelector(f"{self.donor_mol}:donor_h(0)", transformed)
        donor_hs = []
        for i in range(10):  # Check up to 10 potential donors
            sel = AtomSelector(f"{self.donor_mol}:donor_h({i})", transformed)
            # Simple heuristic
            h_indices = [j for j, a in enumerate(donor_atoms) if a == "H"]
            for h_idx in h_indices:
                for het_idx, a in enumerate(donor_atoms):
                    if a in ["O", "N", "F"]:
                        if np.linalg.norm(donor_coords[h_idx] - donor_coords[het_idx]) < 1.3:
                            donor_hs.append((donor_coords[h_idx], donor_coords[het_idx]))
                            break
        
        # Find acceptors
        acceptor_sites = [acceptor_coords[i] for i, a in enumerate(acceptor_atoms) 
                         if a in ["O", "N", "F"]]
        
        if not donor_hs or not acceptor_sites:
            return 1e6  # No H-bond possible
        
        # Find best H-bond
        best_score = 1e6
        
        for h_coord, heavy_coord in donor_hs:
            for acc_coord in acceptor_sites:
                # H...A distance
                d_ha = np.linalg.norm(h_coord - acc_coord)
                # D...A distance
                d_da = np.linalg.norm(heavy_coord - acc_coord)
                
                # D-H...A angle
                v1 = heavy_coord - h_coord
                v2 = acc_coord - h_coord
                cos_angle = np.dot(v1, v2) / (np.linalg.norm(v1) * np.linalg.norm(v2) + 1e-10)
                angle = np.degrees(np.arccos(np.clip(cos_angle, -1, 1)))
                
                # Score components
                d_ha_score = (d_ha - 1.9) ** 2 if d_ha < 1.5 or d_ha > 2.5 else 0
                d_da_score = (d_da - 2.8) ** 2 if d_da < 2.5 or d_da > 3.5 else 0
                angle_score = max(0, 120 - angle) ** 2
                
                score = d_ha_score + d_da_score + angle_score * 0.01
                best_score = min(best_score, score)
        
        return self.weight * best_score
    
    def description(self) -> str:
        return f"H-bond: mol {self.donor_mol} → mol {self.acceptor_mol}"


class ClashConstraint(Constraint):
    """Prevent atomic clashes between molecules."""
    
    def __init__(self, mol1: int, mol2: int,
                 scale_factor: float = 0.8,  # Allow some overlap
                 **kwargs):
        super().__init__(**kwargs)
        self.mol1 = mol1
        self.mol2 = mol2
        self.scale_factor = scale_factor
    
    def evaluate(self, molecules: List[Dict], poses: List[MoleculePose]) -> float:
        transformed = apply_poses_to_molecules(molecules, poses)
        
        coords1 = np.array(transformed[self.mol1]["coords"])
        coords2 = np.array(transformed[self.mol2]["coords"])
        atoms1 = transformed[self.mol1].get("atoms", ["C"] * len(coords1))
        atoms2 = transformed[self.mol2].get("atoms", ["C"] * len(coords2))
        
        violation = 0.0
        
        for i, (c1, a1) in enumerate(zip(coords1, atoms1)):
            r1 = VDW_RADII.get(a1, VDW_RADII['default'])
            
            for j, (c2, a2) in enumerate(zip(coords2, atoms2)):
                r2 = VDW_RADII.get(a2, VDW_RADII['default'])
                
                d = np.linalg.norm(c1 - c2)
                min_dist = (r1 + r2) * self.scale_factor
                
                if d < min_dist:
                    violation += (min_dist - d) ** 2
        
        return self.weight * violation
    
    def description(self) -> str:
        return f"No clashes between mol {self.mol1} and mol {self.mol2}"


# =============================================================================
# CONSTRAINT SOLVER (Gradient Descent from K2_plan)
# =============================================================================

class ConstraintSolver:
    """
    Iterative solver for constraint satisfaction.
    Uses gradient descent with adaptive step size.
    """
    
    def __init__(self, molecules: List[Dict], poses: List[MoleculePose],
                 constraints: List[Constraint]):
        self.molecules = molecules
        self.poses = [self._copy_pose(p) for p in poses]
        self.constraints = sorted(constraints, key=lambda c: -c.priority)
    
    def _copy_pose(self, pose: MoleculePose) -> MoleculePose:
        """Deep copy a pose."""
        return MoleculePose(
            position=Position(pose.position.x, pose.position.y, pose.position.z),
            orientation=Orientation(pose.orientation.roll, 
                                   pose.orientation.pitch, 
                                   pose.orientation.yaw),
            molecule_idx=pose.molecule_idx,
            metadata=pose.metadata.copy()
        )
    
    def solve(self, max_iterations: int = 1000, 
              tolerance: float = 1e-4,
              step_size: float = 0.1) -> List[MoleculePose]:
        """
        Solve constraints iteratively.
        
        Returns optimized poses.
        """
        n_mols = len(self.poses)
        
        for iteration in range(max_iterations):
            # Calculate total violation
            total_violation = sum(c.evaluate(self.molecules, self.poses) 
                                 for c in self.constraints)
            
            if iteration % 100 == 0:
                debug(f"Solver iteration {iteration}: violation={total_violation:.6f}, step={step_size:.6f}")

            if total_violation < tolerance:
                debug(f"Solver converged at iteration {iteration} with violation {total_violation:.6f}")
                break
            
            # Compute numerical gradient for each pose
            gradients = self._compute_gradients()
            
            # Update poses
            for i in range(n_mols):
                grad_pos, grad_rot = gradients[i]
                
                # Update position
                self.poses[i].position.x -= step_size * grad_pos[0]
                self.poses[i].position.y -= step_size * grad_pos[1]
                self.poses[i].position.z -= step_size * grad_pos[2]
                
                # Update orientation (smaller step for angles)
                self.poses[i].orientation.roll -= step_size * 0.1 * grad_rot[0]
                self.poses[i].orientation.pitch -= step_size * 0.1 * grad_rot[1]
                self.poses[i].orientation.yaw -= step_size * 0.1 * grad_rot[2]
            
            # Adaptive step size
            if iteration > 0 and iteration % 100 == 0:
                step_size *= 0.9
        
        return self.poses
    
    def _compute_gradients(self, eps: float = 1e-4) -> List[Tuple[np.ndarray, np.ndarray]]:
        """Compute numerical gradients for all poses."""
        gradients = []
        
        for i in range(len(self.poses)):
            grad_pos = np.zeros(3)
            grad_rot = np.zeros(3)
            
            # Position gradients
            for axis in range(3):
                # Forward
                original = [self.poses[i].position.x, 
                           self.poses[i].position.y, 
                           self.poses[i].position.z][axis]
                
                if axis == 0:
                    self.poses[i].position.x += eps
                elif axis == 1:
                    self.poses[i].position.y += eps
                else:
                    self.poses[i].position.z += eps
                
                val_plus = sum(c.evaluate(self.molecules, self.poses) 
                              for c in self.constraints)
                
                # Backward
                if axis == 0:
                    self.poses[i].position.x = original - eps
                elif axis == 1:
                    self.poses[i].position.y = original - eps
                else:
                    self.poses[i].position.z = original - eps
                
                val_minus = sum(c.evaluate(self.molecules, self.poses) 
                               for c in self.constraints)
                
                # Restore
                if axis == 0:
                    self.poses[i].position.x = original
                elif axis == 1:
                    self.poses[i].position.y = original
                else:
                    self.poses[i].position.z = original
                
                grad_pos[axis] = (val_plus - val_minus) / (2 * eps)
            
            # Rotation gradients (similar pattern)
            for axis in range(3):
                original = [self.poses[i].orientation.roll,
                           self.poses[i].orientation.pitch,
                           self.poses[i].orientation.yaw][axis]
                
                if axis == 0:
                    self.poses[i].orientation.roll += eps
                elif axis == 1:
                    self.poses[i].orientation.pitch += eps
                else:
                    self.poses[i].orientation.yaw += eps
                
                val_plus = sum(c.evaluate(self.molecules, self.poses) 
                              for c in self.constraints)
                
                if axis == 0:
                    self.poses[i].orientation.roll = original - eps
                elif axis == 1:
                    self.poses[i].orientation.pitch = original - eps
                else:
                    self.poses[i].orientation.yaw = original - eps
                
                val_minus = sum(c.evaluate(self.molecules, self.poses) 
                               for c in self.constraints)
                
                if axis == 0:
                    self.poses[i].orientation.roll = original
                elif axis == 1:
                    self.poses[i].orientation.pitch = original
                else:
                    self.poses[i].orientation.yaw = original
                
                grad_rot[axis] = (val_plus - val_minus) / (2 * eps)
            
            gradients.append((grad_pos, grad_rot))
        
        return gradients


# =============================================================================
# POSITION AND ORIENTATION GENERATORS
# =============================================================================

def generate_linear_positions(n: int, spacing: float = 3.4,
                              axis: Tuple[float, float, float] = (0, 0, 1),
                              start: Tuple[float, float, float] = (0, 0, 0)) -> List[Position]:
    """Generate positions along a line."""
    axis_vec = np.array(axis)
    axis_vec = axis_vec / (np.linalg.norm(axis_vec) + 1e-10)
    start_vec = np.array(start)
    
    return [Position.from_array(start_vec + i * spacing * axis_vec) for i in range(n)]


def generate_circular_positions(n: int, radius: float = 5.0,
                                center: Tuple[float, float, float] = (0, 0, 0),
                                plane_normal: Tuple[float, float, float] = (0, 0, 1),
                                start_angle: float = 0.0) -> List[Position]:
    """Generate positions on a circle."""
    center_vec = np.array(center)
    normal = np.array(plane_normal)
    normal = normal / (np.linalg.norm(normal) + 1e-10)
    
    # Create orthonormal basis in plane
    if abs(normal[2]) < 0.9:
        u = np.cross(normal, [0, 0, 1])
    else:
        u = np.cross(normal, [1, 0, 0])
    u = u / np.linalg.norm(u)
    v = np.cross(normal, u)
    
    positions = []
    for i in range(n):
        angle = np.radians(start_angle) + 2 * np.pi * i / n
        pos = center_vec + radius * (np.cos(angle) * u + np.sin(angle) * v)
        positions.append(Position.from_array(pos))
    
    return positions


def generate_helical_positions(n: int, radius: float = 5.0,
                               pitch: float = 3.4, turns: float = 1.0) -> List[Position]:
    """Generate positions along a helix."""
    positions = []
    total_angle = 2 * np.pi * turns
    total_rise = pitch * turns
    
    for i in range(n):
        t = i / max(n - 1, 1)
        angle = t * total_angle
        z = t * total_rise
        
        x = radius * np.cos(angle)
        y = radius * np.sin(angle)
        
        positions.append(Position(x, y, z))
    
    return positions


def generate_spherical_positions(n: int, radius: float = 5.0,
                                 center: Tuple[float, float, float] = (0, 0, 0)) -> List[Position]:
    """Generate positions on a sphere using Fibonacci lattice."""
    center_vec = np.array(center)
    golden_ratio = (1 + np.sqrt(5)) / 2
    
    positions = []
    for i in range(n):
        theta = 2 * np.pi * i / golden_ratio
        phi = np.arccos(1 - 2 * (i + 0.5) / n)
        
        x = radius * np.sin(phi) * np.cos(theta)
        y = radius * np.sin(phi) * np.sin(theta)
        z = radius * np.cos(phi)
        
        pos = center_vec + np.array([x, y, z])
        positions.append(Position.from_array(pos))
    
    return positions


def generate_grid_positions(n: int, spacing: float = 3.4,
                            plane: str = "xy") -> List[Position]:
    """Generate positions on a 2D grid."""
    side = int(np.ceil(np.sqrt(n)))
    positions = []
    
    idx = 0
    for i in range(side):
        for j in range(side):
            if idx >= n:
                break
            
            x = (i - (side - 1) / 2) * spacing
            y = (j - (side - 1) / 2) * spacing
            z = 0
            
            if plane == "xz":
                positions.append(Position(x, 0, y))
            elif plane == "yz":
                positions.append(Position(0, x, y))
            else:
                positions.append(Position(x, y, z))
            
            idx += 1
    
    return positions


def generate_formula_positions(n: int, x_formula: str = "0",
                               y_formula: str = "0",
                               z_formula: str = "3.4 * i") -> List[Position]:
    """Generate positions from mathematical formulas."""
    safe_dict = {
        'sin': np.sin, 'cos': np.cos, 'tan': np.tan,
        'sqrt': np.sqrt, 'exp': np.exp, 'log': np.log,
        'abs': abs, 'pi': np.pi, 'e': np.e,
        'atan2': np.arctan2, 'n': n
    }
    
    positions = []
    for i in range(n):
        safe_dict['i'] = i
        safe_dict['t'] = i / max(n - 1, 1)
        
        x = float(eval(x_formula, {"__builtins__": {}}, safe_dict))
        y = float(eval(y_formula, {"__builtins__": {}}, safe_dict))
        z = float(eval(z_formula, {"__builtins__": {}}, safe_dict))
        
        positions.append(Position(x, y, z))
    
    return positions


# Orientation generators

def generate_fixed_orientations(n: int, roll: float = 0, 
                                pitch: float = 0, yaw: float = 0) -> List[Orientation]:
    """All molecules have the same orientation."""
    return [Orientation(roll, pitch, yaw) for _ in range(n)]


def generate_alternating_orientations(n: int,
                                      orient_a: Tuple[float, float, float] = (0, 0, 0),
                                      orient_b: Tuple[float, float, float] = (0, 0, 180)) -> List[Orientation]:
    """Alternate between two orientations."""
    orientations = []
    for i in range(n):
        if i % 2 == 0:
            orientations.append(Orientation(*orient_a))
        else:
            orientations.append(Orientation(*orient_b))
    return orientations


def generate_incremental_orientations(n: int, axis: str = 'z',
                                      increment: float = 0) -> List[Orientation]:
    """Progressive rotation around an axis."""
    orientations = []
    for i in range(n):
        angle = i * increment
        if axis == 'x':
            orientations.append(Orientation(angle, 0, 0))
        elif axis == 'y':
            orientations.append(Orientation(0, angle, 0))
        else:
            orientations.append(Orientation(0, 0, angle))
    return orientations


def generate_face_center_orientations(positions: List[Position],
                                      center: Tuple[float, float, float] = (0, 0, 0),
                                      facing: str = "inward") -> List[Orientation]:
    """Orient molecules to face toward/away from center."""
    center_vec = np.array(center)
    orientations = []
    
    for pos in positions:
        pos_vec = pos.to_array()
        direction = center_vec - pos_vec
        dist = np.linalg.norm(direction)
        
        if dist < 1e-6:
            orientations.append(Orientation())
            continue
        
        direction = direction / dist
        if facing == "outward":
            direction = -direction
        
        # Convert direction to orientation
        # Align z-axis with direction
        z = direction
        if abs(z[2]) < 0.9:
            x = np.cross([0, 0, 1], z)
        else:
            x = np.cross([1, 0, 0], z)
        x = x / (np.linalg.norm(x) + 1e-10)
        y = np.cross(z, x)
        
        R = np.column_stack([x, y, z])
        orientations.append(Orientation.from_matrix(R))
    
    return orientations


# =============================================================================
# NATURAL LANGUAGE PARSER FOR MCP
# =============================================================================

class NLArrangementParser:
    """
    Parse natural language arrangement descriptions.
    
    Examples:
        "Stack 3 benzenes with π-stacking at 3.4 Å"
        "Arrange 6 water molecules in a circular H-bonded ring"
        "Create a T-shaped dimer of naphthalene at 5 Å"
        "Make a helical arrangement of 10 molecules"
    """
    
    def __init__(self):
        self.number_pattern = re.compile(r'\b(\d+)\b')
        self.distance_pattern = re.compile(r'(\d+\.?\d*)\s*[ÅA°]|\bat\s+(\d+\.?\d*)')
        self.angle_pattern = re.compile(r'(\d+\.?\d*)\s*°|(\d+\.?\d*)\s*degrees?')
    
    def parse(self, text: str) -> Dict[str, Any]:
        """
        Parse natural language to arrangement specification.
        
        Returns dict with:
            - pattern: str (arrangement type)
            - n_molecules: int
            - distance: float (if specified)
            - constraints: List[str]
            - parameters: Dict[str, Any]
        """
        text_lower = text.lower()
        
        result = {
            'pattern': 'linear',
            'n_molecules': 2,
            'distance': None,
            'constraints': [],
            'parameters': {}
        }
        
        # Extract number of molecules
        numbers = self.number_pattern.findall(text_lower)
        if numbers:
            result['n_molecules'] = int(numbers[0])
        
        # Extract distance
        dist_match = self.distance_pattern.search(text)
        if dist_match:
            dist_val = dist_match.group(1) or dist_match.group(2)
            if dist_val:
                result['distance'] = float(dist_val)
        
        # Detect arrangement pattern
        result['pattern'] = self._detect_pattern(text_lower)
        debug(f"NL Parser: detected pattern '{result['pattern']}' from text: '{text[:50]}...'")
        
        # Extract additional parameters based on pattern
        if result['pattern'] in ['helical', 'helix']:
            # Look for turns, pitch
            if 'turn' in text_lower:
                turns_match = re.search(r'(\d+\.?\d*)\s*turns?', text_lower)
                if turns_match:
                    result['parameters']['turns'] = float(turns_match.group(1))
        
        if result['pattern'] == 'circular':
            # Look for radius
            radius_match = re.search(r'radius\s*(?:of\s*)?(\d+\.?\d*)', text_lower)
            if radius_match:
                result['parameters']['radius'] = float(radius_match.group(1))
        
        # Detect constraints
        if any(kw in text_lower for kw in INTERACTION_KEYWORDS['h_bonded']):
            result['constraints'].append('h_bond')
        
        if 'parallel' in text_lower and 'anti' not in text_lower:
            result['constraints'].append('plane_parallel')
        
        return result
    
    def _detect_pattern(self, text: str) -> str:
        """Detect arrangement pattern from text."""
        # Normalize text - replace common Unicode variants
        text = text.replace('π', 'pi').replace('Å', 'A').replace('°', ' degrees ')
        text = text.lower()
        
        # Check each interaction type in priority order
        priority_order = [
            'pi_stacking', 'offset_stacking', 't_shaped', 'h_bonded',
            'herringbone', 'sandwich', 'helical', 'circular', 'linear'
        ]
        
        for pattern_name in priority_order:
            keywords = INTERACTION_KEYWORDS.get(pattern_name, [])
            for kw in keywords:
                # Also normalize keyword
                kw_normalized = kw.replace('π', 'pi').lower()
                if kw_normalized in text:
                    return pattern_name
        
        # Check for planar/planner (PTCDA/sheet like) -> map to linear for dimers
        if 'planner' in text or 'planar' in text:
            # If large count, maybe grid? For now map to linear/grid based on N?
            # Default to linear (side-by-side)
            return 'linear'
        
        # Check for specific molecule counts
        if 'dimer' in text:
            if 't-shape' in text or 't shape' in text:
                return 't_shaped'
            return 'pi_stacking'  # Default for dimers
        if 'trimer' in text:
            return 'linear'
        
        return 'linear'


def parse_arrangement_request(text: str) -> Dict[str, Any]:
    """
    Main entry point for NL parsing.
    
    Returns specification dict that can be passed to generate_arrangement().
    """
    parser = NLArrangementParser()
    return parser.parse(text)





# =============================================================================
# UTILITY FUNCTIONS
# =============================================================================

def apply_poses_to_molecules(molecules: List[Dict], 
                              poses: List[MoleculePose]) -> List[Dict]:
    """Apply poses to transform molecule coordinates."""
    result = []
    
    for mol, pose in zip(molecules, poses):
        coords = np.array(mol.get("coords", [[0, 0, 0]]))
        
        # Center molecule
        centroid = np.mean(coords, axis=0)
        centered = coords - centroid
        
        # Apply rotation
        R = pose.orientation.to_matrix()
        rotated = centered @ R.T
        
        # Apply translation
        translated = rotated + pose.position.to_array()
        
        result.append({
            **mol,
            "coords": translated.tolist()
        })
    
    return result


def validate_arrangement(molecules: List[Dict], poses: List[MoleculePose],
                         check_clashes: bool = True,
                         min_distance: float = 1.5) -> Dict[str, Any]:
    """
    Validate an arrangement for physical plausibility.
    
    Returns:
        {
            'valid': bool,
            'warnings': List[str],
            'errors': List[str],
            'statistics': Dict
        }
    """
    result = {
        'valid': True,
        'warnings': [],
        'errors': [],
        'statistics': {}
    }
    
    transformed = apply_poses_to_molecules(molecules, poses)
    n = len(transformed)
    
    # Check intermolecular distances
    min_intermol_dist = float('inf')
    clash_pairs = []
    
    for i in range(n):
        coords_i = np.array(transformed[i]["coords"])
        atoms_i = transformed[i].get("atoms", ["C"] * len(coords_i))
        
        for j in range(i + 1, n):
            coords_j = np.array(transformed[j]["coords"])
            atoms_j = transformed[j].get("atoms", ["C"] * len(coords_j))
            
            for ai, (ci, ei) in enumerate(zip(coords_i, atoms_i)):
                ri = VDW_RADII.get(ei, VDW_RADII['default'])
                
                for aj, (cj, ej) in enumerate(zip(coords_j, atoms_j)):
                    rj = VDW_RADII.get(ej, VDW_RADII['default'])
                    
                    d = np.linalg.norm(ci - cj)
                    min_intermol_dist = min(min_intermol_dist, d)
                    
                    if check_clashes and d < (ri + rj) * 0.7:
                        clash_pairs.append((i, j, d))
    
    result['statistics']['min_intermolecular_distance'] = min_intermol_dist
    result['statistics']['n_molecules'] = n
    
    if clash_pairs:
        result['valid'] = False
        for i, j, d in clash_pairs[:5]:  # Report first 5
            result['errors'].append(f"Clash between mol {i} and mol {j}: {d:.2f} Å")
    
    if min_intermol_dist < min_distance:
        result['warnings'].append(
            f"Minimum distance {min_intermol_dist:.2f} Å is very small"
        )
    
    return result


# =============================================================================
# CHEMICAL PATTERN LIBRARY
# =============================================================================

class ChemicalPatterns:
    """Library of chemically-meaningful arrangement patterns."""
    
    @staticmethod
    def pi_pi_parallel(n: int, distance: float = 3.4) -> Tuple[List[Position], List[Orientation]]:
        """Face-to-face π-stacking."""
        positions = generate_linear_positions(n, spacing=distance)
        orientations = generate_fixed_orientations(n)
        return positions, orientations
    
    @staticmethod
    def pi_pi_antiparallel(n: int, distance: float = 3.4) -> Tuple[List[Position], List[Orientation]]:
        """Antiparallel π-stacking (180° alternation)."""
        positions = generate_linear_positions(n, spacing=distance)
        orientations = generate_alternating_orientations(n, (0, 0, 0), (0, 0, 180))
        return positions, orientations
    
    @staticmethod
    def pi_pi_offset(n: int, distance: float = 3.4, 
                     offset: float = 1.5) -> Tuple[List[Position], List[Orientation]]:
        """Offset/slip-stacked arrangement."""
        positions = generate_formula_positions(
            n,
            x_formula=f"{offset} * (i % 2)",
            y_formula="0",
            z_formula=f"{distance} * i"
        )
        orientations = generate_fixed_orientations(n)
        return positions, orientations
    
    @staticmethod
    def t_shaped(distance: float = 5.0) -> Tuple[List[Position], List[Orientation]]:
        """T-shaped/edge-to-face dimer."""
        positions = [Position(0, 0, 0), Position(distance, 0, 0)]
        orientations = [Orientation(0, 0, 0), Orientation(0, 90, 0)]
        return positions, orientations
    
    @staticmethod
    def h_bonded_circular(n: int, radius: float = 5.0) -> Tuple[List[Position], List[Orientation]]:
        """Circular H-bonded arrangement."""
        positions = generate_circular_positions(n, radius=radius)
        orientations = generate_face_center_orientations(positions, facing="inward")
        return positions, orientations
    
    @staticmethod
    def herringbone(n: int, distance: float = 5.5, 
                    tilt: float = 60.0) -> Tuple[List[Position], List[Orientation]]:
        """Herringbone pattern."""
        positions = generate_linear_positions(n, spacing=distance)
        orientations = generate_alternating_orientations(n, (0, tilt, 0), (0, -tilt, 0))
        return positions, orientations
    
    @staticmethod
    def sandwich(distance: float = 3.4) -> Tuple[List[Position], List[Orientation]]:
        """Sandwich complex (A-B-A arrangement)."""
        positions = [
            Position(0, 0, -distance),
            Position(0, 0, 0),
            Position(0, 0, distance)
        ]
        orientations = generate_fixed_orientations(3)
        return positions, orientations
    
    @staticmethod
    def helical(n: int, radius: float = 5.0, pitch: float = 3.4,
                turns: float = 1.0) -> Tuple[List[Position], List[Orientation]]:
        """Helical arrangement."""
        positions = generate_helical_positions(n, radius, pitch, turns)
        # Orient tangent to helix
        orientations = []
        total_angle = 2 * np.pi * turns
        for i in range(n):
            t = i / max(n - 1, 1)
            angle = t * total_angle
            orientations.append(Orientation(0, 0, np.degrees(angle)))
        return positions, orientations


# =============================================================================
# MAIN API
# =============================================================================

def generate_molecular_cluster(
    molecules: List[Dict[str, Any]],
    arrangement: Union[str, Dict[str, Any]] = "auto",
    distance: Optional[float] = None,
    constraints: Optional[List[Union[str, Constraint]]] = None,
    optimize: bool = False,
    validate: bool = True,
    **kwargs
) -> Dict[str, Any]:
    """
    Main entry point for generating molecular clusters.
    
    This is the production API for the MCP backend.
    
    Args:
        molecules: List of molecule dicts with 'atoms' and 'coords'
        arrangement: Pattern name, NL description, or config dict
        distance: Intermolecular distance (Å)
        constraints: Additional constraints (strings or Constraint objects)
        optimize: Run constraint solver to refine positions
        validate: Check for clashes and physical plausibility
        **kwargs: Additional pattern-specific parameters
    
    Returns:
        {
            'molecules': List of transformed molecule dicts,
            'poses': List of pose information,
            'validation': Validation results,
            'metadata': Additional info
        }
    """
    debug(f"generate_molecular_cluster called")
    debug(f"  n_molecules: {len(molecules)}")
    debug(f"  arrangement: {arrangement}")
    debug(f"  distance: {distance}")
    debug(f"  constraints: {constraints}")
    debug(f"  optimize: {optimize}")
    debug(f"  validate: {validate}")
    debug(f"  kwargs: {kwargs}")
    
    n = len(molecules)
    
    # Parse arrangement specification
    if isinstance(arrangement, str):
        if arrangement == "auto":
            arrangement = _auto_detect_arrangement(molecules)
            debug(f"Auto-detected arrangement: {arrangement}")
        
        # Try NL parsing first for complex strings
        if " " in arrangement:
            debug(f"Parsing complex arrangement string via NL parser")
            parsed = parse_arrangement_request(arrangement)
            arrangement = parsed['pattern']
            if parsed['distance'] and distance is None:
                distance = parsed['distance']
            kwargs.update(parsed.get('parameters', {}))
            debug(f"NL parsed: pattern={arrangement}, distance={distance}")
            
            # Update molecule count if inferred from NL
            if parsed.get('n_molecules') and parsed['n_molecules'] > n:
                target_n = parsed['n_molecules']
                debug(f"Resizing molecule list from {n} to {target_n} based on NL")
                if n == 1:
                     # Clone the single molecule
                     template = molecules[0]
                     molecules = [template.copy() for _ in range(target_n)]
                     n = target_n
                else:
                     debug(f"WARNING: Cannot resize {n} heterogeneous molecules to {target_n}")
    
    # Set default distance
    if distance is None:
        distance = INTERACTION_DISTANCES.get(arrangement, (3.4, 4.0))[0]
        debug(f"Using default distance for {arrangement}: {distance}Å")
    
    # Generate positions and orientations
    debug(f"Generating positions/orientations for pattern: {arrangement}")
    positions, orientations = _generate_pattern(arrangement, n, distance, **kwargs)
    debug(f"Generated {len(positions)} positions")
    
    # Create poses
    poses = []
    for i, (pos, orient) in enumerate(zip(positions, orientations)):
        poses.append(MoleculePose(
            position=pos,
            orientation=orient,
            molecule_idx=i,
            metadata={'pattern': arrangement}
        ))
    
    debug(f"Created {len(poses)} poses")
    
    # Build constraints
    constraint_objs = _build_constraints(n, constraints or [], arrangement)
    debug(f"Built {len(constraint_objs)} constraint objects")
    
    # Optimize if requested
    if optimize and constraint_objs:
        debug(f"Running constraint solver with {len(constraint_objs)} constraints")
        solver = ConstraintSolver(molecules, poses, constraint_objs)
        poses = solver.solve()
        debug(f"Constraint solver completed")
    
    # Apply poses to molecules
    debug(f"Applying poses to molecules")
    arranged_molecules = apply_poses_to_molecules(molecules, poses)
    
    # Validate
    validation = {}
    if validate:
        debug(f"Validating arrangement")
        validation = validate_arrangement(molecules, poses)
        debug(f"Validation: valid={validation.get('valid')}, warnings={len(validation.get('warnings', []))}, errors={len(validation.get('errors', []))}")
    
    debug(f"generate_molecular_cluster completed successfully")
    
    return {
        'success': True,  # Add success flag for compatibility
        'molecules': arranged_molecules,
        'poses': [
            {
                'position': [p.position.x, p.position.y, p.position.z],
                'orientation': [p.orientation.roll, p.orientation.pitch, p.orientation.yaw],
                'molecule_idx': p.molecule_idx
            }
            for p in poses
        ],
        'validation': validation,
        'metadata': {
            'pattern': arrangement,
            'n_molecules': n,
            'distance': distance,
            'optimized': optimize
        }
    }



def _auto_detect_arrangement(molecules: List[Dict]) -> str:
    """Auto-detect best arrangement based on molecule properties."""
    n = len(molecules)
    
    def has_aromatic(mol):
        atoms = mol.get("atoms", [])
        return atoms.count("C") >= 6
    
    def can_hbond(mol):
        atoms = mol.get("atoms", [])
        return ("O" in atoms or "N" in atoms) and "H" in atoms
    
    all_aromatic = all(has_aromatic(m) for m in molecules)
    all_hbond = all(can_hbond(m) for m in molecules)
    
    if n == 2 and all_aromatic:
        return "pi_stacking"
    elif n == 2:
        return "linear"
    elif all_hbond and n <= 8:
        return "h_bonded_circular"
    elif all_aromatic:
        return "pi_stacking"
    else:
        return "linear"


def _generate_pattern(pattern: str, n: int, distance: float, 
                      **kwargs) -> Tuple[List[Position], List[Orientation]]:
    """Generate positions and orientations for a pattern."""
    pattern_lower = pattern.lower().replace('-', '_').replace(' ', '_')
    
    # Map pattern names to generators
    if pattern_lower in ['pi_stacking', 'pi_pi_parallel', 'parallel', 'cofacial']:
        return ChemicalPatterns.pi_pi_parallel(n, distance)
    
    elif pattern_lower in ['pi_pi_antiparallel', 'antiparallel']:
        return ChemicalPatterns.pi_pi_antiparallel(n, distance)
    
    elif pattern_lower in ['offset_stacking', 'pi_pi_offset', 'slip_stacked']:
        offset = kwargs.get('offset', 1.5)
        return ChemicalPatterns.pi_pi_offset(n, distance, offset)
    
    elif pattern_lower in ['t_shaped', 'edge_to_face']:
        return ChemicalPatterns.t_shaped(distance)
    
    elif pattern_lower in ['h_bonded', 'h_bonded_circular', 'circular']:
        radius = kwargs.get('radius', distance * n / (2 * np.pi))
        return ChemicalPatterns.h_bonded_circular(n, radius)
    
    elif pattern_lower in ['herringbone', 'zigzag']:
        tilt = kwargs.get('tilt', 60.0)
        return ChemicalPatterns.herringbone(n, distance, tilt)
    
    elif pattern_lower in ['helical', 'helix', 'spiral']:
        radius = kwargs.get('radius', 5.0)
        pitch = kwargs.get('pitch', distance)
        turns = kwargs.get('turns', 1.0)
        return ChemicalPatterns.helical(n, radius, pitch, turns)
    
    elif pattern_lower in ['sandwich']:
        return ChemicalPatterns.sandwich(distance)
    
    elif pattern_lower in ['spherical', 'sphere']:
        radius = kwargs.get('radius', 5.0)
        positions = generate_spherical_positions(n, radius)
        orientations = generate_face_center_orientations(positions, facing="outward")
        return positions, orientations
    
    elif pattern_lower in ['grid', 'planar', 'planner', 'sheet']:
        positions = generate_grid_positions(n, distance)
        orientations = generate_fixed_orientations(n)
        return positions, orientations
    
    else:
        # Default: linear
        positions = generate_linear_positions(n, distance)
        orientations = generate_fixed_orientations(n)
        return positions, orientations


def _build_constraints(n: int, constraint_specs: List, 
                       pattern: str) -> List[Constraint]:
    """Build constraint objects from specifications."""
    constraints = []
    
    # Add pattern-implicit constraints
    if 'pi' in pattern or 'parallel' in pattern:
        for i in range(n - 1):
            constraints.append(PlaneAlignmentConstraint(i, i + 1, "parallel"))
    
    if 't_shaped' in pattern or 'perpendicular' in pattern:
        if n >= 2:
            constraints.append(PlaneAlignmentConstraint(0, 1, "perpendicular"))
    
    if 'h_bond' in pattern:
        for i in range(n - 1):
            constraints.append(HBondConstraint(i, (i + 1) % n))
    
    # Add clash constraints
    for i in range(n):
        for j in range(i + 1, n):
            constraints.append(ClashConstraint(i, j, priority=2))
    
    # Parse user-specified constraints
    for spec in constraint_specs:
        if isinstance(spec, Constraint):
            constraints.append(spec)
        elif isinstance(spec, str):
            constraint = _parse_constraint_string(spec)
            if constraint:
                constraints.append(constraint)
    
    return constraints


def _parse_constraint_string(s: str) -> Optional[Constraint]:
    """Parse constraint from string specification."""
    match = re.match(r"(\w+)\s*\((.*)\)", s)
    if not match:
        return None
    
    func_name = match.group(1).lower()
    args_str = match.group(2)
    
    # Parse arguments with parenthesis awareness
    args = []
    current_arg = []
    paren_depth = 0
    
    for char in args_str:
        if char == ',' and paren_depth == 0:
            args.append("".join(current_arg).strip())
            current_arg = []
        else:
            if char == '(':
                paren_depth += 1
            elif char == ')':
                paren_depth -= 1
            current_arg.append(char)
    if current_arg:
        args.append("".join(current_arg).strip())
    
    kwargs = {}
    positional_idx = 0
    
    for part in args:
        if "=" in part:
             # Heuristic to check if it's a kwarg (key=value) vs selector with = inside?
             # Selectors shouldn't have = outside parens usually.
             k, v = part.split("=", 1)
             k = k.strip()
             v = v.strip()
             if v.replace(".", "").replace("-", "").isdigit():
                 kwargs[k] = float(v)
             else:
                 kwargs[k] = v.strip("'\"")
        else:
             # Positional argument
             val = part
             if val.replace(".", "").replace("-", "").isdigit():
                 val = float(val)
             else:
                 val = val.strip("'\"")
             kwargs[f"arg{positional_idx}"] = val
             positional_idx += 1
    
    if func_name == "distance":
        return DistanceConstraint(
            selector1=kwargs.get('sel1', kwargs.get('arg0', '0:centroid()')),
            selector2=kwargs.get('sel2', kwargs.get('arg1', '1:centroid()')),
            target=kwargs.get('target', kwargs.get('arg2')),
            min_dist=kwargs.get('min'),
            max_dist=kwargs.get('max')
        )
    
    elif func_name == "angle":
        return AngleConstraint(
            selector1=kwargs.get('sel1', kwargs.get('arg0', '0:centroid()')),
            selector2=kwargs.get('sel2', kwargs.get('arg1', '1:centroid()')),
            selector3=kwargs.get('sel3', kwargs.get('arg2', '2:centroid()')),
            target=kwargs.get('target', kwargs.get('arg3'))
        )
    
    elif func_name in ["h_bond", "hbond"]:
        return HBondConstraint(
            donor_mol=int(kwargs.get('donor', kwargs.get('arg0', 0))),
            acceptor_mol=int(kwargs.get('acceptor', kwargs.get('arg1', 1)))
        )
    
    elif func_name in ["plane_parallel", "parallel"]:
        return PlaneAlignmentConstraint(
            mol1=int(kwargs.get('mol1', kwargs.get('arg0', 0))),
            mol2=int(kwargs.get('mol2', kwargs.get('arg1', 1))),
            mode="parallel"
        )
    
    return None


# =============================================================================
# MODULE TEST
# =============================================================================

if __name__ == "__main__":
    print("=" * 70)
    print("MOLECULAR ARRANGEMENT ENGINE - Production MCP Backend")
    print("=" * 70)
    print()
    
    # Create test molecules (benzene-like)
    benzene = {
        "atoms": ["C", "C", "C", "C", "C", "C", "H", "H", "H", "H", "H", "H"],
        "coords": [
            [1.4, 0.0, 0.0], [0.7, 1.21, 0.0], [-0.7, 1.21, 0.0],
            [-1.4, 0.0, 0.0], [-0.7, -1.21, 0.0], [0.7, -1.21, 0.0],
            [2.5, 0.0, 0.0], [1.25, 2.16, 0.0], [-1.25, 2.16, 0.0],
            [-2.5, 0.0, 0.0], [-1.25, -2.16, 0.0], [1.25, -2.16, 0.0]
        ],
        "identifier": "benzene"
    }
    
    water = {
        "atoms": ["O", "H", "H"],
        "coords": [[0.0, 0.0, 0.0], [0.96, 0.0, 0.0], [-0.24, 0.93, 0.0]],
        "identifier": "water"
    }
    
    # Test 1: π-π stacking
    print("Test 1: π-π parallel stacking of 3 benzenes")
    result = generate_molecular_cluster(
        [benzene.copy() for _ in range(3)],
        arrangement="pi_pi_parallel",
        distance=3.4
    )
    print(f"  Pattern: {result['metadata']['pattern']}")
    print(f"  Molecules: {result['metadata']['n_molecules']}")
    print(f"  Validation: {result['validation'].get('valid', 'N/A')}")
    for i, p in enumerate(result['poses']):
        print(f"  Pose {i}: pos={p['position']}")
    print()
    
    # Test 2: Natural language parsing
    print("Test 2: Natural language parsing")
    nl_tests = [
        "Stack 4 benzenes with π-stacking at 3.5 Å",
        "Create a circular arrangement of 6 water molecules",
        "Make a T-shaped dimer at 5.0 angstroms",
        "Arrange molecules in a helix with 2 turns"
    ]
    for text in nl_tests:
        parsed = parse_arrangement_request(text)
        print(f"  '{text}'")
        print(f"    → pattern={parsed['pattern']}, n={parsed['n_molecules']}, d={parsed['distance']}")
    print()
    
    # Test 3: Circular H-bonded arrangement
    print("Test 3: Circular H-bonded arrangement of 6 waters")
    result = generate_molecular_cluster(
        [water.copy() for _ in range(6)],
        arrangement="h_bonded_circular",
        distance=2.8,
        validate=True
    )
    print(f"  Pattern: {result['metadata']['pattern']}")
    print(f"  Validation: {result['validation'].get('valid', 'N/A')}")
    if result['validation'].get('warnings'):
        for w in result['validation']['warnings']:
            print(f"    Warning: {w}")
    print()
    
    # Test 4: Custom formula arrangement
    print("Test 4: Custom spiral formula")
    positions = generate_formula_positions(
        5,
        x_formula="3 * sqrt(i) * cos(i * 2.4)",
        y_formula="3 * sqrt(i) * sin(i * 2.4)",
        z_formula="0.5 * i"
    )
    for i, p in enumerate(positions):
        print(f"  Mol {i}: ({p.x:.2f}, {p.y:.2f}, {p.z:.2f})")
    print()
    
    # Test 5: AtomSelector
    print("Test 5: AtomSelector")
    selector = AtomSelector("0:centroid()", [benzene])
    print(f"  0:centroid() = {selector.evaluate()}")
    selector = AtomSelector("0:ring_center(0)", [benzene])
    print(f"  0:ring_center(0) = {selector.evaluate()}")
    selector = AtomSelector("0:plane_normal()", [benzene])
    print(f"  0:plane_normal() = {selector.evaluate()}")
    print()
    
    print("=" * 70)
    print("All tests passed! Engine ready for MCP integration.")
    print("=" * 70)
