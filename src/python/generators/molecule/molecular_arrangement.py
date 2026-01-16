"""
molecular_arrangement.py - Complete Molecular Arrangement System

Unified implementation combining:
- Chemical Intelligence (AtomSelector, Constraints, Solver) from K2_plan
- Rigid Body Engine (Pose, Frame, Relative Placement) from z2_plan  
- Mathematical Generators (Position, Orientation, Formulas) from Original Plan

Entry point: arrange_molecules()
"""

from dataclasses import dataclass, field
from typing import Dict, Any, List, Optional, Tuple, Union, Callable
from abc import ABC, abstractmethod
import numpy as np
import re
import math
import sys
import os

# =============================================================================
# DEBUG LOGGING
# =============================================================================

DEBUG_LOG_FILE = os.environ.get('DEBUG_LOG_FILE', None)

def debug(msg: str) -> None:
    """Write debug message to stderr and optionally to file."""
    full_msg = f"[molecular_arrangement] {msg}"
    sys.stderr.write(full_msg + "\n")
    sys.stderr.flush()
    if DEBUG_LOG_FILE:
        with open(DEBUG_LOG_FILE, 'a') as f:
            f.write(full_msg + "\n")

# =============================================================================
# SECTION 1: Constants & Chemical Data
# =============================================================================

VDW_RADII = {
    'H': 1.20, 'He': 1.40, 'Li': 1.82, 'Be': 1.53, 'B': 1.92, 'C': 1.70,
    'N': 1.55, 'O': 1.52, 'F': 1.47, 'Ne': 1.54, 'Na': 2.27, 'Mg': 1.73,
    'Al': 1.84, 'Si': 2.10, 'P': 1.80, 'S': 1.80, 'Cl': 1.75, 'Ar': 1.88,
    'K': 2.75, 'Ca': 2.31, 'Fe': 2.04, 'Cu': 1.40, 'Zn': 1.39, 'Br': 1.85,
    'I': 1.98,
}

INTERACTION_DISTANCES = {
    'pi_pi_parallel': (3.3, 3.6),
    'pi_pi_offset': (3.3, 3.8),
    't_shaped': (4.5, 6.0),
    'h_bonded': (2.6, 3.2),
    'halogen_bond': (2.8, 3.5),
    'vdw_contact': (3.5, 4.5),
    'linear': (3.5, 5.0),
    'helical': (3.0, 4.0),
    'circular': (3.5, 5.0),
}

PATTERN_ALIASES = {
    'pi-pi': 'pi_pi_parallel', 'pi_stacking': 'pi_pi_parallel',
    'parallel': 'pi_pi_parallel', 'cofacial': 'pi_pi_parallel',
    'stacked': 'pi_pi_parallel', 'face_to_face': 'pi_pi_parallel',
    'antiparallel': 'pi_pi_antiparallel', 'anti_parallel': 'pi_pi_antiparallel',
    'offset': 'pi_pi_offset', 'slip_stacked': 'pi_pi_offset',
    't-shaped': 't_shaped', 'edge_to_face': 't_shaped', 'perpendicular': 't_shaped',
    'h-bonded': 'h_bonded', 'hydrogen_bonded': 'h_bonded', 'hbond': 'h_bonded',
    'helix': 'helical', 'spiral': 'helical',
    'ring': 'circular', 'cycle': 'circular',
    'zigzag': 'herringbone',
    'sheet': 'grid', 'planar': 'grid', 'planner': 'grid',
    'sphere': 'spherical',
}


# =============================================================================
# SECTION 2: Core Data Structures
# =============================================================================

@dataclass
class Position:
    """3D position vector with arithmetic operations."""
    x: float = 0.0
    y: float = 0.0
    z: float = 0.0
    
    def __add__(self, other: 'Position') -> 'Position':
        return Position(self.x + other.x, self.y + other.y, self.z + other.z)
    
    def __sub__(self, other: 'Position') -> 'Position':
        return Position(self.x - other.x, self.y - other.y, self.z - other.z)
    
    def __mul__(self, scalar: float) -> 'Position':
        return Position(self.x * scalar, self.y * scalar, self.z * scalar)
    
    def to_array(self) -> np.ndarray:
        return np.array([self.x, self.y, self.z])
    
    @classmethod
    def from_array(cls, arr) -> 'Position':
        return cls(float(arr[0]), float(arr[1]), float(arr[2]))
    
    def norm(self) -> float:
        return np.sqrt(self.x**2 + self.y**2 + self.z**2)
    
    def normalized(self) -> 'Position':
        n = self.norm()
        if n < 1e-10:
            return Position(0, 0, 1)
        return self * (1.0 / n)


@dataclass
class Orientation:
    """Euler angles (roll, pitch, yaw) in degrees."""
    roll: float = 0.0   # Rotation around x-axis
    pitch: float = 0.0  # Rotation around y-axis
    yaw: float = 0.0    # Rotation around z-axis
    
    def to_matrix(self) -> np.ndarray:
        """Convert to 3x3 rotation matrix (ZYX convention)."""
        r, p, y = np.radians([self.roll, self.pitch, self.yaw])
        cr, sr = np.cos(r), np.sin(r)
        cp, sp = np.cos(p), np.sin(p)
        cy, sy = np.cos(y), np.sin(y)
        
        return np.array([
            [cy*cp, cy*sp*sr - sy*cr, cy*sp*cr + sy*sr],
            [sy*cp, sy*sp*sr + cy*cr, sy*sp*cr - cy*sr],
            [-sp,   cp*sr,            cp*cr]
        ])
    
    @classmethod
    def from_matrix(cls, R: np.ndarray) -> 'Orientation':
        """Extract Euler angles from rotation matrix."""
        pitch = np.arcsin(-R[2, 0])
        if np.abs(np.cos(pitch)) > 1e-6:
            roll = np.arctan2(R[2, 1], R[2, 2])
            yaw = np.arctan2(R[1, 0], R[0, 0])
        else:
            roll = 0
            yaw = np.arctan2(-R[0, 1], R[1, 1])
        return cls(np.degrees(roll), np.degrees(pitch), np.degrees(yaw))


@dataclass
class MolecularFrame:
    """Local coordinate system for a molecule."""
    origin: np.ndarray  # Center of mass or geometric center
    x_axis: np.ndarray  # Primary axis (longest direction)
    y_axis: np.ndarray  # Secondary axis
    z_axis: np.ndarray  # Normal to molecular plane (for planar molecules)
    
    @classmethod
    def from_coords(cls, coords: np.ndarray, method: str = "svd") -> 'MolecularFrame':
        """Build frame from atomic coordinates using SVD."""
        coords = np.array(coords)
        origin = coords.mean(axis=0)
        centered = coords - origin
        
        if len(coords) < 3:
            return cls(origin, np.array([1,0,0]), np.array([0,1,0]), np.array([0,0,1]))
        
        U, S, Vt = np.linalg.svd(centered)
        x_axis = Vt[0] / np.linalg.norm(Vt[0])
        y_axis = Vt[1] / np.linalg.norm(Vt[1]) if len(Vt) > 1 else np.array([0,1,0])
        z_axis = np.cross(x_axis, y_axis)
        z_axis = z_axis / np.linalg.norm(z_axis) if np.linalg.norm(z_axis) > 1e-6 else np.array([0,0,1])
        
        return cls(origin, x_axis, y_axis, z_axis)
    
    def to_rotation_matrix(self) -> np.ndarray:
        """Get rotation matrix that transforms world coords to local frame."""
        return np.column_stack([self.x_axis, self.y_axis, self.z_axis])


@dataclass
class MoleculePose:
    """6-DOF pose specification for a molecule."""
    position: Position = field(default_factory=Position)
    orientation: Orientation = field(default_factory=Orientation)
    molecule_idx: int = 0
    parent_idx: Optional[int] = None  # For relative placement
    local_offset: Optional[Position] = None
    local_orientation: Optional[Orientation] = None
    metadata: Dict[str, Any] = field(default_factory=dict)
    
    def is_relative(self) -> bool:
        return self.parent_idx is not None


def resolve_relative_poses(
    poses: List[MoleculePose], 
    frames: List[MolecularFrame]
) -> List[MoleculePose]:
    """Resolve relative poses to absolute poses using topological sort."""
    n = len(poses)
    resolved = [None] * n
    
    # Find absolute poses first
    for i, pose in enumerate(poses):
        if not pose.is_relative():
            resolved[i] = pose
    
    # Iteratively resolve relative poses
    max_iterations = n
    for _ in range(max_iterations):
        all_resolved = True
        for i, pose in enumerate(poses):
            if resolved[i] is not None:
                continue
            if pose.parent_idx is not None and resolved[pose.parent_idx] is not None:
                parent = resolved[pose.parent_idx]
                parent_frame = frames[pose.parent_idx]
                
                # Transform local offset to world coords
                R = parent.orientation.to_matrix() @ parent_frame.to_rotation_matrix()
                offset = pose.local_offset or Position()
                world_offset = R @ offset.to_array()
                
                new_pos = Position.from_array(parent.position.to_array() + world_offset)
                
                # Combine orientations
                parent_R = parent.orientation.to_matrix()
                local_R = (pose.local_orientation or Orientation()).to_matrix()
                combined_R = parent_R @ local_R
                new_orient = Orientation.from_matrix(combined_R)
                
                resolved[i] = MoleculePose(
                    position=new_pos,
                    orientation=new_orient,
                    molecule_idx=pose.molecule_idx,
                    metadata=pose.metadata
                )
            else:
                all_resolved = False
        if all_resolved:
            break
    
    # Fallback for unresolved
    for i in range(n):
        if resolved[i] is None:
            resolved[i] = poses[i]
    
    return resolved


# =============================================================================
# SECTION 3: AtomSelector DSL
# =============================================================================

class AtomSelector:
    """
    Chemistry-aware atom/point selector with DSL syntax.
    
    Syntax: "mol_idx:target(args)"
    Examples:
        "0:centroid()"       - Center of mass of molecule 0
        "1:ring_center(0)"   - Center of first aromatic ring
        "0:plane_normal()"   - Normal vector to molecular plane
        "0:atom(O)"          - First oxygen atom
        "0:atom(O, 1)"       - Second oxygen atom
        "0:donor_h(0)"       - First H-bond donor hydrogen
        "0:acceptor(0)"      - First H-bond acceptor (O, N, F)
        "0:carbonyl()"       - Carbonyl oxygen (C=O)
        "0:func_group(amide)"- Amide group center
    """
    
    def __init__(self, selector: str, molecules: List[Dict[str, Any]]):
        self.selector = selector
        self.molecules = molecules
        self.mol_idx, self.target, self.args = self._parse(selector)
    
    def _parse(self, selector: str) -> Tuple[int, str, List[str]]:
        match = re.match(r'(\d+):(\w+)\((.*)\)', selector.strip())
        if not match:
            return 0, 'centroid', []
        mol_idx = int(match.group(1))
        target = match.group(2)
        args_str = match.group(3).strip()
        args = [a.strip() for a in args_str.split(',')] if args_str else []
        return mol_idx, target, args
    
    def evaluate(self) -> np.ndarray:
        """Evaluate selector and return 3D point or vector."""
        if self.mol_idx >= len(self.molecules):
            return np.array([0.0, 0.0, 0.0])
        
        mol = self.molecules[self.mol_idx]
        atoms = mol.get('atoms', [])
        coords = np.array(mol.get('coords', []))
        
        if len(coords) == 0:
            return np.array([0.0, 0.0, 0.0])
        
        method = getattr(self, f'_eval_{self.target}', None)
        if method:
            return method(atoms, coords)
        return self._eval_centroid(atoms, coords)
    
    def _eval_centroid(self, atoms, coords) -> np.ndarray:
        return coords.mean(axis=0)
    
    def _eval_atom(self, atoms, coords) -> np.ndarray:
        if not self.args:
            return coords[0] if len(coords) > 0 else np.zeros(3)
        element = self.args[0]
        idx = int(self.args[1]) if len(self.args) > 1 else 0
        matches = [i for i, a in enumerate(atoms) if a == element]
        matches = [i for i, a in enumerate(atoms) if a == element]
        if idx < len(matches):
            return coords[matches[idx]]
        raise ValueError(f"Atom not found: element={element}, index={idx}")
    
    def _eval_ring_center(self, atoms, coords) -> np.ndarray:
        ring_idx = int(self.args[0]) if self.args else 0
        c_indices = [i for i, a in enumerate(atoms) if a == 'C']
        if len(c_indices) >= 6:
            ring_atoms = c_indices[:6]
            return coords[ring_atoms].mean(axis=0)
        raise ValueError(f"Ring center not found: insufficient carbon atoms ({len(c_indices)})")
    
    def _eval_plane_normal(self, atoms, coords) -> np.ndarray:
        if len(coords) < 3:
            return np.array([0.0, 0.0, 1.0])
        centered = coords - coords.mean(axis=0)
        _, _, Vt = np.linalg.svd(centered)
        return Vt[-1]
    
    def _eval_donor_h(self, atoms, coords) -> np.ndarray:
        idx = int(self.args[0]) if self.args else 0
        h_indices = []
        for i, a in enumerate(atoms):
            if a == 'H':
                for j, other in enumerate(atoms):
                    if j != i and other in ('O', 'N') and np.linalg.norm(coords[i] - coords[j]) < 1.2:
                        h_indices.append(i)
                        break
        if idx < len(h_indices):
            return coords[h_indices[idx]]
        raise ValueError(f"Donor hydrogen not found: index={idx}")
    
    def _eval_acceptor(self, atoms, coords) -> np.ndarray:
        idx = int(self.args[0]) if self.args else 0
        acceptors = [i for i, a in enumerate(atoms) if a in ('O', 'N', 'F')]
        if idx < len(acceptors):
            return coords[acceptors[idx]]
        raise ValueError(f"Acceptor not found: index={idx}")
    
    def _eval_carbonyl(self, atoms, coords) -> np.ndarray:
        for i, a in enumerate(atoms):
            if a == 'O':
                for j, b in enumerate(atoms):
                    if b == 'C' and 1.15 < np.linalg.norm(coords[i] - coords[j]) < 1.35:
                        return coords[i]
        return self._eval_atom(atoms, coords)
    
    def _eval_func_group(self, atoms, coords) -> np.ndarray:
        group = self.args[0] if self.args else 'carbonyl'
        if group == 'carbonyl':
            return self._eval_carbonyl(atoms, coords)
        elif group == 'amide':
            o_pos = self._eval_carbonyl(atoms, coords)
            n_indices = [i for i, a in enumerate(atoms) if a == 'N']
            if n_indices:
                return (o_pos + coords[n_indices[0]]) / 2
        return coords.mean(axis=0)


# =============================================================================
# SECTION 4: Constraint System
# =============================================================================

class Constraint(ABC):
    """Abstract base class for arrangement constraints."""
    priority: int = 1
    
    @abstractmethod
    def evaluate(self, molecules: List[Dict], poses: List[MoleculePose]) -> float:
        """Return violation magnitude (0 = satisfied)."""
        pass
    
    def is_satisfied(self, molecules: List[Dict], poses: List[MoleculePose], tol: float = 0.1) -> bool:
        return self.evaluate(molecules, poses) < tol


@dataclass
class DistanceConstraint(Constraint):
    """Constraint on distance between two selected points."""
    selector1: str
    selector2: str
    target: Optional[float] = None
    min_dist: Optional[float] = None
    max_dist: Optional[float] = None
    weight: float = 1.0
    
    def evaluate(self, molecules: List[Dict], poses: List[MoleculePose]) -> float:
        transformed = _apply_poses_to_coords(molecules, poses)
        p1 = AtomSelector(self.selector1, transformed).evaluate()
        p2 = AtomSelector(self.selector2, transformed).evaluate()
        dist = np.linalg.norm(p1 - p2)
        
        if self.target is not None:
            return self.weight * abs(dist - self.target)
        if self.min_dist is not None and dist < self.min_dist:
            return self.weight * (self.min_dist - dist)
        if self.max_dist is not None and dist > self.max_dist:
            return self.weight * (dist - self.max_dist)
        return 0.0


@dataclass
class AngleConstraint(Constraint):
    """Constraint on angle formed by three points."""
    selector1: str
    selector2: str
    selector3: str
    target: float = 180.0
    weight: float = 1.0
    
    def evaluate(self, molecules: List[Dict], poses: List[MoleculePose]) -> float:
        transformed = _apply_poses_to_coords(molecules, poses)
        p1 = AtomSelector(self.selector1, transformed).evaluate()
        p2 = AtomSelector(self.selector2, transformed).evaluate()
        p3 = AtomSelector(self.selector3, transformed).evaluate()
        
        v1 = p1 - p2
        v2 = p3 - p2
        n1, n2 = np.linalg.norm(v1), np.linalg.norm(v2)
        if n1 < 1e-6 or n2 < 1e-6:
            return 0.0
        cos_angle = np.clip(np.dot(v1, v2) / (n1 * n2), -1, 1)
        angle = np.degrees(np.arccos(cos_angle))
        # Proper angular difference (handle wrap-around)
        diff = abs(angle - self.target)
        diff = min(diff, 360 - diff)  # Shortest path on circle
        return self.weight * diff / 180.0


@dataclass
class PlaneAlignmentConstraint(Constraint):
    """Constraint for parallel or perpendicular molecular planes."""
    mol1: int
    mol2: int
    mode: str = "parallel"  # "parallel" or "perpendicular"
    weight: float = 1.0
    
    def evaluate(self, molecules: List[Dict], poses: List[MoleculePose]) -> float:
        if self.mol1 >= len(molecules) or self.mol2 >= len(molecules):
            return 0.0
        transformed = _apply_poses_to_coords(molecules, poses)
        n1 = AtomSelector(f"{self.mol1}:plane_normal()", transformed).evaluate()
        n2 = AtomSelector(f"{self.mol2}:plane_normal()", transformed).evaluate()
        
        dot = abs(np.dot(n1, n2))
        if self.mode == "parallel":
            return self.weight * (1 - dot)
        else:  # perpendicular
            return self.weight * dot


@dataclass
class HBondConstraint(Constraint):
    """Constraint for H-bond geometry."""
    donor_mol: int
    acceptor_mol: int
    distance_range: Tuple[float, float] = (2.6, 3.2)
    angle_range: Tuple[float, float] = (150.0, 180.0)
    weight: float = 1.0
    
    def evaluate(self, molecules: List[Dict], poses: List[MoleculePose]) -> float:
        transformed = _apply_poses_to_coords(molecules, poses)
        h = AtomSelector(f"{self.donor_mol}:donor_h(0)", transformed).evaluate()
        acc = AtomSelector(f"{self.acceptor_mol}:acceptor(0)", transformed).evaluate()
        
        dist = np.linalg.norm(h - acc)
        if dist < self.distance_range[0]:
            return self.weight * (self.distance_range[0] - dist)
        if dist > self.distance_range[1]:
            return self.weight * (dist - self.distance_range[1])
        return 0.0


@dataclass
class ClashConstraint(Constraint):
    """Constraint to prevent atomic overlaps."""
    mol1: int
    mol2: int
    scale_factor: float = 0.8
    priority: int = 2
    
    def evaluate(self, molecules: List[Dict], poses: List[MoleculePose]) -> float:
        if self.mol1 >= len(molecules) or self.mol2 >= len(molecules):
            return 0.0
        transformed = _apply_poses_to_coords(molecules, poses)
        
        atoms1 = transformed[self.mol1].get('atoms', [])
        coords1 = np.array(transformed[self.mol1].get('coords', []))
        atoms2 = transformed[self.mol2].get('atoms', [])
        coords2 = np.array(transformed[self.mol2].get('coords', []))
        
        if len(coords1) == 0 or len(coords2) == 0:
            return 0.0
        
        total_violation = 0.0
        for i, (a1, c1) in enumerate(zip(atoms1, coords1)):
            r1 = VDW_RADII.get(a1, 1.7)
            for j, (a2, c2) in enumerate(zip(atoms2, coords2)):
                r2 = VDW_RADII.get(a2, 1.7)
                dist = np.linalg.norm(c1 - c2)
                min_dist = (r1 + r2) * self.scale_factor
                if dist < min_dist:
                    total_violation += (min_dist - dist) ** 2
        return np.sqrt(total_violation)


def _apply_poses_to_coords(molecules: List[Dict], poses: List[MoleculePose]) -> List[Dict]:
    """Helper to transform molecule coordinates by poses."""
    result = []
    for i, mol in enumerate(molecules):
        coords = np.array(mol.get('coords', []))
        if len(coords) == 0 or i >= len(poses):
            result.append(mol)
            continue
        
        pose = poses[i]
        R = pose.orientation.to_matrix()
        center = coords.mean(axis=0)
        centered = coords - center
        rotated = (R @ centered.T).T
        translated = rotated + pose.position.to_array()
        
        new_mol = mol.copy()
        new_mol['coords'] = translated.tolist()
        result.append(new_mol)
    return result


class ConstraintSolver:
    """Gradient descent solver for constraint satisfaction."""
    
    def __init__(self, molecules: List[Dict], poses: List[MoleculePose], 
                 constraints: List[Constraint]):
        self.molecules = molecules
        self.poses = [MoleculePose(
            position=Position(p.position.x, p.position.y, p.position.z),
            orientation=Orientation(p.orientation.roll, p.orientation.pitch, p.orientation.yaw),
            molecule_idx=p.molecule_idx,
            metadata=p.metadata
        ) for p in poses]
        self.constraints = sorted(constraints, key=lambda c: -c.priority)
    
    def solve(self, max_iterations: int = 100, tolerance: float = 1e-3,
              step_size: float = 0.1, max_time_seconds: float = 10.0) -> List[MoleculePose]:
        """Solve constraints iteratively with time limit."""
        import time
        start_time = time.time()
        n = len(self.poses)
        
        debug(f"ConstraintSolver.solve: {len(self.constraints)} constraints, max_iter={max_iterations}, max_time={max_time_seconds}s")
        
        for iteration in range(max_iterations):
            # Time limit check
            elapsed = time.time() - start_time
            if elapsed > max_time_seconds:
                debug(f"  Solver timeout after {elapsed:.2f}s at iteration {iteration}")
                break
            
            total_violation = sum(c.evaluate(self.molecules, self.poses) for c in self.constraints)
            
            if iteration % 20 == 0:
                debug(f"  iteration {iteration}: violation={total_violation:.4f}")
            
            if total_violation < tolerance:
                debug(f"  Converged at iteration {iteration}")
                break
            
            gradients = self._compute_gradients()
            
            total_change = 0.0
            for i in range(n):
                grad_pos, grad_rot = gradients[i]
                
                dx = step_size * grad_pos[0]
                dy = step_size * grad_pos[1]
                dz = step_size * grad_pos[2]
                # Apply all rotation gradients (roll, pitch, yaw)
                droll = step_size * grad_rot[0] * 10  # Scale factor for degrees
                dpitch = step_size * grad_rot[1] * 10
                dyaw = step_size * grad_rot[2] * 10
                
                total_change += abs(dx) + abs(dy) + abs(dz) + abs(droll) + abs(dpitch) + abs(dyaw)
                
                self.poses[i].position.x -= dx
                self.poses[i].position.y -= dy
                self.poses[i].position.z -= dz
                self.poses[i].orientation.roll -= droll
                self.poses[i].orientation.pitch -= dpitch
                self.poses[i].orientation.yaw -= dyaw
            
            # Escape stagnant states (e.g. parallel when wanting perpendicular)
            if total_violation > 0.1 and total_change < 1e-6:
                debug("    Solver stuck; applying random perturbation")
                for i in range(len(self.poses)):
                    self.poses[i].orientation.roll += np.random.normal(0, 5)
                    self.poses[i].orientation.pitch += np.random.normal(0, 5)
                    self.poses[i].orientation.yaw += np.random.normal(0, 5)
                # Restart step size
                step_size = 0.1
            
            step_size *= 0.99
        
        debug(f"  Solver finished in {time.time() - start_time:.2f}s")
        return self.poses
    
    def _compute_gradients(self) -> List[Tuple[np.ndarray, np.ndarray]]:
        """Compute numerical gradients for both position and rotation."""
        n = len(self.poses)
        gradients = []
        
        # Adaptive delta based on molecular size
        # Larger molecules need larger position deltas, but rotation delta stays small
        max_extent = 3.4  # Default fallback
        for mol in self.molecules:
            coords = mol.get('coords', [[0, 0, 0]])
            if len(coords) > 1:
                extent = np.ptp(np.array(coords), axis=0).max()
                max_extent = max(max_extent, extent)
        
        delta_pos = max(0.01, min(0.1, max_extent / 100))  # 1% of molecule size, clamped
        delta_rot = 0.5  # degrees - small enough for smooth gradients
        
        for i in range(n):
            grad_pos = np.zeros(3)
            grad_rot = np.zeros(3)  # [roll, pitch, yaw]
            
            # Position gradients
            for dim in range(3):
                attr = ['x', 'y', 'z'][dim]
                orig = getattr(self.poses[i].position, attr)
                
                setattr(self.poses[i].position, attr, orig + delta_pos)
                v_plus = sum(c.evaluate(self.molecules, self.poses) for c in self.constraints)
                
                setattr(self.poses[i].position, attr, orig - delta_pos)
                v_minus = sum(c.evaluate(self.molecules, self.poses) for c in self.constraints)
                
                setattr(self.poses[i].position, attr, orig)  # Restore
                grad_pos[dim] = (v_plus - v_minus) / (2 * delta_pos)
            
            # Rotation gradients (roll, pitch, yaw)
            for rdim, attr in enumerate(['roll', 'pitch', 'yaw']):
                orig = getattr(self.poses[i].orientation, attr)
                
                setattr(self.poses[i].orientation, attr, orig + delta_rot)
                v_plus = sum(c.evaluate(self.molecules, self.poses) for c in self.constraints)
                
                setattr(self.poses[i].orientation, attr, orig - delta_rot)
                v_minus = sum(c.evaluate(self.molecules, self.poses) for c in self.constraints)
                
                setattr(self.poses[i].orientation, attr, orig)  # Restore
                grad_rot[rdim] = (v_plus - v_minus) / (2 * delta_rot)
            
            gradients.append((grad_pos, grad_rot))
        
        return gradients

# =============================================================================
# SECTION 5: Position & Orientation Generators
# =============================================================================

def generate_linear_positions(n: int, spacing: float = 3.4, 
                              axis: Tuple[float, float, float] = (0, 0, 1)) -> List[Position]:
    """Generate positions along a line."""
    debug(f"generate_linear_positions: n={n}, spacing={spacing}, axis={axis}")
    axis_arr = np.array(axis)
    axis_arr = axis_arr / np.linalg.norm(axis_arr)
    positions = []
    for i in range(n):
        pos = axis_arr * spacing * i
        positions.append(Position(pos[0], pos[1], pos[2]))
        debug(f"  position[{i}] = ({pos[0]:.3f}, {pos[1]:.3f}, {pos[2]:.3f})")
    return positions


# -----------------------------------------------------------------------------
# Gap-Aware Positioning Helpers (Smart Spacing)
# -----------------------------------------------------------------------------

def _get_molecule_extent_along_axis(mol: Dict[str, Any], axis: np.ndarray) -> Tuple[float, float]:
    """
    Project molecule onto axis and return (min, max) extent including VdW radii.
    
    This is the core helper for gap-based positioning. It calculates how far
    a molecule extends along a given direction vector.
    
    Args:
        mol: Molecule dict with 'atoms' and 'coords'
        axis: Normalized unit vector for projection direction
    
    Returns:
        (min_extent, max_extent) in Angstroms, including VdW radii at extremes
    """
    coords = np.array(mol.get('coords', [[0, 0, 0]]))
    atoms = mol.get('atoms', ['C'] * len(coords))
    
    if len(coords) == 0:
        return 0.0, 0.0
    
    # Project all atoms onto axis: scalar projection = coord · axis
    projections = coords @ axis
    
    # Find extreme atoms and add their VdW radii
    min_idx = int(np.argmin(projections))
    max_idx = int(np.argmax(projections))
    
    min_atom = atoms[min_idx] if min_idx < len(atoms) else 'C'
    max_atom = atoms[max_idx] if max_idx < len(atoms) else 'C'
    
    min_vdw = VDW_RADII.get(min_atom, 1.7)
    max_vdw = VDW_RADII.get(max_atom, 1.7)
    
    return float(projections.min() - min_vdw), float(projections.max() + max_vdw)


def generate_gap_based_positions(
    molecules: List[Dict[str, Any]],
    gaps: Union[float, List[float]],
    axis: Tuple[float, float, float] = (0, 0, 1),
    orientations: Optional[List['Orientation']] = None
) -> List[Position]:
    """
    Generate positions ensuring surface separation (gap) between molecules.
    
    This is the SMART position generator that replaces naive center-to-center
    spacing. It calculates molecular extents and places each molecule so that
    the gap between VdW surfaces equals the requested value.
    
    Args:
        molecules: List of molecule dicts with 'atoms' and 'coords'
        gaps: Single gap value (Å) OR list of gaps between consecutive pairs
        axis: Placement direction vector (will be normalized)
        orientations: Pre-computed orientations (for rotation-aware extent)
    
    Returns:
        List of Position objects with gap-aware spacing
    
    Example:
        # Place 3 PTCDAs with 3.4Å gap between surfaces
        positions = generate_gap_based_positions(
            molecules=[ptcda1, ptcda2, ptcda3],
            gaps=3.4,
            axis=(0, 0, 1)
        )
        # Result: ~18Å center-to-center (14Å molecule + 3.4Å gap)
    """
    debug(f"generate_gap_based_positions: n={len(molecules)}, gaps={gaps}, axis={axis}")
    
    axis_arr = np.array(axis, dtype=float)
    norm = np.linalg.norm(axis_arr)
    if norm < 1e-10:
        axis_arr = np.array([0, 0, 1], dtype=float)
    else:
        axis_arr = axis_arr / norm
    
    n = len(molecules)
    if n == 0:
        return []
    
    # Normalize gaps to a list
    if isinstance(gaps, (int, float)):
        gap_list = [float(gaps)] * max(n - 1, 1)
    else:
        gap_list = [float(g) for g in gaps]
        # Extend if needed
        while len(gap_list) < n - 1:
            gap_list.append(gap_list[-1] if gap_list else 3.4)
    
    positions = []
    current_pos = 0.0
    prev_max = 0.0
    
    for i in range(n):
        mol = molecules[i]
        
        # If orientations provided, apply rotation before calculating extent
        if orientations and i < len(orientations):
            R = orientations[i].to_matrix()
            coords = np.array(mol.get('coords', [[0, 0, 0]]))
            if len(coords) > 0:
                center = coords.mean(axis=0)
                rotated = (R @ (coords - center).T).T
                temp_mol = {**mol, 'coords': rotated.tolist()}
                min_ext, max_ext = _get_molecule_extent_along_axis(temp_mol, axis_arr)
            else:
                min_ext, max_ext = 0.0, 0.0
        else:
            min_ext, max_ext = _get_molecule_extent_along_axis(mol, axis_arr)
        
        if i == 0:
            # First molecule: place so its min_extent is at origin
            current_pos = -min_ext
        else:
            # Subsequent: place so gap exists between prev_max and curr_min
            gap = gap_list[i - 1]
            current_pos = prev_max + gap - min_ext
        
        prev_max = current_pos + max_ext
        
        # Convert scalar position along axis to 3D position
        pos_3d = axis_arr * current_pos
        positions.append(Position(float(pos_3d[0]), float(pos_3d[1]), float(pos_3d[2])))
        
        debug(f"  gap_position[{i}]: pos={current_pos:.3f}, extent=[{min_ext:.2f}, {max_ext:.2f}]")
    
    return positions



def direction_from_angles(azimuth: float = 0.0, elevation: float = 0.0) -> Tuple[float, float, float]:
    """
    Convert spherical angles to unit direction vector.
    
    This enables natural language like "along 30 degree direction" to be
    converted to a proper placement vector.
    
    Args:
        azimuth: Angle in XY plane from X-axis (degrees). 0°=X, 90°=Y
        elevation: Angle from XY plane (degrees). 0°=horizontal, 90°=Z
    
    Returns:
        (x, y, z) unit vector
    
    Examples:
        direction_from_angles(0, 0)   -> (1, 0, 0)  # Along X
        direction_from_angles(90, 0)  -> (0, 1, 0)  # Along Y
        direction_from_angles(0, 90)  -> (0, 0, 1)  # Along Z
        direction_from_angles(30, 0)  -> (0.866, 0.5, 0)  # 30° in XY plane
    """
    az_rad = np.radians(azimuth)
    el_rad = np.radians(elevation)
    x = np.cos(el_rad) * np.cos(az_rad)
    y = np.cos(el_rad) * np.sin(az_rad)
    z = np.sin(el_rad)
    return (float(x), float(y), float(z))


# =============================================================================
# SMART POSITION GENERATOR - GENERIC FOR ALL SCENARIOS
# =============================================================================

def generate_smart_positions(
    molecules: List[Dict[str, Any]],
    distance: Optional[float],
    axis: Tuple[float, float, float],
    orientations: Optional[List['Orientation']] = None,
    allow_clashes: bool = False
) -> Tuple[List[Position], Dict[str, Any]]:
    """
    SMART position placement for ANY scenario.
    
    This is the UNIFIED position generator that handles:
    - Any axis (X, Y, Z, or arbitrary vectors)
    - Asymmetric molecules (calculates extent along placement axis)
    - Clash prevention (uses safe distance when needed)
    - User notification (returns info about adjustments)
    
    Algorithm:
    1. If distance is None/0, treat as requesting minimum safe distance
    2. Calculate molecular extent along the specified axis
    3. Compute minimum_safe_distance = max half-extents + VdW buffer
    4. If requested distance >= safe_distance: use requested distance
       Else: use safe_distance (and mark as adjusted)
    5. Place molecules along axis at that distance
    
    Args:
        molecules: List of molecule dicts with 'atoms' and 'coords'
        distance: User's requested center-to-center distance (or None/0)
        axis: Placement direction vector (will be normalized)
        orientations: Optional pre-computed orientations for extent calculation
        allow_clashes: If True, use distance even if clashes occur
    
    Returns:
        (positions, position_info) where position_info contains:
            - requested_distance: what user asked for
            - actual_distance: what was used
            - was_adjusted: True if we had to change the distance
            - minimum_safe_distance: calculated safe distance
            - axis_used: the normalized axis tuple
            - adjustment_reason: explanation if adjusted
    """
    n = len(molecules)
    if n == 0:
        return [], {'error': 'No molecules provided'}
    
    # Normalize axis
    axis_arr = np.array(axis, dtype=float)
    norm = np.linalg.norm(axis_arr)
    if norm < 1e-10:
        axis_arr = np.array([0, 0, 1], dtype=float)  # Default to Z
    else:
        axis_arr = axis_arr / norm
    
    debug(f"generate_smart_positions: n={n}, distance={distance}, axis={tuple(axis_arr)}")
    
    # =========================================================================
    # PRE-ALIGNMENT FOR EXTENT CALCULATION
    # =========================================================================
    # For stacking of planar molecules, we need to align them so their plane
    # is PERPENDICULAR to the stacking axis. This ensures minimal extent
    # along the stacking direction, preventing clashes.
    # This now works for ANY axis, not just Z.
    # =========================================================================
    
    aligned_molecules = []
    for mol in molecules:
        # Align planar molecules perpendicular to stacking axis
        aligned_mol = _align_planar_molecule_to_axis(mol, axis_arr)
        aligned_molecules.append(aligned_mol)
    
    # Calculate molecular extents along the specified axis
    extents = []
    max_dimensions = []  # Store max dimension of each molecule for fallback
    
    for i, mol in enumerate(aligned_molecules):
        coords = np.array(mol.get('coords', [[0, 0, 0]]))
        
        # Calculate max molecular dimension (diameter/span)
        if len(coords) > 1:
            from scipy.spatial.distance import pdist
            max_dim = float(np.max(pdist(coords))) + 3.4  # Add VdW diameter
        else:
            max_dim = 3.4
        max_dimensions.append(max_dim)
        
        if orientations and i < len(orientations):
            # Apply orientation before calculating extent
            R = orientations[i].to_matrix()
            if len(coords) > 0:
                center = coords.mean(axis=0)
                rotated = (R @ (coords - center).T).T
                temp_mol = {**mol, 'coords': rotated.tolist()}
                min_ext, max_ext = _get_molecule_extent_along_axis(temp_mol, axis_arr)
            else:
                min_ext, max_ext = 0.0, 0.0
        else:
            min_ext, max_ext = _get_molecule_extent_along_axis(mol, axis_arr)
        
        extent_size = max_ext - min_ext
        extents.append({
            'min': min_ext,
            'max': max_ext,
            'size': extent_size,
            'half': extent_size / 2
        })
        debug(f"  molecule[{i}]: extent=[{min_ext:.2f}, {max_ext:.2f}], size={extent_size:.2f}Å, max_dim={max_dim:.2f}Å")
    
    # Calculate minimum safe distance for adjacent molecules
    # Use MAX of: (extent-based safe) and (max_dimension-based safe)
    # This handles asymmetric molecules that might clash perpendicular to stacking axis
    VDW_BUFFER = 0.5  # Small buffer to prevent surface contact
    
    min_safe_distances = []
    for i in range(n - 1):
        # Extent-based: half-extent of mol_i + half-extent of mol_{i+1}
        extent_safe = extents[i]['half'] + extents[i + 1]['half'] + VDW_BUFFER
        
        # Dimension-based: max dimension of larger molecule (fallback for asymmetric)
        dim_safe = max(max_dimensions[i], max_dimensions[i + 1]) / 2 + VDW_BUFFER
        
        # Use the larger of the two for safety
        safe_dist = max(extent_safe, dim_safe)
        min_safe_distances.append(safe_dist)
        debug(f"  safe_distance[{i}→{i+1}] = {safe_dist:.2f}Å (extent={extent_safe:.2f}, dim={dim_safe:.2f})")
    
    max_safe = max(min_safe_distances) if min_safe_distances else 3.4
    
    # Determine actual distance to use
    requested = distance if distance is not None else 0.0
    was_adjusted = False
    adjustment_reason = None
    actual_distance = requested
    
    if requested <= 0:
        # User didn't specify distance - use safe distance
        actual_distance = max_safe
        was_adjusted = True
        adjustment_reason = f"No distance specified, using minimum safe distance ({max_safe:.2f}Å)"
        debug(f"  No distance specified -> using safe distance {max_safe:.2f}Å")
    elif requested < max_safe and not allow_clashes:
        # User's distance is too small - would cause clashes
        actual_distance = max_safe
        was_adjusted = True
        adjustment_reason = f"Requested {requested:.2f}Å would cause clashes, adjusted to {max_safe:.2f}Å"
        debug(f"  Distance {requested}Å too small -> adjusted to {max_safe:.2f}Å")
    else:
        debug(f"  Using requested distance {requested}Å")
    
    # Generate positions along the axis - CUMULATIVE per-pair approach
    # This ensures each pair gets the appropriate distance based on their specific extents
    positions = [Position(0.0, 0.0, 0.0)]  # First molecule at origin
    cumulative_offset = 0.0
    
    for i in range(1, n):
        # Get per-pair safe distance (or max_safe as fallback)
        pair_safe = min_safe_distances[i - 1] if i - 1 < len(min_safe_distances) else max_safe
        
        # Use user's distance if larger than safe, otherwise use safe distance
        if allow_clashes:
            pair_distance = actual_distance
        else:
            pair_distance = max(actual_distance, pair_safe)
        
        cumulative_offset += pair_distance
        pos_3d = axis_arr * cumulative_offset
        positions.append(Position(float(pos_3d[0]), float(pos_3d[1]), float(pos_3d[2])))
        debug(f"  position[{i}] = ({pos_3d[0]:.3f}, {pos_3d[1]:.3f}, {pos_3d[2]:.3f}) [gap={pair_distance:.2f}Å]")
    
    position_info = {
        'requested_distance': requested,
        'actual_distance': actual_distance,
        'was_adjusted': was_adjusted,
        'minimum_safe_distance': max_safe,
        'axis_used': tuple(axis_arr.tolist()),
        'adjustment_reason': adjustment_reason,
        'molecular_extents': [e['size'] for e in extents]
    }
    
    return positions, position_info


# Patterns that should use SMART position calculation
GAP_BASED_PATTERNS = {
    'linear', 'pi_pi_parallel', 'pi_pi_antiparallel', 'pi_pi_offset', 
    'herringbone', 'sandwich', 'stacked'
}



def generate_circular_positions(n: int, radius: float = 5.0,
                                center: Tuple[float, float, float] = (0, 0, 0)) -> List[Position]:
    """Generate positions on a circle in the XY plane."""
    debug(f"generate_circular_positions: n={n}, radius={radius}")
    positions = []
    for i in range(n):
        theta = 2 * np.pi * i / n
        x = center[0] + radius * np.cos(theta)
        y = center[1] + radius * np.sin(theta)
        z = center[2]
        positions.append(Position(x, y, z))
        debug(f"  position[{i}] = ({x:.3f}, {y:.3f}, {z:.3f})")
    return positions


def generate_helical_positions(n: int, radius: float = 5.0, pitch: float = 3.4,
                               turns: float = 1.0) -> List[Position]:
    """Generate positions along a helix."""
    debug(f"generate_helical_positions: n={n}, radius={radius}, pitch={pitch}, turns={turns}")
    positions = []
    total_angle = 2 * np.pi * turns
    for i in range(n):
        t = i / max(n - 1, 1)
        theta = total_angle * t
        x = radius * np.cos(theta)
        y = radius * np.sin(theta)
        z = pitch * turns * t
        positions.append(Position(x, y, z))
        debug(f"  position[{i}] = ({x:.3f}, {y:.3f}, {z:.3f})")
    return positions


def generate_spherical_positions(n: int, radius: float = 5.0) -> List[Position]:
    """Generate positions on a sphere using Fibonacci lattice."""
    debug(f"generate_spherical_positions: n={n}, radius={radius}")
    positions = []
    golden_ratio = (1 + np.sqrt(5)) / 2
    for i in range(n):
        theta = 2 * np.pi * i / golden_ratio
        phi = np.arccos(1 - 2 * (i + 0.5) / n)
        x = radius * np.sin(phi) * np.cos(theta)
        y = radius * np.sin(phi) * np.sin(theta)
        z = radius * np.cos(phi)
        positions.append(Position(x, y, z))
        debug(f"  position[{i}] = ({x:.3f}, {y:.3f}, {z:.3f})")
    return positions


def generate_grid_positions(n: int, spacing: float = 3.4) -> List[Position]:
    """Generate positions in a 2D grid (XY plane)."""
    debug(f"generate_grid_positions: n={n}, spacing={spacing}")
    cols = int(np.ceil(np.sqrt(n)))
    positions = []
    for i in range(n):
        row = i // cols
        col = i % cols
        x = col * spacing - (cols - 1) * spacing / 2
        y = row * spacing - ((n - 1) // cols) * spacing / 2
        positions.append(Position(x, y, 0))
        debug(f"  position[{i}] = ({x:.3f}, {y:.3f}, 0.0)")
    return positions


class FormulaError(Exception):
    """Error in formula evaluation."""
    pass


def _safe_sqrt(x):
    """Safe sqrt that returns 0 for negative values instead of NaN."""
    return 0.0 if x < 0 else np.sqrt(x)


def _eval_formula(formula: str, local_vars: dict, axis: str, i: int) -> float:
    """
    Safely evaluate a formula string.

    Security: eval() is sandboxed with __builtins__={} and only whitelisted
    math functions are available. This is necessary for mathematical formula
    evaluation in this scientific computing context.
    """
    result = eval(formula, {"__builtins__": {}}, local_vars)  # noqa: S307
    return float(result)


def generate_formula_positions(n: int, x_formula: str = "0", y_formula: str = "0",
                               z_formula: str = "i * 3.4") -> List[Position]:
    """
    Generate positions using mathematical formulas.

    Handles errors gracefully:
    - Division by zero returns 0
    - sqrt of negative returns 0
    - Invalid functions raise FormulaError with helpful message
    - NaN/Inf values are replaced with 0
    - Coordinates > 1000 Å are clamped (physically unreasonable)
    """
    debug(f"generate_formula_positions: n={n}")
    debug(f"  x_formula = '{x_formula}'")
    debug(f"  y_formula = '{y_formula}'")
    debug(f"  z_formula = '{z_formula}'")

    if n <= 0:
        return []

    positions = []

    # Whitelisted safe math functions with common aliases
    safe_dict = {
        # Trig functions
        'sin': np.sin, 'cos': np.cos, 'tan': np.tan,
        'radians': np.radians, 'degrees': np.degrees,
        'rad': np.radians, 'deg': np.degrees,
        'asin': np.arcsin, 'acos': np.arccos, 'atan': np.arctan,
        'arcsin': np.arcsin, 'arccos': np.arccos, 'arctan': np.arctan,
        'sinh': np.sinh, 'cosh': np.cosh, 'tanh': np.tanh,
        # Power/log
        'sqrt': _safe_sqrt, 'exp': np.exp, 'log': np.log, 'log10': np.log10,
        'pow': pow,
        # Rounding
        'abs': np.abs, 'floor': np.floor, 'ceil': np.ceil, 'round': np.round,
        # Other
        'atan2': np.arctan2,
        'min': min, 'max': max,
        'int': int, 'float': float,
        # Constants
        'pi': np.pi, 'e': np.e,
        'phi': (1 + np.sqrt(5)) / 2,  # Golden ratio
    }

    # Grid dimensions for 3D decomposition
    side = int(np.ceil(n**(1/3))) if n > 0 else 1

    # Maximum reasonable coordinate (1000 Å = 100 nm)
    MAX_COORD = 1000.0

    for i in range(n):
        # Avoid division by zero in t calculation
        t = i / (n - 1) if n > 1 else 0.0

        # Grid variables for 3D formulas
        iz = i // (side * side)
        rem = i % (side * side)
        iy = rem // side
        ix = rem % side

        local_vars = {
            'i': i, 'n': n, 't': t,
            'j': iy, 'k': iz,
            'ix': ix, 'iy': iy, 'iz': iz,
            **safe_dict
        }

        # Evaluate each formula with error handling
        coords = []
        for axis, formula in [('x', x_formula), ('y', y_formula), ('z', z_formula)]:
            try:
                value = _eval_formula(formula, local_vars, axis, i)
            except ZeroDivisionError:
                debug(f"  WARNING: Division by zero in {axis}_formula at i={i}, using 0")
                value = 0.0
            except NameError as err:
                raise FormulaError(
                    f"Unknown function in {axis}_formula: {err}. "
                    f"Available: sin, cos, tan, radians, degrees, sqrt, exp, log, abs, floor, ceil, "
                    f"asin, acos, atan, atan2, sinh, cosh, tanh, pow, min, max, pi, e, phi"
                )
            except SyntaxError as err:
                raise FormulaError(f"Syntax error in {axis}_formula '{formula}': {err}")
            except Exception as err:
                raise FormulaError(f"Error evaluating {axis}_formula '{formula}' at i={i}: {err}")

            # Check for NaN or Inf
            if np.isnan(value) or np.isinf(value):
                debug(f"  WARNING: {axis}_formula produced {'NaN' if np.isnan(value) else 'Inf'} at i={i}, using 0")
                value = 0.0

            # Clamp unreasonable coordinates
            if abs(value) > MAX_COORD:
                debug(f"  WARNING: {axis}={value:.2e} Å exceeds {MAX_COORD}, clamping")
                value = MAX_COORD if value > 0 else -MAX_COORD

            coords.append(value)

        x, y, z = coords
        positions.append(Position(x, y, z))
        debug(f"  position[{i}] = ({x:.3f}, {y:.3f}, {z:.3f})")

    return positions


def generate_fixed_orientations(n: int, roll: float = 0, pitch: float = 0,
                                yaw: float = 0) -> List[Orientation]:
    """Generate identical orientations for all molecules."""
    debug(f"generate_fixed_orientations: n={n}, roll={roll}, pitch={pitch}, yaw={yaw}")
    return [Orientation(roll, pitch, yaw) for _ in range(n)]


def generate_alternating_orientations(n: int, 
                                      orient1: Tuple[float, float, float] = (0, 0, 0),
                                      orient2: Tuple[float, float, float] = (0, 0, 180)) -> List[Orientation]:
    """Generate alternating orientations (e.g., for antiparallel stacking)."""
    debug(f"generate_alternating_orientations: n={n}, orient1={orient1}, orient2={orient2}")
    orientations = []
    for i in range(n):
        if i % 2 == 0:
            orientations.append(Orientation(*orient1))
        else:
            orientations.append(Orientation(*orient2))
    return orientations


def generate_incremental_orientations(n: int, axis: str = 'z', 
                                      increment: float = 30.0) -> List[Orientation]:
    """Generate orientations with incremental rotation."""
    debug(f"generate_incremental_orientations: n={n}, axis={axis}, increment={increment}")
    orientations = []
    for i in range(n):
        angle = i * increment
        if axis == 'x':
            orientations.append(Orientation(roll=angle))
        elif axis == 'y':
            orientations.append(Orientation(pitch=angle))
        else:
            orientations.append(Orientation(yaw=angle))
    return orientations


def generate_face_center_orientations(positions: List[Position],
                                      facing: str = "inward") -> List[Orientation]:
    """Generate orientations facing toward/away from center."""
    debug(f"generate_face_center_orientations: n={len(positions)}, facing={facing}")
    center = Position(
        sum(p.x for p in positions) / len(positions),
        sum(p.y for p in positions) / len(positions),
        sum(p.z for p in positions) / len(positions)
    )
    orientations = []
    for pos in positions:
        direction = center - pos if facing == "inward" else pos - center
        yaw = np.degrees(np.arctan2(direction.y, direction.x))
        horizontal = np.sqrt(direction.x**2 + direction.y**2)
        pitch = np.degrees(np.arctan2(direction.z, horizontal))
        orientations.append(Orientation(0, pitch, yaw))
    return orientations


def generate_formula_orientations(n: int, positions: List[Position],
                                  roll_formula: str = "0",
                                  pitch_formula: str = "0",
                                  yaw_formula: str = "0") -> List[Orientation]:
    """Generate orientations using formulas."""
    debug(f"generate_formula_orientations: n={n}")
    safe_dict = {
        'sin': np.sin, 'cos': np.cos, 'tan': np.tan,
        'sqrt': np.sqrt, 'pi': np.pi, 'e': np.e,
    }
    orientations = []
    for i in range(n):
        pos = positions[i] if i < len(positions) else Position()
        t = i / max(n - 1, 1)
        local_vars = {'i': i, 'n': n, 't': t, 'x': pos.x, 'y': pos.y, 'z': pos.z, **safe_dict}
        
        roll = float(eval(roll_formula, {"__builtins__": {}}, local_vars))
        pitch = float(eval(pitch_formula, {"__builtins__": {}}, local_vars))
        yaw = float(eval(yaw_formula, {"__builtins__": {}}, local_vars))
        
        orientations.append(Orientation(roll, pitch, yaw))
    return orientations


# =============================================================================
# SECTION 6: Chemical Pattern Library
# =============================================================================

class ChemicalPatterns:
    """Library of chemically-meaningful arrangement patterns."""
    
    @staticmethod
    def pi_pi_parallel(n: int, distance: float = 3.4, axis: Tuple[float, float, float] = (0, 0, 1), **kwargs) -> Tuple[List[Position], List[Orientation]]:
        """Face-to-face π-stacking with all molecules parallel."""
        debug(f"ChemicalPatterns.pi_pi_parallel: n={n}, distance={distance}, axis={axis}")
        positions = generate_linear_positions(n, spacing=distance, axis=axis)
        orientations = generate_fixed_orientations(n)
        return positions, orientations
    
    @staticmethod
    def pi_pi_antiparallel(n: int, distance: float = 3.4, axis: Tuple[float, float, float] = (0, 0, 1), **kwargs) -> Tuple[List[Position], List[Orientation]]:
        """Face-to-face stacking with alternating 180° rotation."""
        debug(f"ChemicalPatterns.pi_pi_antiparallel: n={n}, distance={distance}, axis={axis}")
        positions = generate_linear_positions(n, spacing=distance, axis=axis)
        orientations = generate_alternating_orientations(n, (0, 0, 0), (0, 0, 180))
        return positions, orientations
    
    @staticmethod
    def pi_pi_offset(n: int, distance: float = 3.4,
                     offset: float = 1.5, **kwargs) -> Tuple[List[Position], List[Orientation]]:
        """Slip-stacked arrangement with lateral offset."""
        debug(f"ChemicalPatterns.pi_pi_offset: n={n}, distance={distance}, offset={offset}")
        positions = []
        for i in range(n):
            x = offset * (i % 2)
            z = distance * i
            positions.append(Position(x, 0, z))
        orientations = generate_fixed_orientations(n)
        return positions, orientations
    
    @staticmethod
    def t_shaped(distance: float = 5.0, **kwargs) -> Tuple[List[Position], List[Orientation]]:
        """Edge-to-face perpendicular dimer."""
        debug(f"ChemicalPatterns.t_shaped: distance={distance}")
        positions = [Position(0, 0, 0), Position(distance, 0, 0)]
        orientations = [Orientation(0, 0, 0), Orientation(0, 90, 0)]
        return positions, orientations
    
    @staticmethod
    def herringbone(n: int, distance: float = 5.5,
                    tilt: float = 45.0, **kwargs) -> Tuple[List[Position], List[Orientation]]:
        """Alternating tilt pattern (common in organic crystals)."""
        debug(f"ChemicalPatterns.herringbone: n={n}, distance={distance}, tilt={tilt}")
        positions = generate_linear_positions(n, spacing=distance, axis=(0, 0, 1))
        orientations = []
        for i in range(n):
            if i % 2 == 0:
                orientations.append(Orientation(tilt, 0, 0))
            else:
                orientations.append(Orientation(-tilt, 0, 0))
        return positions, orientations
    
    @staticmethod
    def h_bonded_dimer(distance: float = 2.8, **kwargs) -> Tuple[List[Position], List[Orientation]]:
        """Hydrogen-bonded dimer with optimal geometry."""
        debug(f"ChemicalPatterns.h_bonded_dimer: distance={distance}")
        positions = [Position(0, 0, 0), Position(distance, 0, 0)]
        orientations = [Orientation(0, 0, 0), Orientation(0, 0, 180)]
        return positions, orientations
    
    @staticmethod
    def h_bonded_circular(n: int, radius: float = 5.0, **kwargs) -> Tuple[List[Position], List[Orientation]]:
        """Circular ring with H-bond connectivity."""
        debug(f"ChemicalPatterns.h_bonded_circular: n={n}, radius={radius}")
        positions = generate_circular_positions(n, radius)
        orientations = generate_face_center_orientations(positions, facing="inward")
        return positions, orientations
    
    @staticmethod
    def helical(n: int, radius: float = 5.0, pitch: float = 3.4,
                turns: float = 1.0, **kwargs) -> Tuple[List[Position], List[Orientation]]:
        """Helical (DNA-like) arrangement."""
        debug(f"ChemicalPatterns.helical: n={n}, radius={radius}, pitch={pitch}, turns={turns}")
        positions = generate_helical_positions(n, radius, pitch, turns)
        orientations = generate_face_center_orientations(positions, facing="inward")
        return positions, orientations
    
    @staticmethod
    def sandwich(distance: float = 3.4, **kwargs) -> Tuple[List[Position], List[Orientation]]:
        """A-B-A intercalated sandwich arrangement."""
        debug(f"ChemicalPatterns.sandwich: distance={distance}")
        positions = [
            Position(0, 0, 0),
            Position(0, 0, distance),
            Position(0, 0, 2 * distance)
        ]
        orientations = generate_fixed_orientations(3)
        return positions, orientations
    
    @staticmethod
    def grid(n: int, spacing: float = 3.4, **kwargs) -> Tuple[List[Position], List[Orientation]]:
        """2D grid arrangement in XY plane."""
        debug(f"ChemicalPatterns.grid: n={n}, spacing={spacing}")
        positions = generate_grid_positions(n, spacing)
        orientations = generate_fixed_orientations(n)
        return positions, orientations
    
    @staticmethod
    def spherical(n: int, radius: float = 5.0, **kwargs) -> Tuple[List[Position], List[Orientation]]:
        """Molecules distributed on sphere surface."""
        debug(f"ChemicalPatterns.spherical: n={n}, radius={radius}")
        positions = generate_spherical_positions(n, radius)
        orientations = generate_face_center_orientations(positions, facing="outward")
        return positions, orientations
    
    @staticmethod
    def linear(n: int, distance: float = 5.0, axis: Tuple[float, float, float] = (0, 0, 1), **kwargs) -> Tuple[List[Position], List[Orientation]]:
        """Simple linear arrangement along defined axis with optional per-molecule rotations."""
        debug(f"ChemicalPatterns.linear: n={n}, distance={distance}, axis={axis}, kwargs={list(kwargs.keys())}")
        positions = generate_linear_positions(n, spacing=distance, axis=axis)
        
        # Check for per-molecule rotations in kwargs
        rotations = kwargs.get('rotations')
        if rotations and len(rotations) > 0:
            debug(f"  Using per-molecule rotations: {rotations}")
            orientations = []
            for i in range(n):
                if i < len(rotations):
                    rot = rotations[i]
                    # rot is dict like {'x': 0, 'y': 0, 'z': 90}
                    # Orientation takes (roll, pitch, yaw) = (x, y, z)
                    orientations.append(Orientation(
                        rot.get('x', 0), 
                        rot.get('y', 0), 
                        rot.get('z', 0)
                    ))
                else:
                    orientations.append(Orientation(0, 0, 0))
        else:
            orientations = generate_fixed_orientations(n)
        
        return positions, orientations
    
    @staticmethod
    def swastika(arm_length: float = 10.0, **kwargs) -> Tuple[List[Position], List[Orientation]]:
        """Four molecules in swastika arrangement."""
        debug(f"ChemicalPatterns.swastika: arm_length={arm_length}")
        half = arm_length / 2
        positions = [
            Position(half, 0, 0),
            Position(0, half, 0),
            Position(-half, 0, 0),
            Position(0, -half, 0)
        ]
        orientations = [
            Orientation(0, 0, 0),
            Orientation(0, 0, 90),
            Orientation(0, 0, 180),
            Orientation(0, 0, 270)
        ]
        return positions, orientations


# =============================================================================
# SECTION 7: Main API
# =============================================================================

def _normalize_pattern_name(pattern: str) -> str:
    """Normalize pattern name to canonical form."""
    normalized = pattern.lower().replace('-', '_').replace(' ', '_')
    return PATTERN_ALIASES.get(normalized, normalized)


def _get_pattern_generator(pattern: str) -> Callable:
    """Get the pattern generator function."""
    pattern_map = {
        'pi_pi_parallel': ChemicalPatterns.pi_pi_parallel,
        'pi_pi_antiparallel': ChemicalPatterns.pi_pi_antiparallel,
        'pi_pi_offset': ChemicalPatterns.pi_pi_offset,
        't_shaped': ChemicalPatterns.t_shaped,
        'herringbone': ChemicalPatterns.herringbone,
        'h_bonded': ChemicalPatterns.h_bonded_dimer,
        'h_bonded_circular': ChemicalPatterns.h_bonded_circular,
        'circular': ChemicalPatterns.h_bonded_circular,
        'helical': ChemicalPatterns.helical,
        'sandwich': ChemicalPatterns.sandwich,
        'grid': ChemicalPatterns.grid,
        'spherical': ChemicalPatterns.spherical,
        'linear': ChemicalPatterns.linear,
        'swastika': ChemicalPatterns.swastika,
    }
    generator = pattern_map.get(pattern)
    if generator is None:
        raise ValueError(f"Unknown pattern: '{pattern}'. Available: {list(pattern_map.keys())}")
    return generator


def auto_detect_pattern(molecules: List[Dict]) -> str:
    """Auto-detect best arrangement pattern based on molecule properties."""
    debug(f"auto_detect_pattern: n_molecules={len(molecules)}")
    n = len(molecules)
    
    def has_aromatic(mol):
        atoms = mol.get('atoms', [])
        c_count = atoms.count('C') if isinstance(atoms, list) else 0
        return c_count >= 6
    
    def can_hbond(mol):
        atoms = mol.get('atoms', [])
        has_donor_acceptor = any(a in ('O', 'N') for a in atoms)
        has_h = 'H' in atoms
        return has_donor_acceptor and has_h
    
    all_aromatic = all(has_aromatic(m) for m in molecules)
    all_hbond = all(can_hbond(m) for m in molecules)
    
    if n == 2 and all_aromatic:
        pattern = 'pi_pi_parallel'
    elif n == 2:
        pattern = 'linear'
    elif all_hbond and n <= 8:
        pattern = 'h_bonded_circular'
    elif all_aromatic:
        pattern = 'pi_pi_parallel'
    else:
        pattern = 'linear'
    
    debug(f"  auto-detected pattern: {pattern}")
    return pattern


def _align_planar_molecule_to_axis(mol: Dict[str, Any], target_axis: np.ndarray) -> Dict[str, Any]:
    """
    Align a planar molecule so its plane is PERPENDICULAR to target_axis.
    
    This is the GENERIC version of alignment that works for ANY stacking axis.
    For stacking along axis A, planar molecules should have their plane normal
    parallel to A (so the plane is perpendicular to the stacking direction).
    
    Args:
        mol: Molecule dict with 'atoms' and 'coords'
        target_axis: Normalized axis vector (the stacking direction)
    
    Returns:
        New molecule dict with aligned coordinates
    """
    coords = np.array(mol.get('coords', []))
    if len(coords) < 3:
        return mol  # Can't align with fewer than 3 atoms
    
    # Normalize target axis
    target_axis = np.array(target_axis, dtype=float)
    norm = np.linalg.norm(target_axis)
    if norm < 1e-10:
        target_axis = np.array([0, 0, 1])
    else:
        target_axis = target_axis / norm
    
    # Center the molecule
    center = coords.mean(axis=0)
    centered = coords - center
    
    # SVD to find plane normal (last singular vector = direction of min variance)
    try:
        _, S, Vt = np.linalg.svd(centered)
    except np.linalg.LinAlgError:
        debug("  WARNING: SVD failed in _align_planar_molecule_to_axis")
        return mol
    
    # Check planarity: if 3rd singular value is small compared to others, it's planar
    planarity_ratio = S[2] / S[0] if S[0] > 1e-10 else 1.0
    if planarity_ratio > 0.3:
        debug(f"  Molecule is not planar (ratio={planarity_ratio:.3f}), skipping alignment")
        return mol  # Not planar enough, don't align
    
    # Current plane normal
    current_normal = Vt[2]  # Last row = min variance direction
    
    # Ensure consistent normal direction (prefer alignment with target)
    if np.dot(current_normal, target_axis) < 0:
        current_normal = -current_normal
    
    # Check if already aligned
    alignment = np.abs(np.dot(current_normal, target_axis))
    if alignment > 0.999:
        debug(f"  Molecule already aligned to axis (alignment={alignment:.4f})")
        return mol  # Already aligned
    
    # Rodrigues' rotation formula: R = I + sin(θ)K + (1-cos(θ))K²
    axis = np.cross(current_normal, target_axis)
    axis_norm = np.linalg.norm(axis)
    if axis_norm < 1e-10:
        # Vectors are parallel or anti-parallel
        if np.dot(current_normal, target_axis) < 0:
            # Anti-parallel: rotate 180° around any perpendicular axis
            perp = np.array([1, 0, 0]) if abs(target_axis[0]) < 0.9 else np.array([0, 1, 0])
            perp = perp - np.dot(perp, target_axis) * target_axis
            perp = perp / np.linalg.norm(perp)
            R = 2 * np.outer(perp, perp) - np.eye(3)
        else:
            return mol
    else:
        axis = axis / axis_norm
        angle = np.arccos(np.clip(np.dot(current_normal, target_axis), -1.0, 1.0))
        
        # Skew-symmetric matrix K
        K = np.array([
            [0, -axis[2], axis[1]],
            [axis[2], 0, -axis[0]],
            [-axis[1], axis[0], 0]
        ])
        
        R = np.eye(3) + np.sin(angle) * K + (1 - np.cos(angle)) * (K @ K)
    
    # Apply rotation
    rotated = (R @ centered.T).T
    
    new_mol = mol.copy()
    new_mol['coords'] = rotated.tolist()
    
    rotation_angle = np.degrees(np.arccos(np.clip(alignment, -1.0, 1.0)))
    debug(f"  Aligned planar molecule to axis {tuple(target_axis)}: angle={rotation_angle:.1f}°")
    
    return new_mol


def _align_planar_molecule_to_xy(mol: Dict[str, Any]) -> Dict[str, Any]:
    """
    Align a planar molecule so its plane lies in XY (normal along Z).

    This is critical for pi-stacking: molecules must be flat in XY plane
    before being stacked along Z-axis. Without this, molecules from
    database with arbitrary orientations cause atom clashes.

    Uses SVD to find the plane normal (direction of minimum variance),
    then applies Rodrigues' rotation to align it with Z-axis.

    Args:
        mol: Molecule dict with 'atoms' and 'coords'

    Returns:
        New molecule dict with aligned coordinates
    """
    coords = np.array(mol.get('coords', []))
    if len(coords) < 3:
        return mol  # Can't align with fewer than 3 atoms

    # Center the molecule
    center = coords.mean(axis=0)
    centered = coords - center

    # SVD to find plane normal (last singular vector = direction of min variance)
    try:
        _, S, Vt = np.linalg.svd(centered)
    except np.linalg.LinAlgError:
        debug("  WARNING: SVD failed in _align_planar_molecule_to_xy")
        return mol

    # Check planarity: if 3rd singular value is small compared to others, it's planar
    planarity_ratio = S[2] / S[0] if S[0] > 1e-10 else 1.0
    if planarity_ratio > 0.3:
        debug(f"  Molecule is not planar (planarity_ratio={planarity_ratio:.3f}), skipping alignment")
        return mol  # Not planar enough, don't align

    # Current plane normal
    current_normal = Vt[2]  # Last row = min variance direction

    # Ensure consistent normal direction (prefer positive Z component)
    if current_normal[2] < 0:
        current_normal = -current_normal

    # Target: Z-axis
    target_normal = np.array([0.0, 0.0, 1.0])

    # Check if already aligned
    alignment = np.abs(np.dot(current_normal, target_normal))
    if alignment > 0.999:
        debug(f"  Molecule already aligned to XY (alignment={alignment:.4f})")
        return mol  # Already aligned

    # Rodrigues' rotation formula: R = I + sin(θ)K + (1-cos(θ))K²
    # where K is the skew-symmetric cross-product matrix
    axis = np.cross(current_normal, target_normal)
    axis_norm = np.linalg.norm(axis)
    if axis_norm < 1e-10:
        # Vectors are parallel or anti-parallel
        if np.dot(current_normal, target_normal) < 0:
            # Anti-parallel: rotate 180° around X-axis
            R = np.array([[1, 0, 0], [0, -1, 0], [0, 0, -1]])
        else:
            # Parallel: no rotation needed
            return mol
    else:
        axis = axis / axis_norm  # Normalize rotation axis
        angle = np.arccos(np.clip(np.dot(current_normal, target_normal), -1.0, 1.0))

        # Skew-symmetric matrix K
        K = np.array([
            [0, -axis[2], axis[1]],
            [axis[2], 0, -axis[0]],
            [-axis[1], axis[0], 0]
        ])

        # Rodrigues' formula
        R = np.eye(3) + np.sin(angle) * K + (1 - np.cos(angle)) * (K @ K)

    # Apply rotation (coordinates stay centered at origin after rotation)
    rotated = (R @ centered.T).T

    # Translate back so center is at origin (for subsequent stacking)
    # We keep it centered - the pose system will translate it
    new_mol = mol.copy()
    new_mol['coords'] = rotated.tolist()

    # Compute angle for debug output
    rotation_angle = np.degrees(np.arccos(np.clip(alignment, -1.0, 1.0)))
    debug(f"  Aligned planar molecule: angle={rotation_angle:.1f}°, "
          f"normal ({current_normal[0]:.3f}, {current_normal[1]:.3f}, {current_normal[2]:.3f}) → (0, 0, 1)")

    return new_mol


def _prealign_molecules_for_pattern(molecules: List[Dict[str, Any]],
                                    pattern: str) -> List[Dict[str, Any]]:
    """
    Pre-align molecules for patterns that require specific orientations.

    For pi-stacking patterns, ensures all planar molecules are flat in XY plane.
    This must happen BEFORE pose generation and application.

    Args:
        molecules: List of molecule dicts
        pattern: The arrangement pattern name

    Returns:
        List of pre-aligned molecule dicts
    """
    # Patterns that require planar alignment (molecules flat in XY)
    # Includes grid since planar molecules in 2D grid should be flat
    planar_patterns = {
        'pi_pi_parallel', 'pi_pi_antiparallel', 'pi_pi_offset',
        't_shaped', 'herringbone', 'sandwich', 'grid'
    }

    if pattern not in planar_patterns:
        return molecules  # No pre-alignment needed

    debug(f"Pre-aligning molecules for {pattern} pattern")
    aligned = []
    for mol in molecules:
        aligned_mol = _align_planar_molecule_to_xy(mol)
        aligned.append(aligned_mol)

    return aligned


def apply_poses_to_molecules(molecules: List[Dict],
                             poses: List[MoleculePose]) -> List[Dict]:
    """Apply poses to molecules, transforming their coordinates."""
    debug(f"apply_poses_to_molecules: n_molecules={len(molecules)}, n_poses={len(poses)}")
    result = []

    for i, mol in enumerate(molecules):
        coords = np.array(mol.get('coords', []))
        debug(f"  molecule[{i}]: n_atoms={len(coords)}")
        
        if len(coords) == 0:
            debug(f"    WARNING: molecule {i} has no coordinates")
            result.append(mol)
            continue
        
        if i >= len(poses):
            debug(f"    WARNING: no pose for molecule {i}, using original")
            result.append(mol)
            continue
        
        pose = poses[i]
        debug(f"    pose: position=({pose.position.x:.3f}, {pose.position.y:.3f}, {pose.position.z:.3f})")
        debug(f"    pose: orientation=(roll={pose.orientation.roll:.1f}, pitch={pose.orientation.pitch:.1f}, yaw={pose.orientation.yaw:.1f})")
        
        # Center molecule at origin
        center = coords.mean(axis=0)
        centered = coords - center
        
        # Apply rotation
        R = pose.orientation.to_matrix()
        rotated = (R @ centered.T).T
        
        # Apply translation
        translated = rotated + pose.position.to_array()
        
        new_mol = mol.copy()
        new_mol['coords'] = translated.tolist()
        result.append(new_mol)
        debug(f"    transformed: new_center=({translated.mean(axis=0)[0]:.3f}, {translated.mean(axis=0)[1]:.3f}, {translated.mean(axis=0)[2]:.3f})")
    
    return result


def validate_arrangement(molecules: List[Dict], poses: List[MoleculePose],
                        min_distance: float = 1.5, constraints: Optional[List['Constraint']] = None) -> Dict[str, Any]:
    """
    Validate arrangement for clashes and physical plausibility.

    Returns:
        Dict with keys:
        - valid: bool - True if no errors
        - warnings: List[str] - Non-fatal issues
        - errors: List[str] - Fatal issues
        - clash_pairs: List[Tuple] - Clashing atom pairs
        - physics_violation: bool - True if atoms overlap (< 0.5Å)
    """
    debug(f"validate_arrangement: n_molecules={len(molecules)}, min_distance={min_distance}")

    transformed = apply_poses_to_molecules(molecules, poses)

    warnings: List[str] = []
    errors: List[str] = []
    clash_pairs: List[Tuple] = []
    physics_violation = False
    
    n = len(transformed)
    for i in range(n):
        for j in range(i + 1, n):
            coords_i = np.array(transformed[i].get('coords', []))
            coords_j = np.array(transformed[j].get('coords', []))
            atoms_i = transformed[i].get('atoms', [])
            atoms_j = transformed[j].get('atoms', [])
            
            if len(coords_i) == 0 or len(coords_j) == 0:
                continue
            
            # Check for clashes
            for ai, ci in enumerate(coords_i):
                r_i = VDW_RADII.get(atoms_i[ai] if ai < len(atoms_i) else 'C', 1.7)
                for aj, cj in enumerate(coords_j):
                    r_j = VDW_RADII.get(atoms_j[aj] if aj < len(atoms_j) else 'C', 1.7)
                    dist = np.linalg.norm(np.array(ci) - np.array(cj))
                    min_allowed = (r_i + r_j) * 0.7  # 70% of vdW sum
                    
                    if dist < 0.5:
                        # CRITICAL ERROR: Physics violation (nuclear fusion)
                        debug(f"    CRITICAL: Physics violation between mol {i} atom {ai} and mol {j} atom {aj}: dist={dist:.4f}A")
                        physics_violation = True
                        errors.append(f"Physics violation: Atoms overlapping at {dist:.4f}A between mol {i} atom {ai} and mol {j} atom {aj}")
                        # Return immediately for severe violations - no point checking more
                        if dist < 0.1:
                            return {
                                'valid': False,
                                'warnings': warnings,
                                'errors': errors,
                                'clash_pairs': clash_pairs,
                                'physics_violation': True
                            }

                    
                    if dist < min_allowed:
                        clash_pairs.append((i, j, ai, aj, dist))
                        if dist < min_distance:
                            errors.append(f"Severe clash: mol{i}:atom{ai} - mol{j}:atom{aj} = {dist:.2f}Å")
                        else:
                            warnings.append(f"Close contact: mol{i}:atom{ai} - mol{j}:atom{aj} = {dist:.2f}Å")
                            
    # Generic Constraint Validation
    if constraints:
        debug(f"Checking {len(constraints)} constraints...")
        for i, constraint in enumerate(constraints):
            violation = constraint.evaluate(molecules, poses)
            # Use a generous tolerance (0.5A) to allow for minor optimization residuals
            # but flag significant deviations.
            TOLERANCE = 0.5
            if violation > TOLERANCE:
                debug(f"  Constraint {i} violation: {violation:.4f}")
                errors.append(f"Constraint Violation: {type(constraint).__name__} mismatch by {violation:.2f} (tol={TOLERANCE})")

    # Valid only if no errors AND no physics violations
    valid = len(errors) == 0 and not physics_violation
    debug(f"  validation result: valid={valid}, physics_violation={physics_violation}, warnings={len(warnings)}, errors={len(errors)}")

    return {
        'valid': valid,
        'warnings': warnings,
        'errors': errors,
        'clash_pairs': clash_pairs,
        'physics_violation': physics_violation
    }


def combine_molecules(molecules: List[Dict], vacuum: float = 10.0) -> Dict[str, Any]:
    """Combine multiple molecules into a single structure."""
    debug(f"combine_molecules: n_molecules={len(molecules)}, vacuum={vacuum}")
    
    all_atoms = []
    all_coords = []
    
    for mol in molecules:
        atoms = mol.get('atoms', [])
        coords = mol.get('coords', [])
        all_atoms.extend(atoms)
        all_coords.extend(coords)
    
    all_coords = np.array(all_coords) if all_coords else np.array([])
    
    # Calculate bounding box for lattice
    if len(all_coords) > 0:
        min_coords = all_coords.min(axis=0)
        max_coords = all_coords.max(axis=0)
        size = max_coords - min_coords + 2 * vacuum
    else:
        size = np.array([vacuum, vacuum, vacuum])
    
    debug(f"  total atoms: {len(all_atoms)}")
    debug(f"  bounding box size: ({size[0]:.2f}, {size[1]:.2f}, {size[2]:.2f})")
    
    return {
        'atoms': all_atoms,
        'coords': all_coords.tolist() if len(all_coords) > 0 else [],
        'lattice': {
            'a': float(size[0]),
            'b': float(size[1]),
            'c': float(size[2]),
            'alpha': 90.0,
            'beta': 90.0,
            'gamma': 90.0
        }
    }


def _parse_axis(axis: Union[str, Tuple[float, float, float], List[float]]) -> Tuple[float, float, float]:
    """
    Convert axis specification to normalized vector tuple.
    
    Supports:
    - String: 'x', 'y', 'z' (case-insensitive)
    - Tuple/List: [x, y, z] or (x, y, z)
    - Normalizes non-unit vectors
    - Defaults to Z-axis for invalid/zero inputs
    """
    # String axis: 'x', 'y', 'z'
    if isinstance(axis, str):
        axis_map = {'x': (1.0, 0.0, 0.0), 'y': (0.0, 1.0, 0.0), 'z': (0.0, 0.0, 1.0)}
        return axis_map.get(axis.lower().strip(), (0.0, 0.0, 1.0))
    
    # Vector: [x, y, z] or (x, y, z)
    if isinstance(axis, (list, tuple)) and len(axis) >= 3:
        # Convert to floats explicitly (no try/except)
        if all(isinstance(v, (int, float)) for v in axis[:3]):
            x, y, z = float(axis[0]), float(axis[1]), float(axis[2])
            # Calculate magnitude
            mag_sq = x*x + y*y + z*z
            if mag_sq < 1e-10:  # Zero vector → default to Z
                return (0.0, 0.0, 1.0)
            # Normalize
            mag = mag_sq ** 0.5
            return (x/mag, y/mag, z/mag)
    
    # Invalid input → default to Z
    return (0.0, 0.0, 1.0)




def arrange_molecules(
    molecules: List[Dict[str, Any]],
    pattern: str = "auto",
    distance: Optional[float] = None,
    constraints: Optional[List[str]] = None,
    optimize: bool = False,
    validate: bool = True,
    relative_poses: Optional[List[Dict]] = None,
    formulas: Optional[Dict[str, str]] = None,
    vacuum: float = 10.0,
    **kwargs
) -> Dict[str, Any]:
    """
    Arrange pre-generated molecules in a chemical pattern.
    
    This is the main entry point for molecular arrangement.
    
    Args:
        molecules: List of molecule dicts with 'atoms' and 'coords'
        pattern: Arrangement pattern ("auto", "pi_pi_parallel", "circular", etc.)
        distance: Intermolecular distance in Angstroms
        constraints: List of constraint strings (e.g., "distance(0:centroid(), 1:centroid(), 3.4)")
        optimize: Run constraint solver after initial placement
        validate: Check for clashes and return validation info
        relative_poses: List of relative pose specifications
        formulas: Dict with position formulas {"x": "...", "y": "...", "z": "..."}
        vacuum: Vacuum padding for combined structure
        **kwargs: Pattern-specific parameters (radius, pitch, turns, offset, tilt, etc.)
    
    Returns:
        {
            "success": True/False,
            "molecules": [...],       # Transformed molecule dicts
            "poses": [...],           # Pose info per molecule
            "validation": {...},      # Clash/distance checks
            "combined": {...},        # Single structure (atoms, coords, lattice)
            "structure": {...},       # MCP-compatible format
            "metadata": {...}         # Pattern info
        }
    """
    debug("=" * 70)
    debug("arrange_molecules called")
    debug(f"  n_molecules: {len(molecules)}")
    debug(f"  pattern: {pattern}")
    debug(f"  distance: {distance}")
    debug(f"  constraints: {constraints}")
    debug(f"  optimize: {optimize}")
    debug(f"  formulas: {formulas}")
    debug(f"  kwargs: {kwargs}")
    debug("=" * 70)
    
    n = len(molecules)
    
    # Validate input
    if n == 0:
        debug("ERROR: No molecules provided")
        return {
            'success': False,
            'error': {'code': 'NO_MOLECULES', 'message': 'No molecules provided'}
        }
    
    # Check molecules have required fields
    for i, mol in enumerate(molecules):
        if 'atoms' not in mol or 'coords' not in mol:
            debug(f"ERROR: Molecule {i} missing atoms or coords")
            return {
                'success': False,
                'error': {'code': 'INVALID_MOLECULE', 'message': f'Molecule {i} missing atoms or coords'}
            }
        debug(f"  molecule[{i}]: {len(mol['atoms'])} atoms, identifier={mol.get('identifier', 'unknown')}")
    
    # Normalize pattern name
    pattern_normalized = _normalize_pattern_name(pattern) if pattern != "auto" else "auto"
    debug(f"  normalized pattern: {pattern_normalized}")
    
    # Auto-detect if needed
    if pattern_normalized == "auto":
        # CRITICAL FIX: If axis is explicitly provided, use linear stacking
        # because specifying axis implies linear arrangement along that axis
        axis = kwargs.get('axis')
        if axis is not None and axis != 'z':  # 'z' is default, so explicit axis means user intent
            pattern_normalized = 'linear'
            debug(f"  axis='{axis}' provided, forcing pattern to 'linear'")
        else:
            pattern_normalized = auto_detect_pattern(molecules)
    
    # =========================================================================
    # ROBUST PARAMETER HANDLING
    # =========================================================================
    # Note: We rely on generating a safe distance later if distance is missing/small.
    # We do NOT infer distance from offset_x/y/z as that conflates lateral offset 
    # with longitudinal spacing.

    # Set default distance based on pattern
    if distance is None:
        default_distances = INTERACTION_DISTANCES.get(pattern_normalized)
        distance = default_distances[0] if default_distances else 3.4
        debug(f"  using default distance: {distance}")
        
        # For LINEAR stacking, calculate smart distance based on molecule size
        if pattern_normalized == 'linear' and molecules:
            import numpy as np
            max_extent = 0
            for mol in molecules:
                coords = mol.get('coords', [])
                if coords:
                    coords_arr = np.array(coords)
                    extent = np.max(coords_arr, axis=0) - np.min(coords_arr, axis=0)
                    max_extent = max(max_extent, np.max(extent))
            
            # Smart distance = molecule size + 2Å buffer
            if max_extent > 0:
                smart_distance = max_extent + 2.0
                if smart_distance > distance:
                    debug(f"  LINEAR: auto-calculated distance {smart_distance:.2f}Å based on molecule size (~{max_extent:.1f}Å)")
                    distance = smart_distance
    
    # =========================================================================
    # EXPLICIT DISTANCE VALIDATION (no try/except)
    # =========================================================================
    # Handle non-numeric distance
    if not isinstance(distance, (int, float)):
        distance = 3.4  # Safe default
        debug(f"  WARNING: Invalid distance type, using default 3.4Å")
    
    # Handle negative distance - use absolute value
    if distance < 0:
        debug(f"  INFO: Negative distance {distance}Å → using absolute value {abs(distance)}Å")
        distance = abs(distance)
    
    # Cap extreme molecule counts
    MAX_MOLECULES = 500
    if n > MAX_MOLECULES:
        debug(f"ERROR: n={n} exceeds MAX_MOLECULES={MAX_MOLECULES}")
        return {
            'success': False,
            'error': {
                'code': 'TOO_MANY_MOLECULES', 
                'message': f'Requested {n} molecules exceeds limit of {MAX_MOLECULES}. Please reduce count or use multiple arrangements.'
            }
        }
    # Optimization: Large System Heuristic
    if n > 50 and pattern == "auto" and not constraints:
        debug(f"Optimization: Large cluster (n={n}) detected, auto-selecting 'grid' pattern for efficiency")
        pattern_normalized = "grid"

    # PRE-ALIGNMENT: Align planar molecules for patterns that require it
    # This must happen BEFORE pose generation to ensure proper stacking
    molecules = _prealign_molecules_for_pattern(molecules, pattern_normalized)

    # =========================================================================
    # GENERIC DIRECTION PARSING (supports axis, direction_vector, direction_angle)
    # =========================================================================
    # Priority: direction_vector > direction_angle > axis > offset_x/y/z inference
    if 'direction_vector' in kwargs and kwargs['direction_vector'] is not None:
        # Direct vector: [x, y, z]
        vec = kwargs['direction_vector']
        parsed_axis = _parse_axis(vec)
        debug(f"  Direction from vector: {parsed_axis}")
    elif 'direction_angle' in kwargs and kwargs['direction_angle'] is not None:
        # Angle-based: azimuth + elevation
        azimuth = kwargs['direction_angle']
        elevation = kwargs.get('direction_elevation', 0.0)
        parsed_axis = direction_from_angles(azimuth, elevation)
        debug(f"  Direction from angles: azimuth={azimuth}°, elevation={elevation}° -> {parsed_axis}")
    elif 'axis' in kwargs and kwargs['axis'] is not None:
        # Explicit axis string: 'x', 'y', 'z'
        axis_param = kwargs['axis']
        parsed_axis = _parse_axis(axis_param)
        debug(f"  Direction from axis: {axis_param} -> {parsed_axis}")
    else:
        # Infer axis from offset_x/y/z if present (LLM sometimes uses these instead of axis)
        inferred_axis = 'z'  # Default
        if kwargs.get('offset_x', 0) != 0:
            inferred_axis = 'x'
            debug(f"  Inferred axis=x from offset_x={kwargs.get('offset_x')}")
        elif kwargs.get('offset_y', 0) != 0:
            inferred_axis = 'y'
            debug(f"  Inferred axis=y from offset_y={kwargs.get('offset_y')}")
        else:
            debug(f"  Using default axis=z (no axis specified)")
        parsed_axis = _parse_axis(inferred_axis)

    # Generate positions and orientations
    positions = []
    orientations = []
    position_info = {}  # Will be filled by SMART positioning if used
    
    if formulas:
        # Use formula-based generation
        debug("Using formula-based position generation")
        x_formula = formulas.get('x', '0')
        y_formula = formulas.get('y', '0')
        z_formula = formulas.get('z', f'{distance} * i')
        try:
            positions = generate_formula_positions(n, x_formula, y_formula, z_formula)
        except FormulaError as err:
            debug(f"FormulaError: {err}")
            return {
                'success': False,
                'error': {
                    'code': 'FORMULA_ERROR',
                    'message': str(err),
                    'formulas': formulas
                }
            }
        orientations = generate_fixed_orientations(n)
    else:
        # Use pattern-based generation
        debug(f"Using pattern-based generation: {pattern_normalized}")
        generator = _get_pattern_generator(pattern_normalized)
        
        # Build kwargs for generator
        gen_kwargs = kwargs.copy()
        gen_kwargs.update({'n': n, 'distance': distance, 'axis': parsed_axis})
        
        # Handle different generator signatures
        if pattern_normalized in ['t_shaped', 'h_bonded', 'sandwich']:
            # t_shaped is strictly a dimer - reject n > 2
            if pattern_normalized == 't_shaped' and n > 2:
                return {
                    'success': False,
                    'error': {
                        'code': 'INVALID_COUNT',
                        'message': f"t_shaped pattern is a dimer (exactly 2 molecules). Requested {n}. Use 'herringbone' for alternating tilts or 'linear' for simple stacking."
                    }
                }
            
            positions, orientations = generator(distance=distance)
            
            # Extend to n molecules if needed (for sandwich and h_bonded)
            if len(positions) < n:
                if pattern_normalized == 'sandwich':
                    # Sandwich is A-B-A pattern - extend with alternating orientations
                    for i in range(len(positions), n):
                        z_offset = distance * i
                        positions.append(Position(0, 0, z_offset))
                        # Alternate orientation for intercalated pattern
                        orientations.append(Orientation(0, 0, (i % 2) * 180))
                else:
                    # h_bonded: extend as pairs maintaining H-bond geometry
                    base_count = len(positions)
                    for i in range(base_count, n):
                        pair_idx = i % base_count
                        z_offset = distance * (i // base_count)
                        base_pos = positions[pair_idx]
                        positions.append(Position(base_pos.x, base_pos.y, base_pos.z + z_offset * base_count))
                        orientations.append(orientations[pair_idx])
        elif pattern_normalized in ['swastika']:
            positions, orientations = generator(arm_length=kwargs.get('arm_length', 10.0))
        elif pattern_normalized in ['circular', 'h_bonded_circular']:
            # These use radius, not distance
            radius = kwargs.get('radius', distance * n / (2 * 3.14159))  # Approximate circle
            debug(f"  circular pattern: using radius={radius:.2f}")
            positions, orientations = generator(n=n, radius=radius)
        elif pattern_normalized in ['helical']:
            # Helical uses radius and pitch
            radius = kwargs.get('radius', 5.0)
            pitch = kwargs.get('pitch', distance)
            turns = kwargs.get('turns', 1.0)
            debug(f"  helical pattern: radius={radius}, pitch={pitch}, turns={turns}")
            positions, orientations = generator(n=n, radius=radius, pitch=pitch, turns=turns)
        elif pattern_normalized in ['spherical']:
            # Spherical uses radius
            radius = kwargs.get('radius', distance * 2)
            debug(f"  spherical pattern: using radius={radius:.2f}")
            positions, orientations = generator(n=n, radius=radius)
        elif pattern_normalized in ['grid']:
            # Grid uses spacing
            spacing = kwargs.get('spacing', distance)
            debug(f"  grid pattern: using spacing={spacing:.2f}")
            positions, orientations = generator(n=n, spacing=spacing)
        else:
            # =========================================================
            # CENTER-TO-CENTER POSITIONING (Chemical convention)
            # =========================================================
            # When user says "separated by 3.4Å", they mean center-to-center
            # distance, which is the standard chemical convention for
            # intermolecular distances (especially in pi-stacking).
            # =========================================================
            if pattern_normalized in GAP_BASED_PATTERNS:
                debug(f"  Using CENTER-TO-CENTER positioning for {pattern_normalized}")
                
                # Generate orientations first (pattern-dependent)
                user_rotations = kwargs.get('rotations')
                if user_rotations and isinstance(user_rotations, list) and len(user_rotations) > 0:
                    debug(f"  Using {len(user_rotations)} explicit per-molecule rotations from input")
                    orientations = []
                    for i in range(n):
                        rot_data = user_rotations[i] if i < len(user_rotations) else {}
                        # Handle various rotation formats
                        if isinstance(rot_data, dict):
                            rx = rot_data.get('x', 0)
                            ry = rot_data.get('y', 0)
                            rz = rot_data.get('z', 0)
                            orientations.append(Orientation(rx, ry, rz))
                        else:
                            orientations.append(Orientation(0, 0, 0))
                elif pattern_normalized == 'pi_pi_antiparallel':
                    orientations = generate_alternating_orientations(n, (0, 0, 0), (0, 0, 180))
                elif pattern_normalized == 'herringbone':
                    tilt = kwargs.get('tilt', 45.0)
                    orientations = []
                    for i in range(n):
                        if i % 2 == 0:
                            orientations.append(Orientation(tilt, 0, 0))
                        else:
                            orientations.append(Orientation(-tilt, 0, 0))
                elif pattern_normalized == 'pi_pi_offset':
                    # For offset, we use fixed orientations but the offset is in positions
                    orientations = generate_fixed_orientations(n)
                else:
                    orientations = generate_fixed_orientations(n)
                
                # =========================================================
                # SMART POSITION CALCULATION - Generic for ALL scenarios
                # =========================================================
                # Uses generate_smart_positions which:
                # 1. Works for ANY axis (X, Y, Z, or arbitrary vectors)
                # 2. Calculates molecular extent along the placement axis
                # 3. Uses safe distance if requested distance would cause clashes
                # 4. Returns position_info for user notification
                # =========================================================
                positions, position_info = generate_smart_positions(
                    molecules=molecules,
                    distance=distance,  # May be None/0 -> use safe distance
                    axis=parsed_axis,   # CRITICAL: Use parsed axis!
                    orientations=orientations,
                    allow_clashes=kwargs.get('allow_clashes', False)
                )
                
                if position_info.get('was_adjusted'):
                    debug(f"  ⚠ Distance adjusted: {position_info.get('adjustment_reason')}")
                
                # Handle lateral offset for pi_pi_offset pattern
                if pattern_normalized == 'pi_pi_offset':
                    offset = kwargs.get('offset', 1.5)
                    for i in range(len(positions)):
                        if i % 2 == 1:
                            positions[i] = Position(
                                positions[i].x + offset,
                                positions[i].y,
                                positions[i].z
                            )
            else:
                # Fallback to legacy generator for unknown patterns
                debug(f"  Using legacy generator for {pattern_normalized}")
                positions, orientations = generator(**gen_kwargs)
                position_info = {}
    
    # Trim to actual number of molecules
    positions = positions[:n]
    orientations = orientations[:n]
    
    debug(f"  generated {len(positions)} positions and {len(orientations)} orientations")
    
    # Create poses
    poses = []
    for i in range(n):
        pose = MoleculePose(
            position=positions[i] if i < len(positions) else Position(),
            orientation=orientations[i] if i < len(orientations) else Orientation(),
            molecule_idx=i,
            metadata={'pattern': pattern_normalized}
        )
        poses.append(pose)
    
    debug(f"  created {len(poses)} poses")
    
    # Handle relative poses if specified
    if relative_poses:
        debug("Processing relative poses")
        frames = [MolecularFrame.from_coords(np.array(mol.get('coords', [[0,0,0]]))) for mol in molecules]
        for rel_spec in relative_poses:
            idx = rel_spec.get('molecule_idx', 0)
            if idx < len(poses):
                poses[idx].parent_idx = rel_spec.get('parent_idx')
                if 'offset' in rel_spec:
                    poses[idx].local_offset = Position(*rel_spec['offset'])
                if 'orientation' in rel_spec:
                    poses[idx].local_orientation = Orientation(*rel_spec['orientation'])
        poses = resolve_relative_poses(poses, frames)
    
    # Build constraint objects
    constraint_objs = []
    if constraints:
        debug(f"Parsing {len(constraints)} constraints")
        for spec in constraints:
            constraint = _parse_constraint_string(spec)
            if constraint:
                constraint_objs.append(constraint)
                debug(f"  parsed constraint: {type(constraint).__name__}")
    
    # Add implicit pattern constraints
    if 'pi' in pattern_normalized:
        debug("Adding implicit plane alignment constraints for pi-stacking")
        for i in range(n - 1):
            constraint_objs.append(PlaneAlignmentConstraint(i, i + 1, "parallel"))
    
    # Add clash constraints
    debug("Adding clash constraints")
    for i in range(n):
        for j in range(i + 1, n):
            constraint_objs.append(ClashConstraint(i, j, priority=2))
    
    # Optimize if requested (but skip if explicit rotations provided to preserve user intent)
    explicit_rotations = kwargs.get('rotations') and len(kwargs.get('rotations')) > 0
    if optimize and constraint_objs and not explicit_rotations:
        debug(f"Running constraint solver with {len(constraint_objs)} constraints")
        solver = ConstraintSolver(molecules, poses, constraint_objs)
        poses = solver.solve(max_iterations=1000, max_time_seconds=30.0)
        debug("Constraint solver completed")
    
    # Apply poses to molecules
    debug("Applying poses to molecules")
    arranged_molecules = apply_poses_to_molecules(molecules, poses)
    
    # =========================================================================
    # SMART VALIDATION WITH AUTO-RELAXATION
    # =========================================================================
    validation = {}
    original_distance = distance
    relaxation_applied = False
    relaxation_attempts = 0
    max_relax_attempts = 5
    relax_factor = 1.2
    
    if validate:
        debug("Validating arrangement with auto-relaxation")
        
        for attempt in range(max_relax_attempts + 1):
            validation = validate_arrangement(molecules, poses, min_distance=1.5, constraints=constraint_objs)
            debug(f"  Attempt {attempt}: valid={validation.get('valid')}, physics_violation={validation.get('physics_violation')}")
            
            # Physics violation is unrecoverable (atoms at same position)
            if validation.get('physics_violation'):
                debug("  FATAL: Physics violation detected, cannot auto-relax")
                break
            
            # If valid, we're done
            if validation.get('valid'):
                break
            
            # If not valid and we can try again, auto-relax by expanding spacing
            if attempt < max_relax_attempts:
                relaxation_attempts = attempt + 1
                relaxation_applied = True
                distance = distance * relax_factor
                debug(f"  Auto-relaxing: increasing distance to {distance:.2f}Å (attempt {relaxation_attempts})")
                
                # Re-generate positions with increased gap
                if pattern_normalized in GAP_BASED_PATTERNS:
                    positions = generate_gap_based_positions(
                        molecules=molecules,
                        gaps=distance,
                        axis=parsed_axis,
                        orientations=orientations
                    )
                    # Update poses with new positions
                    for i, pos in enumerate(positions):
                        if i < len(poses):
                            poses[i].position = pos
                    
                    # Re-apply poses
                    arranged_molecules = apply_poses_to_molecules(molecules, poses)
        
        if relaxation_applied:
            debug(f"  Relaxation complete: {original_distance:.2f}Å -> {distance:.2f}Å ({relaxation_attempts} attempts)")
        
        # Return early with clear error for physics violations
        if validation.get('physics_violation'):
            errors = validation.get('errors', [])
            return {
                'success': False,
                'error': {
                    'code': 'PHYSICS_VIOLATION',
                    'message': f"Atoms overlapping: {errors[0] if errors else 'molecules placed at same position'}",
                    'details': errors
                },
                'validation': validation
            }

    # Combine into single structure
    debug("Combining molecules into single structure")
    combined = combine_molecules(arranged_molecules, vacuum)
    
    # Build MCP-compatible structure
    structure = {
        'atoms': combined['atoms'],
        'coords': combined['coords'],
        'lattice': combined['lattice'],
        'species': [{'element': a, 'occupation': 1.0} for a in combined['atoms']],
        'sites': [
            {'species': [{'element': combined['atoms'][i], 'occupation': 1.0}], 
             'abc': [0.5, 0.5, 0.5],  # Placeholder fractional coords
             'xyz': combined['coords'][i]}
            for i in range(len(combined['atoms']))
        ]
    }
    
    # Determine success
    success = validation.get('valid', True) if validate else True

    debug(f"arrange_molecules completed: success={success}")
    debug("=" * 70)

    result = {
        'success': success,
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
        'combined': combined,
        'structure': structure,
        'n_molecules': n,
        'n_atoms': len(combined['atoms']),
        'atoms': combined['atoms'],
        'metadata': {
            'pattern': pattern_normalized,
            'n_molecules': n,
            'position_info': position_info,
            'validation': validation,
            # Surface warnings to user
            'warnings': (
                ([position_info.get('adjustment_reason')] if position_info.get('was_adjusted') else []) +
                validation.get('warnings', [])
            )
        }
    }


    # Add error field when validation fails (BUG #2 fix)
    if not success and validate:
        errors = validation.get('errors', [])
        if validation.get('physics_violation'):
            result['error'] = {
                'code': 'PHYSICS_VIOLATION',
                'message': errors[0] if errors else 'Atoms overlapping',
                'details': errors
            }
        elif errors:
            result['error'] = {
                'code': 'VALIDATION_FAILED',
                'message': errors[0],
                'details': errors
            }
        else:
            result['error'] = {
                'code': 'VALIDATION_FAILED',
                'message': 'Arrangement validation failed',
                'details': validation.get('warnings', [])
            }

    return result


def _parse_constraint_string(s: str) -> Optional[Constraint]:
    """Parse constraint from string specification."""
    debug(f"_parse_constraint_string: '{s}'")
    
    match = re.match(r"(\w+)\s*\((.*)\)", s)
    if not match:
        debug(f"  failed to parse constraint string")
        return None
    
    func_name = match.group(1).lower()
    args_str = match.group(2)
    
    # Parse arguments
    args = []
    kwargs = {}
    current = []
    paren_depth = 0
    
    for char in args_str + ',':
        if char == ',' and paren_depth == 0:
            part = ''.join(current).strip()
            if part:
                if '=' in part and '(' not in part.split('=')[0]:
                    k, v = part.split('=', 1)
                    k = k.strip()
                    v = v.strip().strip('"\'')
                    # Try to convert to number
                    if v.replace('.', '').replace('-', '').isdigit():
                        kwargs[k] = float(v)
                    else:
                        kwargs[k] = v
                else:
                    args.append(part)
            current = []
        else:
            if char == '(':
                paren_depth += 1
            elif char == ')':
                paren_depth -= 1
            current.append(char)
    
    debug(f"  func_name={func_name}, args={args}, kwargs={kwargs}")
    
    if func_name == 'distance':
        return DistanceConstraint(
            selector1=args[0] if len(args) > 0 else kwargs.get('sel1', '0:centroid()'),
            selector2=args[1] if len(args) > 1 else kwargs.get('sel2', '1:centroid()'),
            target=float(args[2]) if len(args) > 2 else kwargs.get('target'),
            min_dist=kwargs.get('min'),
            max_dist=kwargs.get('max')
        )
    elif func_name == 'angle':
        return AngleConstraint(
            selector1=args[0] if len(args) > 0 else kwargs.get('sel1', '0:centroid()'),
            selector2=args[1] if len(args) > 1 else kwargs.get('sel2', '1:centroid()'),
            selector3=args[2] if len(args) > 2 else kwargs.get('sel3', '2:centroid()'),
            target=float(args[3]) if len(args) > 3 else kwargs.get('target', 180.0),
            weight=5.0  # Boost weight to ensure angular compliance
        )
    elif func_name in ['h_bond', 'hbond']:
        return HBondConstraint(
            donor_mol=int(args[0]) if len(args) > 0 else int(kwargs.get('donor', 0)),
            acceptor_mol=int(args[1]) if len(args) > 1 else int(kwargs.get('acceptor', 1))
        )
    elif func_name in ['plane_parallel', 'parallel']:
        return PlaneAlignmentConstraint(
            mol1=int(args[0]) if len(args) > 0 else int(kwargs.get('mol1', 0)),
            mol2=int(args[1]) if len(args) > 1 else int(kwargs.get('mol2', 1)),
            mode="parallel"
        )
    elif func_name in ['plane_perpendicular', 'perpendicular']:
        return PlaneAlignmentConstraint(
            mol1=int(args[0]) if len(args) > 0 else int(kwargs.get('mol1', 0)),
            mol2=int(args[1]) if len(args) > 1 else int(kwargs.get('mol2', 1)),
            mode="perpendicular"
        )
    
    debug(f"  unknown constraint type: {func_name}")
    return None


# =============================================================================
# MODULE TEST
# =============================================================================

if __name__ == "__main__":
    print("=" * 70)
    print("MOLECULAR ARRANGEMENT - Unified Engine Test")
    print("=" * 70)
    
    # Test molecule
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
    
    # Test 1: Pi-pi stacking
    print("\nTest 1: π-π parallel stacking of 3 benzenes")
    result = arrange_molecules(
        [benzene.copy() for _ in range(3)],
        pattern="pi_pi_parallel",
        distance=3.4
    )
    print(f"  Success: {result['success']}")
    print(f"  Molecules: {result['n_molecules']}")
    print(f"  Atoms: {result['n_atoms']}")
    print(f"  Valid: {result['validation'].get('valid', 'N/A')}")
    
    # Test 2: Circular
    print("\nTest 2: Circular arrangement of 6 benzenes")
    result = arrange_molecules(
        [benzene.copy() for _ in range(6)],
        pattern="circular",
        kwargs={'radius': 5.0}
    )
    print(f"  Success: {result['success']}")
    print(f"  Pattern: {result['metadata']['pattern']}")
    
    # Test 3: Formula
    print("\nTest 3: Formula-based spiral")
    result = arrange_molecules(
        [benzene.copy() for _ in range(5)],
        formulas={
            'x': '3 * cos(i * 1.5)',
            'y': '3 * sin(i * 1.5)',
            'z': 'i * 2'
        }
    )
    print(f"  Success: {result['success']}")
    for i, pose in enumerate(result['poses']):
        print(f"  Pose {i}: {pose['position']}")
    
    print("\n" + "=" * 70)
    print("All tests completed!")
    print("=" * 70)
