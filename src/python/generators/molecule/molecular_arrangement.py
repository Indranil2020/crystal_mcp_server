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
        return self.weight * abs(angle - self.target) / 180.0


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
                dyaw = step_size * grad_rot[2] * 10
                
                total_change += abs(dx) + abs(dy) + abs(dz) + abs(dyaw)
                
                self.poses[i].position.x -= dx
                self.poses[i].position.y -= dy
                self.poses[i].position.z -= dz
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
        """Compute numerical gradients."""
        n = len(self.poses)
        gradients = []
        delta = 0.01
        
        for i in range(n):
            grad_pos = np.zeros(3)
            grad_rot = np.zeros(3)
            
            for dim in range(3):
                # Position gradient
                orig = [self.poses[i].position.x, self.poses[i].position.y, self.poses[i].position.z][dim]
                
                if dim == 0:
                    self.poses[i].position.x = orig + delta
                elif dim == 1:
                    self.poses[i].position.y = orig + delta
                else:
                    self.poses[i].position.z = orig + delta
                
                v_plus = sum(c.evaluate(self.molecules, self.poses) for c in self.constraints)
                
                if dim == 0:
                    self.poses[i].position.x = orig - delta
                elif dim == 1:
                    self.poses[i].position.y = orig - delta
                else:
                    self.poses[i].position.z = orig - delta
                
                v_minus = sum(c.evaluate(self.molecules, self.poses) for c in self.constraints)
                
                if dim == 0:
                    self.poses[i].position.x = orig
                elif dim == 1:
                    self.poses[i].position.y = orig
                else:
                    self.poses[i].position.z = orig
                
                grad_pos[dim] = (v_plus - v_minus) / (2 * delta)
            
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


def generate_formula_positions(n: int, x_formula: str = "0", y_formula: str = "0",
                               z_formula: str = "i * 3.4") -> List[Position]:
    """Generate positions using mathematical formulas."""
    debug(f"generate_formula_positions: n={n}")
    debug(f"  x_formula = '{x_formula}'")
    debug(f"  y_formula = '{y_formula}'")
    debug(f"  z_formula = '{z_formula}'")
    
    positions = []
    safe_dict = {
        'sin': np.sin, 'cos': np.cos, 'tan': np.tan,
        'sqrt': np.sqrt, 'exp': np.exp, 'log': np.log,
        'abs': np.abs, 'pi': np.pi, 'e': np.e,
        'atan2': np.arctan2, 'floor': np.floor, 'ceil': np.ceil, 'round': np.round,
        'int': int
    }
    
    for i in range(n):
        t = i / max(n - 1, 1)  # Normalized [0, 1]
        local_vars = {'i': i, 'n': n, 't': t, **safe_dict}
        
        x = float(eval(x_formula, {"__builtins__": {}}, local_vars))
        y = float(eval(y_formula, {"__builtins__": {}}, local_vars))
        z = float(eval(z_formula, {"__builtins__": {}}, local_vars))
        
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
    def pi_pi_parallel(n: int, distance: float = 3.4) -> Tuple[List[Position], List[Orientation]]:
        """Face-to-face π-stacking with all molecules parallel."""
        debug(f"ChemicalPatterns.pi_pi_parallel: n={n}, distance={distance}")
        positions = generate_linear_positions(n, spacing=distance, axis=(0, 0, 1))
        orientations = generate_fixed_orientations(n)
        return positions, orientations
    
    @staticmethod
    def pi_pi_antiparallel(n: int, distance: float = 3.4) -> Tuple[List[Position], List[Orientation]]:
        """Face-to-face stacking with alternating 180° rotation."""
        debug(f"ChemicalPatterns.pi_pi_antiparallel: n={n}, distance={distance}")
        positions = generate_linear_positions(n, spacing=distance, axis=(0, 0, 1))
        orientations = generate_alternating_orientations(n, (0, 0, 0), (0, 0, 180))
        return positions, orientations
    
    @staticmethod
    def pi_pi_offset(n: int, distance: float = 3.4, 
                     offset: float = 1.5) -> Tuple[List[Position], List[Orientation]]:
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
    def t_shaped(distance: float = 5.0) -> Tuple[List[Position], List[Orientation]]:
        """Edge-to-face perpendicular dimer."""
        debug(f"ChemicalPatterns.t_shaped: distance={distance}")
        positions = [Position(0, 0, 0), Position(distance, 0, 0)]
        orientations = [Orientation(0, 0, 0), Orientation(0, 90, 0)]
        return positions, orientations
    
    @staticmethod
    def herringbone(n: int, distance: float = 5.5, 
                    tilt: float = 45.0) -> Tuple[List[Position], List[Orientation]]:
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
    def h_bonded_dimer(distance: float = 2.8) -> Tuple[List[Position], List[Orientation]]:
        """Hydrogen-bonded dimer with optimal geometry."""
        debug(f"ChemicalPatterns.h_bonded_dimer: distance={distance}")
        positions = [Position(0, 0, 0), Position(distance, 0, 0)]
        orientations = [Orientation(0, 0, 0), Orientation(0, 0, 180)]
        return positions, orientations
    
    @staticmethod
    def h_bonded_circular(n: int, radius: float = 5.0) -> Tuple[List[Position], List[Orientation]]:
        """Circular ring with H-bond connectivity."""
        debug(f"ChemicalPatterns.h_bonded_circular: n={n}, radius={radius}")
        positions = generate_circular_positions(n, radius)
        orientations = generate_face_center_orientations(positions, facing="inward")
        return positions, orientations
    
    @staticmethod
    def helical(n: int, radius: float = 5.0, pitch: float = 3.4,
                turns: float = 1.0) -> Tuple[List[Position], List[Orientation]]:
        """Helical (DNA-like) arrangement."""
        debug(f"ChemicalPatterns.helical: n={n}, radius={radius}, pitch={pitch}, turns={turns}")
        positions = generate_helical_positions(n, radius, pitch, turns)
        orientations = generate_face_center_orientations(positions, facing="inward")
        return positions, orientations
    
    @staticmethod
    def sandwich(distance: float = 3.4) -> Tuple[List[Position], List[Orientation]]:
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
    def grid(n: int, spacing: float = 3.4) -> Tuple[List[Position], List[Orientation]]:
        """2D grid arrangement in XY plane."""
        debug(f"ChemicalPatterns.grid: n={n}, spacing={spacing}")
        positions = generate_grid_positions(n, spacing)
        orientations = generate_fixed_orientations(n)
        return positions, orientations
    
    @staticmethod
    def spherical(n: int, radius: float = 5.0) -> Tuple[List[Position], List[Orientation]]:
        """Molecules distributed on sphere surface."""
        debug(f"ChemicalPatterns.spherical: n={n}, radius={radius}")
        positions = generate_spherical_positions(n, radius)
        orientations = generate_face_center_orientations(positions, facing="outward")
        return positions, orientations
    
    @staticmethod
    def linear(n: int, distance: float = 5.0) -> Tuple[List[Position], List[Orientation]]:
        """Simple linear arrangement along Z axis."""
        debug(f"ChemicalPatterns.linear: n={n}, distance={distance}")
        positions = generate_linear_positions(n, spacing=distance, axis=(0, 0, 1))
        orientations = generate_fixed_orientations(n)
        return positions, orientations
    
    @staticmethod
    def swastika(arm_length: float = 10.0) -> Tuple[List[Position], List[Orientation]]:
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
                        min_distance: float = 1.5) -> Dict[str, Any]:
    """Validate arrangement for clashes and physical plausibility."""
    debug(f"validate_arrangement: n_molecules={len(molecules)}, min_distance={min_distance}")
    
    transformed = apply_poses_to_molecules(molecules, poses)
    
    warnings = []
    warnings = []
    errors = []
    clash_pairs = []
    validation = {'valid': True, 'warnings': warnings, 'errors': errors}
    
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
                        validation['valid'] = False
                        if 'errors' not in validation: validation['errors'] = []
                        validation['errors'].append(f"Physics violation: Atom overlap at {dist:.2f}A (<0.5A)")
                        # In strict mode, we might raise here, but for now we mark invalid
                        # Actually, let's raise for blatant 0.0 case
                        if dist < 0.1:
                            raise ValueError(f"Physics violation: Atoms overlapping at {dist:.4f}A")

                    
                    if dist < min_allowed:
                        clash_pairs.append((i, j, ai, aj, dist))
                        if dist < min_distance:
                            errors.append(f"Severe clash: mol{i}:atom{ai} - mol{j}:atom{aj} = {dist:.2f}Å")
                        else:
                            warnings.append(f"Close contact: mol{i}:atom{ai} - mol{j}:atom{aj} = {dist:.2f}Å")
    
    valid = len(errors) == 0
    debug(f"  validation result: valid={valid}, warnings={len(warnings)}, errors={len(errors)}")
    
    return {
        'valid': valid,
        'warnings': warnings,
        'errors': errors,
        'clash_pairs': clash_pairs
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
        pattern_normalized = auto_detect_pattern(molecules)
    
    # Set default distance based on pattern
    if distance is None:
        default_distances = INTERACTION_DISTANCES.get(pattern_normalized)
        distance = default_distances[0] if default_distances else 3.4
        debug(f"  using default distance: {distance}")
    
    # Generate positions and orientations
    positions = []
    orientations = []
    
    if formulas:
        # Use formula-based generation
        debug("Using formula-based position generation")
        x_formula = formulas.get('x', '0')
        y_formula = formulas.get('y', '0')
        z_formula = formulas.get('z', f'{distance} * i')
        positions = generate_formula_positions(n, x_formula, y_formula, z_formula)
        orientations = generate_fixed_orientations(n)
    else:
        # Use pattern-based generation
        debug(f"Using pattern-based generation: {pattern_normalized}")
        generator = _get_pattern_generator(pattern_normalized)
        
        # Build kwargs for generator
        gen_kwargs = {'n': n, 'distance': distance}
        gen_kwargs.update(kwargs)
        
        # Handle different generator signatures
        if pattern_normalized in ['t_shaped', 'h_bonded', 'sandwich']:
            positions, orientations = generator(distance=distance)
            # Extend to n molecules if needed
            while len(positions) < n:
                last_pos = positions[-1]
                positions.append(Position(last_pos.x, last_pos.y, last_pos.z + distance))
                orientations.append(orientations[-1] if orientations else Orientation())
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
            # Standard generators (linear, pi_pi_*) take n and distance
            positions, orientations = generator(n=n, distance=distance, **{k: v for k, v in kwargs.items() if k not in ['n', 'distance']})
    
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
    
    # Optimize if requested
    if optimize and constraint_objs:
        debug(f"Running constraint solver with {len(constraint_objs)} constraints")
        solver = ConstraintSolver(molecules, poses, constraint_objs)
        poses = solver.solve(max_iterations=1000, max_time_seconds=30.0)
        debug("Constraint solver completed")
    
    # Apply poses to molecules
    debug("Applying poses to molecules")
    arranged_molecules = apply_poses_to_molecules(molecules, poses)
    
    # Validate
    validation = {}
    if validate:
        debug("Validating arrangement")
        validation = validate_arrangement(molecules, poses)
        debug(f"Validation: valid={validation.get('valid')}, warnings={len(validation.get('warnings', []))}")
    
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
    
    debug("arrange_molecules completed successfully")
    debug("=" * 70)
    
    return {
        'success': True,
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
        'coords': combined['coords'],
        'metadata': {
            'pattern': pattern_normalized,
            'distance': distance,
            'optimized': optimize
        }
    }


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
