/**
 * Crystal Structure Type Definitions
 *
 * This module defines all data structures for crystal structures,
 * lattices, atoms, symmetry information, and related types.
 */
/**
 * 3x3 matrix representing lattice vectors or transformations
 */
export type Matrix3x3 = readonly [
    readonly [number, number, number],
    readonly [number, number, number],
    readonly [number, number, number]
];
/**
 * 3D vector for positions, forces, etc.
 */
export type Vector3 = readonly [number, number, number];
/**
 * Crystal system classifications
 */
export type CrystalSystem = "triclinic" | "monoclinic" | "orthorhombic" | "tetragonal" | "trigonal" | "hexagonal" | "cubic";
/**
 * Lattice parameters in crystallographic convention
 */
export interface LatticeParameters {
    readonly a: number;
    readonly b: number;
    readonly c: number;
    readonly alpha: number;
    readonly beta: number;
    readonly gamma: number;
}
/**
 * Complete lattice information
 */
export interface Lattice extends LatticeParameters {
    readonly matrix: Matrix3x3;
    readonly volume: number;
    readonly reciprocal_matrix?: Matrix3x3;
}
/**
 * Space group information
 */
export interface SpaceGroup {
    readonly number: number;
    readonly symbol: string;
    readonly hall_symbol?: string;
    readonly point_group: string;
    readonly crystal_system: CrystalSystem;
    readonly international_symbol?: string;
}
/**
 * Wyckoff position information
 */
export interface WyckoffPosition {
    readonly letter: string;
    readonly multiplicity: number;
    readonly site_symmetry: string;
    readonly coordinates?: Vector3;
}
/**
 * Atomic site information
 */
export interface AtomicSite {
    readonly element: string;
    readonly coords: Vector3;
    readonly cartesian: Vector3;
    readonly wyckoff: string;
    readonly multiplicity: number;
    readonly site_symmetry: string;
    readonly occupancy?: number;
    readonly charge?: number;
}
/**
 * Crystal structure metadata
 */
export interface CrystalMetadata {
    readonly formula: string;
    readonly natoms: number;
    readonly volume: number;
    readonly density: number;
    readonly packing_fraction?: number;
    readonly formation_energy?: number;
}
/**
 * Complete crystal structure
 */
export interface CrystalStructure {
    readonly lattice: Lattice;
    readonly atoms: readonly AtomicSite[];
    readonly space_group: SpaceGroup;
    readonly metadata: CrystalMetadata;
    readonly symmetry_operations?: {
        readonly rotations: readonly Matrix3x3[];
        readonly translations: readonly Vector3[];
        readonly n_operations: number;
    };
}
/**
 * Input parameters for crystal generation
 */
export interface CrystalGenerationParams {
    readonly composition: readonly string[];
    readonly space_group: number | string;
    readonly num_atoms?: number;
    readonly lattice_params?: Partial<LatticeParameters>;
    readonly volume_factor?: number;
    readonly min_distance?: Record<string, number>;
    readonly wyckoff_positions?: readonly {
        readonly element: string;
        readonly wyckoff: string;
        readonly coords?: Vector3;
    }[];
    readonly seed?: number;
    readonly max_attempts?: number;
}
/**
 * Validation issues found in a structure
 */
export interface StructureValidation {
    readonly valid: boolean;
    readonly issues: readonly string[];
    readonly warnings: readonly string[];
    readonly metrics?: {
        readonly min_distance: number;
        readonly avg_distance: number;
        readonly min_angle?: number;
        readonly max_angle?: number;
    };
}
/**
 * Structure generation result with validation
 */
export interface StructureGenerationResult {
    readonly structure: CrystalStructure;
    readonly validation: StructureValidation;
    readonly generation_metadata: {
        readonly attempts: number;
        readonly seed?: number;
        readonly generation_time_ms: number;
    };
    readonly files?: Record<string, string>;
    readonly file_paths?: Record<string, string>;
}
/**
 * Supercell transformation parameters
 */
export interface SupercellParams {
    readonly structure: CrystalStructure;
    readonly scaling_matrix: Matrix3x3 | Vector3;
    readonly wrap_atoms?: boolean;
    readonly preserve_symmetry?: boolean;
    readonly min_distance_check?: boolean;
}
/**
 * Slab generation parameters
 */
export interface SlabParams {
    readonly structure: CrystalStructure;
    readonly miller_indices: Vector3;
    readonly thickness: number;
    readonly vacuum: number;
    readonly center_slab?: boolean;
    readonly fix_bottom_layers?: number;
    readonly fix_atoms?: readonly number[];
    readonly orthogonalize?: boolean;
    readonly min_slab_size?: number;
    readonly symmetric?: boolean;
}
/**
 * MLFF optimization parameters
 */
export interface MLFFOptimizationParams {
    readonly structure: CrystalStructure;
    readonly mlff_model: string;
    readonly optimizer: string;
    readonly fmax: number;
    readonly steps: number;
    readonly constrain_symmetry?: boolean;
    readonly fix_lattice?: boolean;
    readonly fix_volume?: boolean;
    readonly pressure?: number;
    readonly trajectory_file?: string;
}
/**
 * MLFF optimization result
 */
export interface MLFFOptimizationResult {
    readonly optimized_structure: CrystalStructure;
    readonly initial_energy: number;
    readonly final_energy: number;
    readonly energy_change: number;
    readonly max_force_initial: number;
    readonly max_force_final: number;
    readonly n_steps: number;
    readonly converged: boolean;
    readonly preserved_symmetry: boolean;
    readonly final_space_group?: number;
    readonly trajectory?: readonly {
        readonly step: number;
        readonly energy: number;
        readonly max_force: number;
        readonly structure: CrystalStructure;
    }[];
    readonly timing: {
        readonly setup_time_ms: number;
        readonly optimization_time_ms: number;
        readonly total_time_ms: number;
    };
}
/**
 * Energy calculation result
 */
export interface EnergyCalculationResult {
    readonly energy: number;
    readonly energy_per_atom: number;
    readonly forces?: readonly Vector3[];
    readonly stress?: Matrix3x3;
    readonly metadata: {
        readonly model: string;
        readonly model_version?: string;
        readonly calculation_time_ms: number;
    };
}
/**
 * Ground state search parameters
 */
export interface GroundStateSearchParams {
    readonly composition: readonly string[];
    readonly space_groups?: readonly number[];
    readonly num_structures_per_group: number;
    readonly mlff_model: string;
    readonly optimization_settings: {
        readonly optimizer: string;
        readonly fmax: number;
        readonly steps: number;
        readonly constrain_symmetry: boolean;
    };
    readonly parallel: boolean;
    readonly n_workers?: number;
    readonly output_directory?: string;
    readonly save_trajectories?: boolean;
}
/**
 * Ground state search result
 */
export interface GroundStateSearchResult {
    readonly ground_state: {
        readonly structure: CrystalStructure;
        readonly energy: number;
        readonly space_group: number;
        readonly formation_energy?: number;
    };
    readonly all_results: readonly {
        readonly space_group: number;
        readonly attempt: number;
        readonly initial_structure: CrystalStructure;
        readonly optimized_structure: CrystalStructure;
        readonly initial_energy: number;
        readonly final_energy: number;
        readonly converged: boolean;
        readonly preserved_symmetry: boolean;
    }[];
    readonly energy_ranking: readonly {
        readonly rank: number;
        readonly space_group: number;
        readonly energy: number;
        readonly energy_above_ground_state: number;
    }[];
    readonly statistics: {
        readonly total_structures: number;
        readonly converged_structures: number;
        readonly unique_space_groups: number;
        readonly energy_range: readonly [number, number];
        readonly most_stable_crystal_system: CrystalSystem;
    };
}
/**
 * File export formats
 */
export type FileFormat = "cif" | "poscar" | "contcar" | "xyz" | "extxyz" | "json" | "pdb" | "pwscf" | "castep" | "lammps";
/**
 * Export structure parameters
 */
export interface ExportParams {
    readonly structure: CrystalStructure;
    readonly formats: readonly FileFormat[];
    readonly include_metadata?: boolean;
    readonly dft_software?: "vasp" | "qe" | "cp2k" | "castep";
    readonly output_directory?: string;
}
/**
 * Export result
 */
export interface ExportResult {
    readonly files: Record<FileFormat, string>;
    readonly file_paths?: Record<FileFormat, string>;
}
//# sourceMappingURL=crystal.d.ts.map