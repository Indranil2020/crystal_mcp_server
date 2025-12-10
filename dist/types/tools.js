/**
 * Tool Schema Definitions for MCP Tools
 *
 * This module defines input/output schemas and metadata for all MCP tools
 * exposed by the Crystal Structure Generation server.
 */
import { z } from "zod";
/**
 * Common schema components
 */
const Vector3Schema = z.tuple([z.number(), z.number(), z.number()]);
const Matrix3x3Schema = z.tuple([
    z.tuple([z.number(), z.number(), z.number()]),
    z.tuple([z.number(), z.number(), z.number()]),
    z.tuple([z.number(), z.number(), z.number()])
]);
const CrystalSystemSchema = z.enum([
    "triclinic",
    "monoclinic",
    "orthorhombic",
    "tetragonal",
    "trigonal",
    "hexagonal",
    "cubic"
]);
/**
 * Schema for generate_crystal tool
 */
export const GenerateCrystalSchema = z.object({
    composition: z.array(z.string()).min(1)
        .describe("Chemical composition as array of element symbols, e.g., ['Si', 'Si'] or ['Na', 'Cl']"),
    space_group: z.union([
        z.number().int().min(1).max(230),
        z.string()
    ]).describe("Space group number (1-230) or Hermann-Mauguin symbol"),
    num_atoms: z.number().int().positive().optional()
        .describe("Total number of atoms in the unit cell"),
    lattice_params: z.object({
        a: z.number().positive().optional(),
        b: z.number().positive().optional(),
        c: z.number().positive().optional(),
        alpha: z.number().min(0).max(180).optional(),
        beta: z.number().min(0).max(180).optional(),
        gamma: z.number().min(0).max(180).optional()
    }).optional()
        .describe("Optional lattice parameters in Angstroms and degrees"),
    volume_factor: z.number().positive().default(1.0)
        .describe("Relative volume factor (1.0 = standard, >1.0 = larger)"),
    min_distance: z.record(z.string(), z.number().positive()).optional()
        .describe("Minimum distances between element pairs, e.g., {'Si-Si': 2.0, 'Si-O': 1.6}"),
    wyckoff_positions: z.array(z.object({
        element: z.string(),
        wyckoff: z.string(),
        coords: Vector3Schema.optional()
    })).optional()
        .describe("Specify exact Wyckoff positions for atoms"),
    seed: z.number().int().optional()
        .describe("Random seed for reproducibility"),
    max_attempts: z.number().int().positive().default(100)
        .describe("Maximum attempts to generate valid structure")
});
/**
 * Schema for space_group_scan tool
 */
export const SpaceGroupScanSchema = z.object({
    composition: z.array(z.string()).min(1)
        .describe("Chemical composition as array of element symbols"),
    space_groups: z.array(z.number().int().min(1).max(230)).optional()
        .describe("Specific space groups to test, or all 230 if not provided"),
    space_group_range: z.tuple([
        z.number().int().min(1).max(230),
        z.number().int().min(1).max(230)
    ]).optional()
        .describe("Range of space groups to test, e.g., [1, 230]"),
    crystal_systems: z.array(CrystalSystemSchema).optional()
        .describe("Filter by crystal system"),
    num_atoms: z.number().int().positive().optional(),
    volume_factor: z.number().positive().default(1.0),
    parallel: z.boolean().default(false)
        .describe("Enable parallel generation"),
    output_directory: z.string().optional(),
    naming_scheme: z.string().optional()
        .describe("Template for filenames, e.g., '{formula}_{spg}.cif'")
});
/**
 * Schema for make_supercell tool
 */
export const MakeSupercellSchema = z.object({
    structure: z.union([z.string(), z.any()])
        .describe("Crystal structure as CIF file path, JSON, or structure object"),
    scaling_matrix: z.union([
        Matrix3x3Schema,
        Vector3Schema
    ]).describe("3x3 transformation matrix or [nx, ny, nz] scaling factors"),
    wrap_atoms: z.boolean().default(true)
        .describe("Wrap atoms back into unit cell"),
    preserve_symmetry: z.boolean().default(false)
        .describe("Attempt to preserve space group symmetry"),
    min_distance_check: z.boolean().default(true)
        .describe("Check for overlapping atoms after transformation")
});
/**
 * Schema for generate_slab tool
 */
export const GenerateSlabSchema = z.object({
    structure: z.union([z.string(), z.any()])
        .describe("Crystal structure as input"),
    miller_indices: Vector3Schema
        .describe("Miller indices [h, k, l] defining the surface plane"),
    thickness: z.number().int().positive()
        .describe("Number of layers in the slab"),
    vacuum: z.number().positive()
        .describe("Vacuum thickness in Angstroms"),
    center_slab: z.boolean().default(true)
        .describe("Center the slab in the cell"),
    fix_bottom_layers: z.number().int().nonnegative().optional()
        .describe("Number of bottom layers to fix for DFT relaxation"),
    fix_atoms: z.array(z.number().int().nonnegative()).optional()
        .describe("Specific atom indices to fix"),
    orthogonalize: z.boolean().default(false)
        .describe("Create an orthogonal cell"),
    min_slab_size: z.number().positive().optional()
        .describe("Minimum slab size in Angstroms"),
    symmetric: z.boolean().default(true)
        .describe("Make a symmetric slab")
});
/**
 * Schema for analyze_symmetry tool
 */
export const AnalyzeSymmetrySchema = z.object({
    structure: z.union([z.string(), z.any()])
        .describe("Crystal structure to analyze"),
    symprec: z.number().positive().default(1e-3)
        .describe("Symmetry precision tolerance"),
    angle_tolerance: z.number().positive().default(5.0)
        .describe("Angle tolerance in degrees"),
    detect_primitive: z.boolean().default(false)
        .describe("Find and return primitive cell"),
    standardize: z.boolean().default(false)
        .describe("Standardize to conventional cell setting")
});
/**
 * Schema for validate_structure tool
 */
export const ValidateStructureSchema = z.object({
    structure: z.union([z.string(), z.any()])
        .describe("Crystal structure to validate"),
    checks: z.array(z.enum([
        "distances",
        "symmetry",
        "stoichiometry",
        "overlaps",
        "cell_parameters",
        "charge_balance"
    ])).default(["distances", "symmetry", "overlaps"])
        .describe("Validation checks to perform"),
    min_distance: z.number().positive().optional()
        .describe("Minimum allowed interatomic distance"),
    max_coordination: z.number().int().positive().optional()
        .describe("Maximum allowed coordination number"),
    expected_space_group: z.number().int().min(1).max(230).optional()
        .describe("Expected space group number for validation")
});
/**
 * Schema for optimize_structure_mlff tool
 */
export const OptimizeStructureMLFFSchema = z.object({
    structure: z.union([z.string(), z.any()])
        .describe("Crystal structure to optimize"),
    mlff_model: z.enum(["chgnet", "m3gnet", "mace"])
        .describe("Machine learning force field model to use"),
    optimizer: z.enum(["BFGS", "FIRE", "LBFGS", "GPMin"]).default("BFGS")
        .describe("Optimization algorithm"),
    fmax: z.number().positive().default(0.01)
        .describe("Force convergence criterion in eV/Angstrom"),
    steps: z.number().int().positive().default(500)
        .describe("Maximum number of optimization steps"),
    constrain_symmetry: z.boolean().default(false)
        .describe("Preserve space group symmetry during optimization"),
    fix_lattice: z.boolean().default(false)
        .describe("Fix lattice parameters, only optimize atomic positions"),
    fix_volume: z.boolean().default(false)
        .describe("Fix volume but allow cell shape to change"),
    pressure: z.number().default(0.0)
        .describe("External pressure in GPa"),
    trajectory_file: z.string().optional()
        .describe("Path to save optimization trajectory")
});
/**
 * Schema for calculate_energy_mlff tool
 */
export const CalculateEnergyMLFFSchema = z.object({
    structure: z.union([z.string(), z.any()])
        .describe("Crystal structure for energy calculation"),
    mlff_model: z.enum(["chgnet", "m3gnet", "mace"])
        .describe("Machine learning force field model"),
    calculate_forces: z.boolean().default(true)
        .describe("Calculate atomic forces"),
    calculate_stress: z.boolean().default(true)
        .describe("Calculate stress tensor")
});
/**
 * Schema for ground_state_search tool
 */
export const GroundStateSearchSchema = z.object({
    composition: z.array(z.string()).min(1)
        .describe("Chemical composition"),
    space_groups: z.array(z.number().int().min(1).max(230)).optional()
        .describe("Specific space groups to search, or all 230 if not provided"),
    num_structures_per_group: z.number().int().positive().default(3)
        .describe("Number of random structures to generate per space group"),
    mlff_model: z.enum(["chgnet", "m3gnet", "mace"])
        .describe("MLFF model for energy calculations"),
    optimization_settings: z.object({
        optimizer: z.enum(["BFGS", "FIRE", "LBFGS"]).default("BFGS"),
        fmax: z.number().positive().default(0.01),
        steps: z.number().int().positive().default(500),
        constrain_symmetry: z.boolean().default(true)
    }),
    parallel: z.boolean().default(false)
        .describe("Enable parallel computation"),
    n_workers: z.number().int().positive().optional()
        .describe("Number of parallel workers"),
    output_directory: z.string().optional()
        .describe("Directory to save results"),
    save_trajectories: z.boolean().default(false)
        .describe("Save optimization trajectories")
});
/**
 * Schema for export_structure tool
 */
export const ExportStructureSchema = z.object({
    structure: z.union([z.string(), z.any()])
        .describe("Crystal structure to export"),
    formats: z.array(z.enum([
        "cif",
        "poscar",
        "contcar",
        "xyz",
        "extxyz",
        "json",
        "pdb",
        "pwscf",
        "castep",
        "lammps"
    ])).min(1)
        .describe("Output file formats"),
    include_metadata: z.boolean().default(true)
        .describe("Include metadata in output files"),
    dft_software: z.enum(["vasp", "qe", "cp2k", "castep"]).optional()
        .describe("Target DFT software for format-specific options"),
    output_directory: z.string().optional()
        .describe("Directory to write output files")
});
/**
 * All tool definitions
 */
export const TOOL_DEFINITIONS = [
    {
        name: "generate_crystal",
        description: "Generate a crystal structure with specified composition and space group using PyXtal",
        inputSchema: GenerateCrystalSchema,
        annotations: {
            readOnlyHint: false,
            destructiveHint: false,
            idempotentHint: true,
            openWorldHint: true
        }
    },
    {
        name: "generate_space_group_scan",
        description: "Generate structures across multiple space groups for the same composition",
        inputSchema: SpaceGroupScanSchema,
        annotations: {
            readOnlyHint: false,
            destructiveHint: false,
            idempotentHint: true,
            openWorldHint: true
        }
    },
    {
        name: "make_supercell",
        description: "Create a supercell from an existing structure using scaling matrix",
        inputSchema: MakeSupercellSchema,
        annotations: {
            readOnlyHint: false,
            destructiveHint: false,
            idempotentHint: true,
            openWorldHint: false
        }
    },
    {
        name: "generate_slab",
        description: "Generate a surface slab from a bulk structure with specified Miller indices",
        inputSchema: GenerateSlabSchema,
        annotations: {
            readOnlyHint: false,
            destructiveHint: false,
            idempotentHint: true,
            openWorldHint: false
        }
    },
    {
        name: "analyze_symmetry",
        description: "Analyze and detect symmetry properties of a crystal structure",
        inputSchema: AnalyzeSymmetrySchema,
        annotations: {
            readOnlyHint: true,
            destructiveHint: false,
            idempotentHint: true,
            openWorldHint: false
        }
    },
    {
        name: "validate_structure",
        description: "Validate a crystal structure for physical and chemical correctness",
        inputSchema: ValidateStructureSchema,
        annotations: {
            readOnlyHint: true,
            destructiveHint: false,
            idempotentHint: true,
            openWorldHint: false
        }
    },
    {
        name: "optimize_structure_mlff",
        description: "Optimize crystal structure using machine learning force fields (CHGNet, M3GNet, MACE)",
        inputSchema: OptimizeStructureMLFFSchema,
        annotations: {
            readOnlyHint: false,
            destructiveHint: false,
            idempotentHint: false,
            openWorldHint: false
        }
    },
    {
        name: "calculate_energy_mlff",
        description: "Calculate energy, forces, and stress using MLFF without optimization",
        inputSchema: CalculateEnergyMLFFSchema,
        annotations: {
            readOnlyHint: true,
            destructiveHint: false,
            idempotentHint: true,
            openWorldHint: false
        }
    },
    {
        name: "ground_state_search",
        description: "Search for ground state structure across multiple space groups using MLFF",
        inputSchema: GroundStateSearchSchema,
        annotations: {
            readOnlyHint: false,
            destructiveHint: false,
            idempotentHint: false,
            openWorldHint: false
        }
    },
    {
        name: "export_structure",
        description: "Export crystal structure to multiple file formats (CIF, POSCAR, XYZ, etc.)",
        inputSchema: ExportStructureSchema,
        annotations: {
            readOnlyHint: true,
            destructiveHint: false,
            idempotentHint: true,
            openWorldHint: false
        }
    }
];
//# sourceMappingURL=tools.js.map