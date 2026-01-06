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
const IntVector3Schema = z.tuple([z.number().int(), z.number().int(), z.number().int()]);
const Matrix3x3Schema = z.tuple([Vector3Schema, Vector3Schema, Vector3Schema]);
const ElementSymbolSchema = z.string().regex(/^[A-Z][a-z]?$/, "Invalid element symbol");
const ElementPairKeySchema = z.string().regex(/^[A-Z][a-z]?-[A-Z][a-z]?$/, "Expected element pair like 'Si-O'");
const StructureInputSchema = z.union([z.string(), z.any()])
    .describe('Structure input: either a prototype object like {"prototype": "diamond", "elements": {"A": "Si"}, "lattice_constant": 5.43} or {"prototype": "fcc", "elements": {"A": "Al"}, "lattice_constant": 4.05} or {"prototype": "rocksalt", "elements": {"A": "Na", "B": "Cl"}, "lattice_constant": 5.64}');
const SpaceGroupSchema = z.union([
    z.number().int().min(1).max(230),
    z.string().min(1)
]);
const LatticeParamsSchema = z.object({
    a: z.number().positive(),
    b: z.number().positive(),
    c: z.number().positive(),
    alpha: z.number().min(0).max(180),
    beta: z.number().min(0).max(180),
    gamma: z.number().min(0).max(180)
});
const CompositionSchema = z.array(ElementSymbolSchema).min(1);
const StrainTensorSchema = z.union([
    Matrix3x3Schema,
    z.array(z.number()).length(9)
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
    composition: CompositionSchema
        .describe('Chemical composition as array: ["Si", "Si"] for Si2, ["Ga", "As"] for GaAs, ["Ba", "Ti", "O", "O", "O"] for BaTiO3'),
    space_group: SpaceGroupSchema
        .describe("Space group number: 227 for diamond-Si, 216 for zincblende-GaAs, 225 for rocksalt-NaCl, 221 for perovskite"),
    num_atoms: z.number().int().positive().optional()
        .describe("Total number of atoms in the unit cell"),
    lattice_params: LatticeParamsSchema.optional()
        .describe("Optional lattice parameters in Angstroms and degrees"),
    volume_factor: z.number().positive().default(1.0)
        .describe("Relative volume factor (1.0 = standard, >1.0 = larger)"),
    min_distance: z.record(ElementPairKeySchema, z.number().positive()).optional()
        .describe("Minimum distances between element pairs, e.g., {'Si-Si': 2.0, 'Si-O': 1.6}"),
    wyckoff_positions: z.array(z.object({
        element: ElementSymbolSchema,
        wyckoff: z.string().min(1),
        coords: Vector3Schema.optional()
    })).optional()
        .describe("Specify exact Wyckoff positions for atoms"),
    seed: z.number().int().optional()
        .describe("Random seed for reproducibility"),
    dimensionality: z.union([
        z.literal(0),
        z.literal(1),
        z.literal(2),
        z.literal(3)
    ]).default(3)
        .describe("Crystal dimensionality: 0=cluster, 1=rod, 2=slab/layer, 3=bulk"),
    max_attempts: z.number().int().positive().default(100)
        .describe("Maximum attempts to generate valid structure"),
    output_directory: z.string().optional()
        .describe("Directory to save structure files (CIF, POSCAR, XYZ, JSON). If not provided, files are returned in response but not saved to disk.")
});
/**
 * Schema for comprehensive_generate tool
 *
 * Unified entry point for all 51+ generator operations across 18 categories:
 * bulk, two_d, surface, molecule, twist, defect, electronic, thermoelectric,
 * battery, catalyst, adsorption, magnetic, nanotube, quantum, photonic,
 * quality_control, high_pressure, external_fields
 */
export const ComprehensiveGenerateSchema = z.object({
    operation: z.string()
        .describe("Generator operation name. Use 'list_all' to see all available operations, 'list_category' with category param, or 'operation_info' with operation_name param."),
    category: z.string().optional()
        .describe("Category name when using operation='list_category'. Categories: bulk, two_d, surface, molecule, twist, defect, electronic, thermoelectric, battery, catalyst, adsorption, magnetic, nanotube, quantum, photonic, quality_control, high_pressure, external_fields"),
    operation_name: z.string().optional()
        .describe("Operation name to get info for when using operation='operation_info'"),
    // Common parameters - passed through to the generator function
    spacegroup: z.number().int().min(1).max(230).optional()
        .describe("Space group number (1-230) for bulk generation"),
    elements: z.union([
        z.array(z.string()),
        z.record(z.string(), z.string())
    ]).optional()
        .describe("Elements for structure generation"),
    composition: z.array(z.number().int()).optional()
        .describe("Composition (number of each element)"),
    material: z.string().optional()
        .describe("Material name/type for specific generators"),
    supercell: z.array(z.number().int()).optional()
        .describe("Supercell dimensions [a, b, c]"),
    prototype: z.string().optional()
        .describe("Prototype name: rocksalt, perovskite, zincblende, wurtzite, rutile, spinel, bcc, fcc, hcp, diamond"),
    framework: z.string().optional()
        .describe("Zeolite framework type: MFI, FAU, LTA, CHA, BEA, MOR"),
    clathrate: z.string().optional()
        .describe("Clathrate type: Ba8Ga16Ge30, Na8Si46, etc."),
    materials: z.array(z.string()).optional()
        .describe("List of materials for heterostructures"),
    twist_angle: z.number().optional()
        .describe("Twist angle in degrees for bilayer structures"),
    n_atoms: z.number().int().positive().optional()
        .describe("Number of atoms for various generators"),
    size_nm: z.number().positive().optional()
        .describe("Size in nanometers for nanostructures"),
    thickness_QL: z.number().int().positive().optional()
        .describe("Thickness in quintuple layers for topological materials"),
    pressure_GPa: z.number().positive().optional()
        .describe("Pressure in GPa for high-pressure phases"),
    // Additional parameters are passed through dynamically
    list_available: z.boolean().optional()
        .describe("Include list of available options in response")
}).passthrough(); // Allow additional operation-specific parameters to pass through
/**
 * Schema for space_group_scan tool
 */
export const SpaceGroupScanSchema = z.object({
    composition: CompositionSchema
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
    structure: StructureInputSchema
        .describe("Crystal structure as CIF file path, JSON, or structure object"),
    scaling_matrix: z.union([
        z.string().describe("Preset matrix name (e.g., 'sqrt3', 'root2', '2x2x2')"),
        IntVector3Schema.describe("Scaling vector [nx, ny, nz]"),
        Matrix3x3Schema.describe("3x3 scaling matrix")
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
    structure: StructureInputSchema
        .describe('Bulk structure as prototype object. Examples: {"prototype": "diamond", "elements": {"A": "Si"}, "lattice_constant": 5.43} for Si, {"prototype": "fcc", "elements": {"A": "Al"}, "lattice_constant": 4.05} for Al, {"prototype": "rocksalt", "elements": {"A": "Na", "B": "Cl"}, "lattice_constant": 5.64} for NaCl'),
    miller_indices: Vector3Schema
        .describe("Miller indices [h, k, l] as array, e.g. [1,0,0] for (100), [1,1,0] for (110), [1,1,1] for (111)"),
    thickness: z.number().int().positive()
        .describe("Number of atomic layers in the slab (integer, e.g. 3, 4, 5)"),
    vacuum: z.number().positive()
        .describe("Vacuum thickness in Angstroms"),
    center_slab: z.boolean().default(true)
        .describe("Center the slab in the cell"),
    fix_bottom_layers: z.number().int().nonnegative().optional()
        .describe("Number of bottom layers to fix for DFT relaxation"),
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
    structure: StructureInputSchema
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
    structure: StructureInputSchema
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
    structure: StructureInputSchema
        .describe("Crystal structure to optimize"),
    mlff_model: z.enum(["chgnet", "m3gnet", "mace"])
        .describe("Machine learning force field model to use"),
    optimizer: z.enum(["BFGS", "FIRE", "LBFGS"]).default("BFGS")
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
    structure: StructureInputSchema
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
    composition: CompositionSchema
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
    structure: StructureInputSchema
        .describe("Crystal structure to export"),
    formats: z.array(z.enum([
        "cif",
        "poscar",
        "xyz",
        "json"
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
/**
 * Schema for generate_visualization tool
 */
export const VisualizationSchema = z.object({
    structure: StructureInputSchema.describe("Structure object returned by generate_crystal"),
    format: z.enum(["html", "png"]).optional().default("html").describe("Output format: 'html' for interactive 3D, 'png' for static image"),
    output_file: z.string().optional().describe("Path to save the visualization file")
});
/**
 * Schema for create_defect tool
 */
export const CreateDefectSchema = z.object({
    structure: StructureInputSchema
        .describe("Crystal structure to modify"),
    defect_type: z.enum(["vacancy", "substitution", "interstitial"])
        .describe("Type of defect to create"),
    defect_site: z.number().int().nonnegative()
        .describe("Index of the atom site to modify (0-indexed)"),
    defect_species: ElementSymbolSchema.optional()
        .describe("Element symbol for substitution or interstitial (required for these types)"),
    concentration: z.number().positive().default(1.0)
        .describe("Defect concentration (currently unused)")
});
/**
 * Schema for generate_molecular_crystal tool
 */
export const GenerateMolecularCrystalSchema = z.object({
    molecules: z.array(z.string()).min(1)
        .describe("List of molecules (formulas or names), e.g., ['H2O']"),
    space_group: SpaceGroupSchema
        .describe("Space group number (1-230) or Hermann-Mauguin symbol"),
    num_molecules: z.number().int().positive().optional()
        .describe("Total number of molecules in the unit cell"),
    lattice_params: LatticeParamsSchema.optional()
        .describe("Optional lattice parameters"),
    volume_factor: z.number().positive().default(1.0)
        .describe("Relative volume factor"),
    min_distance: z.record(ElementPairKeySchema, z.number().positive()).optional()
        .describe("Minimum distances"),
    seed: z.number().int().optional()
        .describe("Random seed")
});
/**
 * Schema for generate_nanostructure tool
 */
const NanotubeParamsSchema = z.object({
    n: z.number().int().positive().optional().describe("Chiral index n"),
    m: z.number().int().nonnegative().optional().describe("Chiral index m"),
    length: z.number().int().positive().optional().describe("Number of unit cells along tube"),
    bond: z.number().positive().optional().describe("Bond length (default 1.42)"),
    vacuum: z.number().positive().optional().describe("Vacuum padding in Angstroms")
});
const GrapheneParamsSchema = z.object({
    formula: z.string().optional().describe("Element symbol or formula (e.g. 'C')"),
    size: IntVector3Schema.optional().describe("Supercell size [x, y, z]"),
    a: z.number().positive().optional().describe("Lattice constant a in Angstroms"),
    c: z.number().positive().optional().describe("Out-of-plane separation in Angstroms"),
    vacuum: z.number().positive().optional().describe("Vacuum padding in Angstroms")
});
const NanoribbonParamsSchema = z.object({
    width: z.number().int().positive().optional(),
    length: z.number().int().positive().optional(),
    type: z.enum(["armchair", "zigzag"]).optional(),
    ribbon_type: z.enum(["armchair", "zigzag"]).optional(),
    saturated: z.boolean().optional(),
    vacuum: z.number().positive().optional()
}).refine(data => !(data.type && data.ribbon_type), {
    message: "Use either type or ribbon_type, not both"
});
const FullereneParamsSchema = z.object({
    name: z.string().optional().describe("Molecule name (e.g. 'C60')")
});
const Mos2ParamsSchema = z.object({
    formula: z.string().optional().describe("Chemical formula (e.g. 'MoS2')"),
    kind: z.string().optional().describe("Polytype (e.g. '2H')"),
    a: z.number().positive().optional().describe("Lattice constant a in Angstroms"),
    thickness: z.number().positive().optional().describe("Layer thickness in Angstroms"),
    vacuum: z.number().positive().optional().describe("Vacuum padding in Angstroms"),
    size: IntVector3Schema.optional().describe("Supercell size [x, y, z]")
});
const NanowireParamsSchema = z.object({
    formula: z.string().optional().describe("Bulk formula to carve from"),
    radius: z.number().positive().optional().describe("Wire radius in Angstroms"),
    length: z.number().int().positive().optional().describe("Number of unit cells along axis"),
    vacuum: z.number().positive().optional().describe("Vacuum padding in Angstroms")
});
export const GenerateNanostructureSchema = z.discriminatedUnion("type", [
    z.object({
        type: z.literal("nanotube").describe("Carbon nanotube with chiral indices (n,m)"),
        params: NanotubeParamsSchema.describe("Nanotube params: {n: 5, m: 5, length: 3, vacuum: 15}")
    }),
    z.object({
        type: z.literal("graphene").describe("Graphene monolayer sheet"),
        params: GrapheneParamsSchema.describe("Graphene params: {size: [4, 4, 1], vacuum: 15.0}")
    }),
    z.object({
        type: z.literal("nanoribbon").describe("Graphene nanoribbon"),
        params: NanoribbonParamsSchema
    }),
    z.object({
        type: z.literal("fullerene").describe("Fullerene cage"),
        params: FullereneParamsSchema
    }),
    z.object({
        type: z.literal("mos2").describe("MoS2 2D sheet"),
        params: Mos2ParamsSchema
    }),
    z.object({
        type: z.literal("nanowire").describe("Nanowire carved from bulk"),
        params: NanowireParamsSchema
    })
]);
/**
 * Schema for explore_symmetry_relations tool
 */
export const ExploreSymmetryRelationsSchema = z.object({
    operation: z.enum(["get_subgroups", "get_path"]).describe("Operation to perform: find subgroups or find path between groups"),
    space_group: z.number().int().min(1).max(230).optional().describe("Space group number for 'get_subgroups' operation"),
    start_group: z.number().int().min(1).max(230).optional().describe("Starting space group for 'get_path' operation"),
    end_group: z.number().int().min(1).max(230).optional().describe("Target space group for 'get_path' operation"),
    max_depth: z.number().int().min(1).max(10).default(5).optional().describe("Maximum search depth for path finding")
}).refine(data => {
    if (data.operation === "get_subgroups" && !data.space_group)
        return false;
    if (data.operation === "get_path" && (!data.start_group || !data.end_group))
        return false;
    return true;
}, {
    message: "Missing required parameters for selected operation"
});
/**
 * Schema for build_molecule tool
 *
 * Universal molecule generation supporting:
 * - Common names: H2O, CO2, aspirin, caffeine, PTCDA, benzene
 * - SMILES strings: c1ccccc1 (benzene), CC(=O)OC1=CC=CC=C1C(=O)O (aspirin)
 * - IUPAC names: perylene-3,4,9,10-tetracarboxylic dianhydride
 * - PubChem CIDs: 2244 for aspirin
 *
 * Priority: Local DB → Aliases → RDKit from SMILES → PubChem API → OPSIN
 */
export const BuildMoleculeSchema = z.object({
    name: z.string()
        .describe("Molecule identifier - pass the identifier directly as this value. Accepts: " +
        "common names (aspirin, caffeine, PTCDA, benzene), " +
        "SMILES strings (pass the SMILES directly like 'c1ccccc1' for benzene, 'CCO' for ethanol), " +
        "IUPAC names (perylene-3,4,9,10-tetracarboxylic dianhydride), " +
        "or PubChem CIDs (2244 for aspirin). " +
        "For SMILES: set name='c1ccc2ccccc2c1' (the SMILES itself, not a description)."),
    input_type: z.enum(["auto", "name", "smiles", "iupac", "cid"]).default("auto").optional()
        .describe("Input type hint: 'auto' (default, auto-detect), 'name' (common name), 'smiles' (SMILES string), 'iupac' (IUPAC systematic name), 'cid' (PubChem CID)"),
    optimize: z.boolean().default(true).optional()
        .describe("Optimize 3D geometry using MMFF94/UFF force field (default: true)"),
    vacuum: z.number().default(10.0).optional()
        .describe("Vacuum padding around the molecule in Angstroms (default: 10.0)")
});
/**
 * Schema for build_molecular_cluster tool
 *
 * Generate molecular clusters for quantum chemistry:
 * - Homo/hetero dimers, trimers, n-mers
 * - π-π stacking (parallel, antiparallel, offset)
 * - T-shaped (edge-to-face) arrangements
 * - H-bonded clusters
 * - Custom arrangements with full rotation control
 */
const MoleculeSpecSchema = z.object({
    identifier: z.string()
        .describe("Molecule identifier (name, SMILES, IUPAC, CID, ChEMBL ID)"),
    count: z.number().int().min(1).default(1).optional()
        .describe("Number of copies of this molecule (default: 1)"),
    input_type: z.enum(["auto", "name", "smiles", "iupac", "cid"]).default("auto").optional()
        .describe("Input type hint")
});
const Position3DSchema = z.object({
    x: z.number().default(0),
    y: z.number().default(0),
    z: z.number().default(0)
}).describe("3D position in Angstroms");
const Rotation3DSchema = z.object({
    x: z.number().default(0).describe("Rotation around x-axis in degrees"),
    y: z.number().default(0).describe("Rotation around y-axis in degrees"),
    z: z.number().default(0).describe("Rotation around z-axis in degrees")
}).describe("3D rotation in degrees");
export const BuildMolecularClusterSchema = z.object({
    // Molecules list
    molecules: z.array(MoleculeSpecSchema).min(1)
        .describe("List of molecules to combine. Example: [{\"identifier\": \"benzene\", \"count\": 2}] for benzene dimer, or [{\"identifier\": \"water\"}, {\"identifier\": \"benzene\"}] for hetero-dimer"),
    // Stacking/arrangement type  
    stacking: z.enum([
        "auto", // Auto-detect based on molecule types
        "pi_pi_parallel", // Face-to-face π-stacking (3.4Å)
        "parallel", // Alias for pi_pi_parallel
        "stacked", // Alias for pi_pi_parallel
        "pi_pi_antiparallel", // π-stacking with 180° rotation
        "antiparallel", // Alias for pi_pi_antiparallel
        "pi_pi_offset", // Offset/slip-stacked
        "offset", // Alias for pi_pi_offset
        "slip_stacked", // Alias for pi_pi_offset
        "t_shaped", // Edge-to-face perpendicular
        "edge_to_face", // Alias for t_shaped
        "herringbone", // Alternating tilted (organic crystals)
        "h_bonded", // Hydrogen bonded (2.8Å)
        "hydrogen_bonded", // Alias for h_bonded
        "van_der_waals", // General vdW contact
        "vdw", // Alias for van_der_waals
        "linear", // In a line along axis
        "circular", // Ring arrangement
        "ring", // Alias for circular
        "spherical", // 3D spherical distribution
        "swastika", // 4-molecule cross pattern
        "swastic", // Alias for swastika (typo tolerance)
        "custom" // User-defined positions/rotations
    ]).default("auto").optional()
        .describe("Arrangement type. IMPORTANT: Use 'linear' if user specifies a direction (e.g. 'along x') or distance. Use 'auto' only for general stacking."),
    // Distance control
    intermolecular_distance: z.number().positive().optional()
        .describe("Distance between molecule centers in Angstroms. Defaults: π-stacking=3.4Å, H-bonded=2.8Å, vdW=3.5Å"),
    // Offsets for stacking
    offset_x: z.number().default(0).optional()
        .describe("Lateral offset in x-direction (Angstroms)"),
    offset_y: z.number().default(0).optional()
        .describe("Lateral offset in y-direction (Angstroms)"),
    // Global rotation (applied to entire cluster)
    rotation_x: z.number().default(0).optional()
        .describe("Rotate entire cluster around x-axis (degrees)"),
    rotation_y: z.number().default(0).optional()
        .describe("Rotate entire cluster around y-axis (degrees)"),
    rotation_z: z.number().default(0).optional()
        .describe("Rotate entire cluster around z-axis (degrees)"),
    // Per-molecule rotation increment (for stacking)
    rotation_per_molecule: z.number().default(0).optional()
        .describe("Incremental rotation per stacked molecule (degrees). E.g., 45° for spiral stacking."),
    // Axis for linear arrangement
    axis: z.enum(["x", "y", "z"]).default("z").optional()
        .describe("Axis for linear arrangement (default: z)"),
    // Custom positions (for 'custom' stacking)
    positions: z.array(Position3DSchema).optional()
        .describe("Custom positions for each molecule when stacking='custom'"),
    // Custom rotations (for 'custom' stacking)
    rotations: z.array(Rotation3DSchema).optional()
        .describe("Custom rotations for each molecule when stacking='custom'"),
    // Optimization
    optimize: z.boolean().default(false).optional()
        .describe("Optimize cluster geometry with force field (can be slow for large clusters)"),
    // Vacuum box
    vacuum: z.number().default(10.0).optional()
        .describe("Vacuum padding around the cluster in Angstroms")
});
/**
 * Schema for create_alloy tool
 */
export const CreateAlloySchema = z.object({
    structure: StructureInputSchema
        .describe("Base crystal structure"),
    substitutions: z.record(ElementSymbolSchema, z.object({
        element: ElementSymbolSchema.describe("Element to substitute with"),
        concentration: z.number().min(0).max(1).describe("Concentration of substitution (0-1)")
    })).describe("Substitution rules, e.g. {'Si': {element: 'Ge', concentration: 0.5}}"),
    seed: z.number().int().optional()
        .describe("Random seed for reproducibility")
});
/**
 * Schema for create_heterostructure tool
 */
export const CreateHeterostructureSchema = z.object({
    substrate: StructureInputSchema
        .describe("Substrate structure"),
    overlayer: StructureInputSchema
        .describe("Overlayer structure to stack on top"),
    interface_distance: z.number().positive().default(3.0)
        .describe("Distance between substrate and overlayer (Angstroms)"),
    vacuum: z.number().positive().default(10.0)
        .describe("Vacuum padding (Angstroms)")
});
/**
 * Schema for add_adsorbate tool
 */
export const AddAdsorbateSchema = z.object({
    structure: StructureInputSchema
        .describe("Surface structure"),
    molecule: StructureInputSchema
        .describe("Molecule to adsorb (name or structure object)"),
    site_index: z.number().int().nonnegative()
        .describe("Index of the surface atom to adsorb on"),
    distance: z.number().positive().default(2.0)
        .describe("Height above the surface site (Angstroms)")
});
/**
 * Schema for apply_strain tool
 */
export const ApplyStrainSchema = z.object({
    structure: StructureInputSchema
        .describe("Crystal structure to strain"),
    strain_tensor: StrainTensorSchema.optional()
        .describe("3x3 strain tensor matrix or flattened array [e11, e12, e13, e21...]"),
    strain_type: z.enum(["biaxial", "uniaxial", "hydrostatic"]).optional(),
    strain_value: z.number().optional().describe("Strain percentage (e.g. 0.05 for 5%)")
}).refine(data => {
    if (data.strain_tensor) {
        return true;
    }
    return data.strain_type !== undefined && data.strain_value !== undefined;
}, {
    message: "Must provide either strain_tensor or both strain_type and strain_value"
});
// ============================================================================
// ADVANCED STRUCTURE SCHEMAS
// ============================================================================
/**
 * Generate prototype structure (rocksalt, perovskite, etc.)
 */
export const GeneratePrototypeSchema = z.object({
    prototype: z.enum([
        "rocksalt", "zincblende", "wurtzite", "fluorite", "antifluorite",
        "perovskite", "spinel", "heusler", "rutile", "diamond", "bcc", "fcc", "hcp"
    ]).describe("Prototype type: diamond for Si/Ge/C, fcc for Al/Cu/Au, bcc for Fe/W, rocksalt for NaCl/MgO, zincblende for GaAs/ZnS, perovskite for BaTiO3/CaTiO3"),
    elements: z.record(z.string(), ElementSymbolSchema)
        .describe('Element mapping: {"A": "Si"} for diamond/fcc/bcc, {"A": "Na", "B": "Cl"} for rocksalt/zincblende, {"A": "Ba", "B": "Ti", "X": "O"} for perovskite'),
    lattice_constant: z.number().positive().optional()
        .describe("Lattice constant in Angstroms: Si=5.43, Ge=5.66, Al=4.05, Cu=3.61, Fe=2.87, NaCl=5.64, GaAs=5.65"),
    c_over_a: z.number().positive().default(1.0)
        .describe("c/a ratio for non-cubic systems")
});
/**
 * Generate twisted bilayer structure
 */
export const GenerateTwistedBilayerSchema = z.object({
    material: z.enum(["graphene", "MoS2", "WS2", "hBN"]).default("graphene")
        .describe("Base material: graphene, MoS2, WS2, or hBN"),
    twist_angle: z.number().min(0).max(60)
        .describe("Twist angle in degrees (e.g. 1.1 for magic angle, 21.8, 30)"),
    layers: z.number().int().min(2).default(2)
        .describe("Number of layers"),
    stacking: z.enum(["AA", "AB"]).default("AB")
        .describe("Initial stacking"),
    interlayer_distance: z.number().positive().default(3.35)
        .describe("Interlayer spacing in Angstroms"),
    vacuum: z.number().positive().default(15.0)
        .describe("Vacuum padding in Angstroms")
});
/**
 * Generate high-entropy alloy structure
 */
export const GenerateHighEntropyAlloySchema = z.object({
    elements: z.array(ElementSymbolSchema).min(4)
        .describe("List of 4+ elements for HEA"),
    concentrations: z.array(z.number()).optional()
        .describe("Concentrations (default: equimolar)"),
    structure_type: z.enum(["fcc", "bcc", "hcp"]).default("fcc")
        .describe("Base crystal structure"),
    supercell: IntVector3Schema.default([3, 3, 3])
        .describe("Supercell size"),
    lattice_constant: z.number().positive().optional()
        .describe("Lattice constant (estimated if not provided)"),
    seed: z.number().int().optional()
        .describe("Random seed for reproducibility")
});
/**
 * Generate 2D material structure
 */
export const Generate2DMaterialSchema = z.object({
    material: z.enum(["hBN", "MoS2", "WS2", "MoSe2", "WSe2", "phosphorene", "silicene", "MXene"])
        .describe("2D material: hBN, MoS2, WS2, MoSe2, WSe2, phosphorene, silicene, or MXene"),
    size: IntVector3Schema.default([1, 1, 1])
        .describe("Supercell size as [nx, ny, nz] array, e.g. [3, 3, 1] for 3x3 supercell"),
    vacuum: z.number().positive().default(15.0)
        .describe("Vacuum padding"),
    extra_params: z.record(z.string(), z.any()).optional()
        .describe("Material-specific parameters")
});
/**
 * Generate MOF structure
 */
export const GenerateMOFSchema = z.object({
    mof_type: z.enum(["MOF-5", "HKUST-1", "UiO-66", "ZIF-8"])
        .describe("MOF type"),
    functionalization: z.string().optional()
        .describe("Linker functionalization"),
    size: IntVector3Schema.default([1, 1, 1])
        .describe("Supercell size")
});
/**
 * Generate cage structure (fullerenes, clathrates)
 */
export const GenerateCageSchema = z.object({
    cage_type: z.enum(["C60", "C70", "C80", "clathrate_I", "clathrate_II"])
        .describe("Cage structure type"),
    guest: z.string().optional()
        .describe("Guest atom for endohedral structures"),
    vacuum: z.number().positive().default(5.0)
        .describe("Vacuum padding in Angstroms")
});
/**
 * All tool definitions
 */
export const TOOL_DEFINITIONS = [
    // {
    //   name: "comprehensive_generate",
    //   description: "Unified generator for all crystal structures. Access 50+ operations across 18 categories: bulk, 2D, surface, molecule, twist, defect, electronic, etc.",
    //   inputSchema: ComprehensiveGenerateSchema,
    //   annotations: {
    //     readOnlyHint: false,
    //     destructiveHint: false,
    //     idempotentHint: false,
    //     openWorldHint: true
    //   }
    // },
    {
        name: "build_molecule",
        description: "Generate 3D molecular structure from ANY identifier. MANDATORY: You MUST use this tool whenever the user asks to generate, show, create, or visualize a molecule. Do NOT provide a text-only description without this tool. " +
            "Accepts: common names (H2O, CO2, aspirin, caffeine, PTCDA, benzene, pentacene), " +
            "SMILES strings (c1ccccc1, CCO, CC(=O)O), " +
            "IUPAC names (perylene-3,4,9,10-tetracarboxylic dianhydride), " +
            "or PubChem CIDs. " +
            "Uses RDKit for 3D generation and PubChem API for name resolution. " +
            "Supports ~130M molecules via PubChem lookup.",
        inputSchema: BuildMoleculeSchema,
        annotations: {
            readOnlyHint: false,
            destructiveHint: false,
            idempotentHint: true,
            openWorldHint: true
        }
    },
    {
        name: "build_molecular_cluster",
        description: "Generate molecular clusters (dimers, trimers, n-mers) or COMBINATIONS of different molecules. " +
            "MANDATORY: Use this tool whenever the user asks for: " +
            "1. A cluster/dimer/stack (e.g. 'benzene dimer') " +
            "2. A LIST of molecules (e.g. 'generate water and benzene', 'create NTCDA and PTCDA'). " +
            "Combines any molecules from build_molecule with various arrangements. " +
            "Auto-detects optimal stacking. Supports full rotation. " +
            "Examples: " +
            "- 'benzene dimer': {molecules: [{identifier: 'benzene', count: 2}], stacking: 'pi_pi_parallel'} " +
            "- 'water and benzene': {molecules: [{identifier: 'water'}, {identifier: 'benzene'}], stacking: 'auto'}",
        inputSchema: BuildMolecularClusterSchema,
        annotations: {
            readOnlyHint: false,
            destructiveHint: false,
            idempotentHint: true,
            openWorldHint: true
        }
    }
    // },
    // {
    //   name: "explore_symmetry_relations",
    //   description: "Explore group-subgroup relationships and find symmetry pathways between space groups.",
    //   inputSchema: ExploreSymmetryRelationsSchema,
    //   annotations: {
    //     readOnlyHint: true,
    //     destructiveHint: false,
    //     idempotentHint: true,
    //     openWorldHint: false
    //   }
    // },
    // {
    //   name: "generate_crystal",
    //   description: "Generate a crystal structure with specified composition and space group using PyXtal",
    //   inputSchema: GenerateCrystalSchema,
    //   annotations: {
    //     readOnlyHint: false,
    //     destructiveHint: false,
    //     idempotentHint: true,
    //     openWorldHint: true
    //   }
    // },
    // {
    //   name: "generate_molecular_crystal",
    //   description: "Generate a molecular crystal structure using PyXtal",
    //   inputSchema: GenerateMolecularCrystalSchema,
    //   annotations: {
    //     readOnlyHint: false,
    //     destructiveHint: false,
    //     idempotentHint: true,
    //     openWorldHint: true
    //   }
    // },
    // {
    //   name: "generate_nanostructure",
    //   description: "Generate nanostructures (nanotubes, graphene, ribbons, fullerenes)",
    //   inputSchema: GenerateNanostructureSchema,
    //   annotations: {
    //     readOnlyHint: false,
    //     destructiveHint: false,
    //     idempotentHint: true,
    //     openWorldHint: true
    //   }
    // },
    // {
    //   name: "generate_space_group_scan",
    //   description: "Generate structures across multiple space groups for the same composition",
    //   inputSchema: SpaceGroupScanSchema,
    //   annotations: {
    //     readOnlyHint: false,
    //     destructiveHint: false,
    //     idempotentHint: true,
    //     openWorldHint: true
    //   }
    // },
    // {
    //   name: "make_supercell",
    //   description: "Create a supercell from an existing structure using scaling matrix",
    //   inputSchema: MakeSupercellSchema,
    //   annotations: {
    //     readOnlyHint: false,
    //     destructiveHint: false,
    //     idempotentHint: true,
    //     openWorldHint: false
    //   }
    // },
    // {
    //   name: "generate_slab",
    //   description: "Generate a surface slab from a bulk structure with specified Miller indices",
    //   inputSchema: GenerateSlabSchema,
    //   annotations: {
    //     readOnlyHint: false,
    //     destructiveHint: false,
    //     idempotentHint: true,
    //     openWorldHint: false
    //   }
    // },
    // {
    //   name: "create_defect",
    //   description: "Create point defects (vacancy, substitution, interstitial) in a crystal structure",
    //   inputSchema: CreateDefectSchema,
    //   annotations: {
    //     readOnlyHint: false,
    //     destructiveHint: false,
    //     idempotentHint: true,
    //     openWorldHint: false
    //   }
    // },
    // {
    //   name: "create_alloy",
    //   description: "Create a substitutional alloy with random mixing",
    //   inputSchema: CreateAlloySchema,
    //   annotations: {
    //     readOnlyHint: false,
    //     destructiveHint: false,
    //     idempotentHint: true,
    //     openWorldHint: false
    //   }
    // },
    // {
    //   name: "create_heterostructure",
    //   description: "Create a vertical heterostructure by stacking two layers",
    //   inputSchema: CreateHeterostructureSchema,
    //   annotations: {
    //     readOnlyHint: false,
    //     destructiveHint: false,
    //     idempotentHint: true,
    //     openWorldHint: false
    //   }
    // },
    // {
    //   name: "add_adsorbate",
    //   description: "Add an adsorbate molecule to a surface structure",
    //   inputSchema: AddAdsorbateSchema,
    //   annotations: {
    //     readOnlyHint: false,
    //     destructiveHint: false,
    //     idempotentHint: true,
    //     openWorldHint: false
    //   }
    // },
    // {
    //   name: "apply_strain",
    //   description: "Apply strain tensor to a structure",
    //   inputSchema: ApplyStrainSchema,
    //   annotations: {
    //     readOnlyHint: false,
    //     destructiveHint: false,
    //     idempotentHint: true,
    //     openWorldHint: false
    //   }
    // },
    // {
    //   name: "analyze_symmetry",
    //   description: "Analyze and detect symmetry properties of a crystal structure",
    //   inputSchema: AnalyzeSymmetrySchema,
    //   annotations: {
    //     readOnlyHint: true,
    //     destructiveHint: false,
    //     idempotentHint: true,
    //     openWorldHint: false
    //   }
    // },
    // {
    //   name: "validate_structure",
    //   description: "Validate a crystal structure for physical and chemical correctness",
    //   inputSchema: ValidateStructureSchema,
    //   annotations: {
    //     readOnlyHint: true,
    //     destructiveHint: false,
    //     idempotentHint: true,
    //     openWorldHint: false
    //   }
    // },
    // {
    //   name: "optimize_structure_mlff",
    //   description: "Optimize crystal structure using machine learning force fields (CHGNet, M3GNet, MACE)",
    //   inputSchema: OptimizeStructureMLFFSchema,
    //   annotations: {
    //     readOnlyHint: false,
    //     destructiveHint: false,
    //     idempotentHint: false,
    //     openWorldHint: false
    //   }
    // },
    // {
    //   name: "calculate_energy_mlff",
    //   description: "Calculate energy, forces, and stress using MLFF without optimization",
    //   inputSchema: CalculateEnergyMLFFSchema,
    //   annotations: {
    //     readOnlyHint: true,
    //     destructiveHint: false,
    //     idempotentHint: true,
    //     openWorldHint: false
    //   }
    // },
    // {
    //   name: "ground_state_search",
    //   description: "Search for ground state structure across multiple space groups using MLFF",
    //   inputSchema: GroundStateSearchSchema,
    //   annotations: {
    //     readOnlyHint: false,
    //     destructiveHint: false,
    //     idempotentHint: false,
    //     openWorldHint: false
    //   }
    // },
    // {
    //   name: "export_structure",
    //   description: "Export crystal structure to multiple file formats (CIF, POSCAR, XYZ, etc.)",
    //   inputSchema: ExportStructureSchema,
    //   annotations: {
    //     readOnlyHint: true,
    //     destructiveHint: false,
    //     idempotentHint: true,
    //     openWorldHint: false
    //   }
    // },
    // {
    //   name: "generate_visualization",
    //   description: "Generate interactive HTML (3Dmol.js) or static PNG visualization of a crystal structure",
    //   inputSchema: VisualizationSchema,
    //   annotations: {
    //     readOnlyHint: false,
    //     destructiveHint: false,
    //     idempotentHint: true,
    //     openWorldHint: false
    //   }
    // },
    // Advanced structure tools
    // {
    //   name: "generate_prototype",
    //   description: "Generate common prototype structures (rocksalt, perovskite, zincblende, wurtzite, etc.)",
    //   inputSchema: GeneratePrototypeSchema,
    //   annotations: {
    //     readOnlyHint: false,
    //     destructiveHint: false,
    //     idempotentHint: true,
    //     openWorldHint: true
    //   }
    // },
    // {
    //   name: "generate_twisted_bilayer",
    //   description: "Generate twisted bilayer/multilayer structures (graphene, MoS2, hBN) with specified twist angle",
    //   inputSchema: GenerateTwistedBilayerSchema,
    //   annotations: {
    //     readOnlyHint: false,
    //     destructiveHint: false,
    //     idempotentHint: true,
    //     openWorldHint: true
    //   }
    // },
    // {
    //   name: "generate_high_entropy_alloy",
    //   description: "Generate high-entropy alloy structures with 4+ elements in FCC/BCC/HCP lattices",
    //   inputSchema: GenerateHighEntropyAlloySchema,
    //   annotations: {
    //     readOnlyHint: false,
    //     destructiveHint: false,
    //     idempotentHint: true,
    //     openWorldHint: true
    //   }
    // },
    // {
    //   name: "generate_2d_material",
    //   description: "Generate 2D materials (hBN, MoS2, WS2, phosphorene, silicene, MXene)",
    //   inputSchema: Generate2DMaterialSchema,
    //   annotations: {
    //     readOnlyHint: false,
    //     destructiveHint: false,
    //     idempotentHint: true,
    //     openWorldHint: true
    //   }
    // },
    // {
    //   name: "generate_mof",
    //   description: "Generate metal-organic framework structures (MOF-5, HKUST-1, UiO-66, ZIF-8)",
    //   inputSchema: GenerateMOFSchema,
    //   annotations: {
    //     readOnlyHint: false,
    //     destructiveHint: false,
    //     idempotentHint: true,
    //     openWorldHint: true
    //   }
    // },
    // {
    //   name: "generate_cage",
    //   description: "Generate cage structures (fullerenes C60/C70/C80, clathrates)",
    //   inputSchema: GenerateCageSchema,
    //   annotations: {
    //     readOnlyHint: false,
    //     destructiveHint: false,
    //     idempotentHint: true,
    //     openWorldHint: true
    //   }
    // }
];
//# sourceMappingURL=tools.js.map