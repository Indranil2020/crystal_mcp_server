/**
 * Tool Schema Definitions for MCP Tools
 * 
 * This module defines input/output schemas and metadata for all MCP tools
 * exposed by the Crystal Structure Generation server.
 */

import { z } from "zod";
import type { CrystalSystem } from "./crystal.js";

/**
 * Tool annotations for MCP protocol
 */
export interface ToolAnnotations {
  readonly readOnlyHint: boolean;
  readonly destructiveHint: boolean;
  readonly idempotentHint: boolean;
  readonly openWorldHint: boolean;
}

/**
 * Common schema components
 */
const Vector3Schema = z.tuple([z.number(), z.number(), z.number()]);
const IntVector3Schema = z.tuple([z.number().int(), z.number().int(), z.number().int()]);
const Matrix3x3Schema = z.tuple([Vector3Schema, Vector3Schema, Vector3Schema]);

const ElementSymbolSchema = z.string().regex(/^[A-Z][a-z]?$/, "Invalid element symbol");
const ElementPairKeySchema = z.string().regex(
  /^[A-Z][a-z]?-[A-Z][a-z]?$/,
  "Expected element pair like 'Si-O'"
);

const StructureInputSchema = z.union([z.string(), z.any()]);

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

const CrystalSystemSchema: z.ZodType<CrystalSystem> = z.enum([
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
    .describe("Chemical composition as array of element symbols, e.g., ['Si', 'Si'] or ['Na', 'Cl']"),

  space_group: SpaceGroupSchema
    .describe("Space group number (1-230) or Hermann-Mauguin symbol"),

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

export type GenerateCrystalInput = z.infer<typeof GenerateCrystalSchema>;

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
}).passthrough();  // Allow additional operation-specific parameters to pass through

export type ComprehensiveGenerateInput = z.infer<typeof ComprehensiveGenerateSchema>;

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

export type SpaceGroupScanInput = z.infer<typeof SpaceGroupScanSchema>;

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

export type MakeSupercellInput = z.infer<typeof MakeSupercellSchema>;

/**
 * Schema for generate_slab tool
 */
export const GenerateSlabSchema = z.object({
  structure: StructureInputSchema
    .describe("Crystal structure as input"),

  miller_indices: Vector3Schema
    .describe("Miller indices [h, k, l] defining the surface plane"),

  thickness: z.number().int().positive()
    .describe("Slab thickness - used as minimum slab size in Angstroms (see min_slab_size for explicit Angstrom control)"),

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

export type GenerateSlabInput = z.infer<typeof GenerateSlabSchema>;

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

export type AnalyzeSymmetryInput = z.infer<typeof AnalyzeSymmetrySchema>;

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

export type ValidateStructureInput = z.infer<typeof ValidateStructureSchema>;

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

export type OptimizeStructureMLFFInput = z.infer<typeof OptimizeStructureMLFFSchema>;

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

export type CalculateEnergyMLFFInput = z.infer<typeof CalculateEnergyMLFFSchema>;

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

export type GroundStateSearchInput = z.infer<typeof GroundStateSearchSchema>;

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

export type ExportStructureInput = z.infer<typeof ExportStructureSchema>;

/**
 * Tool metadata for registration
 */
export interface ToolMetadata {
  readonly name: string;
  readonly description: string;
  readonly inputSchema: z.ZodSchema;
  readonly annotations: ToolAnnotations;
}

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

export type VisualizationInput = z.infer<typeof VisualizationSchema>;

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

export type CreateDefectInput = z.infer<typeof CreateDefectSchema>;


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

export type GenerateMolecularCrystalInput = z.infer<typeof GenerateMolecularCrystalSchema>;

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
    type: z.literal("nanotube").describe("Carbon nanotube"),
    params: NanotubeParamsSchema
  }),
  z.object({
    type: z.literal("graphene").describe("Graphene sheet"),
    params: GrapheneParamsSchema
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

export type GenerateNanostructureInput = z.infer<typeof GenerateNanostructureSchema>;

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
  if (data.operation === "get_subgroups" && !data.space_group) return false;
  if (data.operation === "get_path" && (!data.start_group || !data.end_group)) return false;
  return true;
}, {
  message: "Missing required parameters for selected operation"
});

export type ExploreSymmetryRelationsInput = z.infer<typeof ExploreSymmetryRelationsSchema>;

export const BuildMoleculeSchema = z.object({
  name: z.string().describe("Name of the molecule (e.g., 'H2O', 'C60', 'Benzene')"),
  vacuum: z.number().default(10.0).optional().describe("Vacuum padding around the molecule (Angstroms)")
});

export type BuildMoleculeInput = z.infer<typeof BuildMoleculeSchema>;


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

export type CreateAlloyInput = z.infer<typeof CreateAlloySchema>;

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

export type CreateHeterostructureInput = z.infer<typeof CreateHeterostructureSchema>;

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

export type AddAdsorbateInput = z.infer<typeof AddAdsorbateSchema>;

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

export type ApplyStrainInput = z.infer<typeof ApplyStrainSchema>;

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
  ]).describe("Prototype structure type"),
  elements: z.record(z.string(), ElementSymbolSchema)
    .describe("Mapping of site labels to elements, e.g., {A: 'Ca', B: 'Ti', X: 'O'}"),
  lattice_constant: z.number().positive().optional()
    .describe("Lattice constant 'a' in Angstroms"),
  c_over_a: z.number().positive().default(1.0)
    .describe("c/a ratio for non-cubic systems")
});
export type GeneratePrototypeInput = z.infer<typeof GeneratePrototypeSchema>;

/**
 * Generate twisted bilayer structure
 */
export const GenerateTwistedBilayerSchema = z.object({
  material: z.enum(["graphene", "MoS2", "WS2", "hBN"]).default("graphene")
    .describe("Base 2D material"),
  twist_angle: z.number().min(0).max(60)
    .describe("Twist angle in degrees"),
  layers: z.number().int().min(2).default(2)
    .describe("Number of layers"),
  stacking: z.enum(["AA", "AB"]).default("AB")
    .describe("Initial stacking"),
  interlayer_distance: z.number().positive().default(3.35)
    .describe("Interlayer spacing in Angstroms"),
  vacuum: z.number().positive().default(15.0)
    .describe("Vacuum padding in Angstroms")
});
export type GenerateTwistedBilayerInput = z.infer<typeof GenerateTwistedBilayerSchema>;

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
export type GenerateHighEntropyAlloyInput = z.infer<typeof GenerateHighEntropyAlloySchema>;

/**
 * Generate 2D material structure
 */
export const Generate2DMaterialSchema = z.object({
  material: z.enum(["hBN", "MoS2", "WS2", "MoSe2", "WSe2", "phosphorene", "silicene", "MXene"])
    .describe("2D material type"),
  size: IntVector3Schema.default([1, 1, 1])
    .describe("Supercell size"),
  vacuum: z.number().positive().default(15.0)
    .describe("Vacuum padding"),
  extra_params: z.record(z.string(), z.any()).optional()
    .describe("Material-specific parameters")
});
export type Generate2DMaterialInput = z.infer<typeof Generate2DMaterialSchema>;

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
export type GenerateMOFInput = z.infer<typeof GenerateMOFSchema>;

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
export type GenerateCageInput = z.infer<typeof GenerateCageSchema>;

/**
 * All tool definitions
 */
export const TOOL_DEFINITIONS: readonly ToolMetadata[] = [
  {
    name: "comprehensive_generate",
    description: "Unified generator for all crystal structures. Access 50+ operations across 18 categories: bulk, 2D, surface, molecule, twist, defect, electronic, etc.",
    inputSchema: ComprehensiveGenerateSchema,
    annotations: {
      readOnlyHint: false,
      destructiveHint: false,
      idempotentHint: false,
      openWorldHint: true
    }
  },
  {
    name: "build_molecule",
    description: "Generate an isolated molecular structure from a name (e.g., H2O, C60, Benzene).",
    inputSchema: BuildMoleculeSchema,
    annotations: {
      readOnlyHint: false,
      destructiveHint: false,
      idempotentHint: true,
      openWorldHint: true
    }
  },
  {
    name: "explore_symmetry_relations",
    description: "Explore group-subgroup relationships and find symmetry pathways between space groups.",
    inputSchema: ExploreSymmetryRelationsSchema,
    annotations: {
      readOnlyHint: true,
      destructiveHint: false,
      idempotentHint: true,
      openWorldHint: false
    }
  },
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
    name: "generate_molecular_crystal",
    description: "Generate a molecular crystal structure using PyXtal",
    inputSchema: GenerateMolecularCrystalSchema,
    annotations: {
      readOnlyHint: false,
      destructiveHint: false,
      idempotentHint: true,
      openWorldHint: true
    }
  },
  {
    name: "generate_nanostructure",
    description: "Generate nanostructures (nanotubes, graphene, ribbons, fullerenes)",
    inputSchema: GenerateNanostructureSchema,
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
    name: "create_defect",
    description: "Create point defects (vacancy, substitution, interstitial) in a crystal structure",
    inputSchema: CreateDefectSchema,
    annotations: {
      readOnlyHint: false,
      destructiveHint: false,
      idempotentHint: true,
      openWorldHint: false
    }
  },
  {
    name: "create_alloy",
    description: "Create a substitutional alloy with random mixing",
    inputSchema: CreateAlloySchema,
    annotations: {
      readOnlyHint: false,
      destructiveHint: false,
      idempotentHint: true,
      openWorldHint: false
    }
  },
  {
    name: "create_heterostructure",
    description: "Create a vertical heterostructure by stacking two layers",
    inputSchema: CreateHeterostructureSchema,
    annotations: {
      readOnlyHint: false,
      destructiveHint: false,
      idempotentHint: true,
      openWorldHint: false
    }
  },
  {
    name: "add_adsorbate",
    description: "Add an adsorbate molecule to a surface structure",
    inputSchema: AddAdsorbateSchema,
    annotations: {
      readOnlyHint: false,
      destructiveHint: false,
      idempotentHint: true,
      openWorldHint: false
    }
  },
  {
    name: "apply_strain",
    description: "Apply strain tensor to a structure",
    inputSchema: ApplyStrainSchema,
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
  },
  {
    name: "generate_visualization",
    description: "Generate interactive HTML (3Dmol.js) or static PNG visualization of a crystal structure",
    inputSchema: VisualizationSchema,
    annotations: {
      readOnlyHint: false,
      destructiveHint: false,
      idempotentHint: true,
      openWorldHint: false
    }
  },
  // Advanced structure tools
  {
    name: "generate_prototype",
    description: "Generate common prototype structures (rocksalt, perovskite, zincblende, wurtzite, etc.)",
    inputSchema: GeneratePrototypeSchema,
    annotations: {
      readOnlyHint: false,
      destructiveHint: false,
      idempotentHint: true,
      openWorldHint: true
    }
  },
  {
    name: "generate_twisted_bilayer",
    description: "Generate twisted bilayer/multilayer structures (graphene, MoS2, hBN) with specified twist angle",
    inputSchema: GenerateTwistedBilayerSchema,
    annotations: {
      readOnlyHint: false,
      destructiveHint: false,
      idempotentHint: true,
      openWorldHint: true
    }
  },
  {
    name: "generate_high_entropy_alloy",
    description: "Generate high-entropy alloy structures with 4+ elements in FCC/BCC/HCP lattices",
    inputSchema: GenerateHighEntropyAlloySchema,
    annotations: {
      readOnlyHint: false,
      destructiveHint: false,
      idempotentHint: true,
      openWorldHint: true
    }
  },
  {
    name: "generate_2d_material",
    description: "Generate 2D materials (hBN, MoS2, WS2, phosphorene, silicene, MXene)",
    inputSchema: Generate2DMaterialSchema,
    annotations: {
      readOnlyHint: false,
      destructiveHint: false,
      idempotentHint: true,
      openWorldHint: true
    }
  },
  {
    name: "generate_mof",
    description: "Generate metal-organic framework structures (MOF-5, HKUST-1, UiO-66, ZIF-8)",
    inputSchema: GenerateMOFSchema,
    annotations: {
      readOnlyHint: false,
      destructiveHint: false,
      idempotentHint: true,
      openWorldHint: true
    }
  },
  {
    name: "generate_cage",
    description: "Generate cage structures (fullerenes C60/C70/C80, clathrates)",
    inputSchema: GenerateCageSchema,
    annotations: {
      readOnlyHint: false,
      destructiveHint: false,
      idempotentHint: true,
      openWorldHint: true
    }
  }
] as const;
