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
 * Schema for generate_crystal tool
 */
export declare const GenerateCrystalSchema: z.ZodObject<{
    composition: z.ZodArray<z.ZodString, "many">;
    space_group: z.ZodUnion<[z.ZodNumber, z.ZodString]>;
    num_atoms: z.ZodOptional<z.ZodNumber>;
    lattice_params: z.ZodOptional<z.ZodObject<{
        a: z.ZodOptional<z.ZodNumber>;
        b: z.ZodOptional<z.ZodNumber>;
        c: z.ZodOptional<z.ZodNumber>;
        alpha: z.ZodOptional<z.ZodNumber>;
        beta: z.ZodOptional<z.ZodNumber>;
        gamma: z.ZodOptional<z.ZodNumber>;
    }, "strip", z.ZodTypeAny, {
        a?: number | undefined;
        b?: number | undefined;
        c?: number | undefined;
        alpha?: number | undefined;
        beta?: number | undefined;
        gamma?: number | undefined;
    }, {
        a?: number | undefined;
        b?: number | undefined;
        c?: number | undefined;
        alpha?: number | undefined;
        beta?: number | undefined;
        gamma?: number | undefined;
    }>>;
    volume_factor: z.ZodDefault<z.ZodNumber>;
    min_distance: z.ZodOptional<z.ZodRecord<z.ZodString, z.ZodNumber>>;
    wyckoff_positions: z.ZodOptional<z.ZodArray<z.ZodObject<{
        element: z.ZodString;
        wyckoff: z.ZodString;
        coords: z.ZodOptional<z.ZodTuple<[z.ZodNumber, z.ZodNumber, z.ZodNumber], null>>;
    }, "strip", z.ZodTypeAny, {
        element: string;
        wyckoff: string;
        coords?: [number, number, number] | undefined;
    }, {
        element: string;
        wyckoff: string;
        coords?: [number, number, number] | undefined;
    }>, "many">>;
    seed: z.ZodOptional<z.ZodNumber>;
    dimensionality: z.ZodDefault<z.ZodUnion<[z.ZodLiteral<0>, z.ZodLiteral<1>, z.ZodLiteral<2>, z.ZodLiteral<3>]>>;
    max_attempts: z.ZodDefault<z.ZodNumber>;
}, "strip", z.ZodTypeAny, {
    composition: string[];
    space_group: string | number;
    volume_factor: number;
    dimensionality: 0 | 3 | 2 | 1;
    max_attempts: number;
    num_atoms?: number | undefined;
    lattice_params?: {
        a?: number | undefined;
        b?: number | undefined;
        c?: number | undefined;
        alpha?: number | undefined;
        beta?: number | undefined;
        gamma?: number | undefined;
    } | undefined;
    min_distance?: Record<string, number> | undefined;
    wyckoff_positions?: {
        element: string;
        wyckoff: string;
        coords?: [number, number, number] | undefined;
    }[] | undefined;
    seed?: number | undefined;
}, {
    composition: string[];
    space_group: string | number;
    num_atoms?: number | undefined;
    lattice_params?: {
        a?: number | undefined;
        b?: number | undefined;
        c?: number | undefined;
        alpha?: number | undefined;
        beta?: number | undefined;
        gamma?: number | undefined;
    } | undefined;
    volume_factor?: number | undefined;
    min_distance?: Record<string, number> | undefined;
    wyckoff_positions?: {
        element: string;
        wyckoff: string;
        coords?: [number, number, number] | undefined;
    }[] | undefined;
    seed?: number | undefined;
    dimensionality?: 0 | 3 | 2 | 1 | undefined;
    max_attempts?: number | undefined;
}>;
export type GenerateCrystalInput = z.infer<typeof GenerateCrystalSchema>;
/**
 * Schema for comprehensive_generate tool
 *
 * Unified entry point for all 51+ generator operations across 18 categories:
 * bulk, two_d, surface, molecule, twist, defect, electronic, thermoelectric,
 * battery, catalyst, adsorption, magnetic, nanotube, quantum, photonic,
 * quality_control, high_pressure, external_fields
 */
export declare const ComprehensiveGenerateSchema: z.ZodObject<{
    operation: z.ZodString;
    category: z.ZodOptional<z.ZodString>;
    spacegroup: z.ZodOptional<z.ZodNumber>;
    elements: z.ZodOptional<z.ZodUnion<[z.ZodArray<z.ZodString, "many">, z.ZodRecord<z.ZodString, z.ZodString>]>>;
    composition: z.ZodOptional<z.ZodArray<z.ZodNumber, "many">>;
    material: z.ZodOptional<z.ZodString>;
    supercell: z.ZodOptional<z.ZodArray<z.ZodNumber, "many">>;
    prototype: z.ZodOptional<z.ZodString>;
    framework: z.ZodOptional<z.ZodString>;
    clathrate: z.ZodOptional<z.ZodString>;
    materials: z.ZodOptional<z.ZodArray<z.ZodString, "many">>;
    twist_angle: z.ZodOptional<z.ZodNumber>;
    n_atoms: z.ZodOptional<z.ZodNumber>;
    size_nm: z.ZodOptional<z.ZodNumber>;
    thickness_QL: z.ZodOptional<z.ZodNumber>;
    pressure_GPa: z.ZodOptional<z.ZodNumber>;
    list_available: z.ZodOptional<z.ZodBoolean>;
}, "strip", z.ZodTypeAny, {
    operation: string;
    composition?: number[] | undefined;
    category?: string | undefined;
    spacegroup?: number | undefined;
    elements?: Record<string, string> | string[] | undefined;
    material?: string | undefined;
    supercell?: number[] | undefined;
    prototype?: string | undefined;
    framework?: string | undefined;
    clathrate?: string | undefined;
    materials?: string[] | undefined;
    twist_angle?: number | undefined;
    n_atoms?: number | undefined;
    size_nm?: number | undefined;
    thickness_QL?: number | undefined;
    pressure_GPa?: number | undefined;
    list_available?: boolean | undefined;
}, {
    operation: string;
    composition?: number[] | undefined;
    category?: string | undefined;
    spacegroup?: number | undefined;
    elements?: Record<string, string> | string[] | undefined;
    material?: string | undefined;
    supercell?: number[] | undefined;
    prototype?: string | undefined;
    framework?: string | undefined;
    clathrate?: string | undefined;
    materials?: string[] | undefined;
    twist_angle?: number | undefined;
    n_atoms?: number | undefined;
    size_nm?: number | undefined;
    thickness_QL?: number | undefined;
    pressure_GPa?: number | undefined;
    list_available?: boolean | undefined;
}>;
export type ComprehensiveGenerateInput = z.infer<typeof ComprehensiveGenerateSchema>;
/**
 * Schema for space_group_scan tool
 */
export declare const SpaceGroupScanSchema: z.ZodObject<{
    composition: z.ZodArray<z.ZodString, "many">;
    space_groups: z.ZodOptional<z.ZodArray<z.ZodNumber, "many">>;
    space_group_range: z.ZodOptional<z.ZodTuple<[z.ZodNumber, z.ZodNumber], null>>;
    crystal_systems: z.ZodOptional<z.ZodArray<z.ZodType<CrystalSystem, z.ZodTypeDef, CrystalSystem>, "many">>;
    num_atoms: z.ZodOptional<z.ZodNumber>;
    volume_factor: z.ZodDefault<z.ZodNumber>;
    parallel: z.ZodDefault<z.ZodBoolean>;
    output_directory: z.ZodOptional<z.ZodString>;
    naming_scheme: z.ZodOptional<z.ZodString>;
}, "strip", z.ZodTypeAny, {
    composition: string[];
    volume_factor: number;
    parallel: boolean;
    num_atoms?: number | undefined;
    space_groups?: number[] | undefined;
    space_group_range?: [number, number] | undefined;
    crystal_systems?: CrystalSystem[] | undefined;
    output_directory?: string | undefined;
    naming_scheme?: string | undefined;
}, {
    composition: string[];
    num_atoms?: number | undefined;
    volume_factor?: number | undefined;
    space_groups?: number[] | undefined;
    space_group_range?: [number, number] | undefined;
    crystal_systems?: CrystalSystem[] | undefined;
    parallel?: boolean | undefined;
    output_directory?: string | undefined;
    naming_scheme?: string | undefined;
}>;
export type SpaceGroupScanInput = z.infer<typeof SpaceGroupScanSchema>;
/**
 * Schema for make_supercell tool
 */
export declare const MakeSupercellSchema: z.ZodObject<{
    structure: z.ZodUnion<[z.ZodString, z.ZodAny]>;
    matrix: z.ZodUnion<[z.ZodString, z.ZodArray<z.ZodArray<z.ZodNumber, "many">, "many">]>;
    wrap_atoms: z.ZodDefault<z.ZodBoolean>;
    preserve_symmetry: z.ZodDefault<z.ZodBoolean>;
    min_distance_check: z.ZodDefault<z.ZodBoolean>;
}, "strip", z.ZodTypeAny, {
    matrix: string | number[][];
    wrap_atoms: boolean;
    preserve_symmetry: boolean;
    min_distance_check: boolean;
    structure?: any;
}, {
    matrix: string | number[][];
    structure?: any;
    wrap_atoms?: boolean | undefined;
    preserve_symmetry?: boolean | undefined;
    min_distance_check?: boolean | undefined;
}>;
export type MakeSupercellInput = z.infer<typeof MakeSupercellSchema>;
/**
 * Schema for generate_slab tool
 */
export declare const GenerateSlabSchema: z.ZodObject<{
    structure: z.ZodUnion<[z.ZodString, z.ZodAny]>;
    miller_indices: z.ZodTuple<[z.ZodNumber, z.ZodNumber, z.ZodNumber], null>;
    thickness: z.ZodNumber;
    vacuum: z.ZodNumber;
    center_slab: z.ZodDefault<z.ZodBoolean>;
    fix_bottom_layers: z.ZodOptional<z.ZodNumber>;
    fix_atoms: z.ZodOptional<z.ZodArray<z.ZodNumber, "many">>;
    orthogonalize: z.ZodDefault<z.ZodBoolean>;
    min_slab_size: z.ZodOptional<z.ZodNumber>;
    symmetric: z.ZodDefault<z.ZodBoolean>;
}, "strip", z.ZodTypeAny, {
    miller_indices: [number, number, number];
    thickness: number;
    vacuum: number;
    center_slab: boolean;
    orthogonalize: boolean;
    symmetric: boolean;
    structure?: any;
    fix_bottom_layers?: number | undefined;
    fix_atoms?: number[] | undefined;
    min_slab_size?: number | undefined;
}, {
    miller_indices: [number, number, number];
    thickness: number;
    vacuum: number;
    structure?: any;
    center_slab?: boolean | undefined;
    fix_bottom_layers?: number | undefined;
    fix_atoms?: number[] | undefined;
    orthogonalize?: boolean | undefined;
    min_slab_size?: number | undefined;
    symmetric?: boolean | undefined;
}>;
export type GenerateSlabInput = z.infer<typeof GenerateSlabSchema>;
/**
 * Schema for analyze_symmetry tool
 */
export declare const AnalyzeSymmetrySchema: z.ZodObject<{
    structure: z.ZodUnion<[z.ZodString, z.ZodAny]>;
    symprec: z.ZodDefault<z.ZodNumber>;
    angle_tolerance: z.ZodDefault<z.ZodNumber>;
    detect_primitive: z.ZodDefault<z.ZodBoolean>;
    standardize: z.ZodDefault<z.ZodBoolean>;
}, "strip", z.ZodTypeAny, {
    symprec: number;
    angle_tolerance: number;
    detect_primitive: boolean;
    standardize: boolean;
    structure?: any;
}, {
    structure?: any;
    symprec?: number | undefined;
    angle_tolerance?: number | undefined;
    detect_primitive?: boolean | undefined;
    standardize?: boolean | undefined;
}>;
export type AnalyzeSymmetryInput = z.infer<typeof AnalyzeSymmetrySchema>;
/**
 * Schema for validate_structure tool
 */
export declare const ValidateStructureSchema: z.ZodObject<{
    structure: z.ZodUnion<[z.ZodString, z.ZodAny]>;
    checks: z.ZodDefault<z.ZodArray<z.ZodEnum<["distances", "symmetry", "stoichiometry", "overlaps", "cell_parameters", "charge_balance"]>, "many">>;
    min_distance: z.ZodOptional<z.ZodNumber>;
    max_coordination: z.ZodOptional<z.ZodNumber>;
    expected_space_group: z.ZodOptional<z.ZodNumber>;
}, "strip", z.ZodTypeAny, {
    checks: ("distances" | "symmetry" | "stoichiometry" | "overlaps" | "cell_parameters" | "charge_balance")[];
    min_distance?: number | undefined;
    structure?: any;
    max_coordination?: number | undefined;
    expected_space_group?: number | undefined;
}, {
    min_distance?: number | undefined;
    structure?: any;
    checks?: ("distances" | "symmetry" | "stoichiometry" | "overlaps" | "cell_parameters" | "charge_balance")[] | undefined;
    max_coordination?: number | undefined;
    expected_space_group?: number | undefined;
}>;
export type ValidateStructureInput = z.infer<typeof ValidateStructureSchema>;
/**
 * Schema for optimize_structure_mlff tool
 */
export declare const OptimizeStructureMLFFSchema: z.ZodObject<{
    structure: z.ZodUnion<[z.ZodString, z.ZodAny]>;
    mlff_model: z.ZodEnum<["chgnet", "m3gnet", "mace"]>;
    optimizer: z.ZodDefault<z.ZodEnum<["BFGS", "FIRE", "LBFGS", "GPMin"]>>;
    fmax: z.ZodDefault<z.ZodNumber>;
    steps: z.ZodDefault<z.ZodNumber>;
    constrain_symmetry: z.ZodDefault<z.ZodBoolean>;
    fix_lattice: z.ZodDefault<z.ZodBoolean>;
    fix_volume: z.ZodDefault<z.ZodBoolean>;
    pressure: z.ZodDefault<z.ZodNumber>;
    trajectory_file: z.ZodOptional<z.ZodString>;
}, "strip", z.ZodTypeAny, {
    mlff_model: "chgnet" | "m3gnet" | "mace";
    optimizer: "BFGS" | "FIRE" | "LBFGS" | "GPMin";
    fmax: number;
    steps: number;
    constrain_symmetry: boolean;
    fix_lattice: boolean;
    fix_volume: boolean;
    pressure: number;
    structure?: any;
    trajectory_file?: string | undefined;
}, {
    mlff_model: "chgnet" | "m3gnet" | "mace";
    structure?: any;
    optimizer?: "BFGS" | "FIRE" | "LBFGS" | "GPMin" | undefined;
    fmax?: number | undefined;
    steps?: number | undefined;
    constrain_symmetry?: boolean | undefined;
    fix_lattice?: boolean | undefined;
    fix_volume?: boolean | undefined;
    pressure?: number | undefined;
    trajectory_file?: string | undefined;
}>;
export type OptimizeStructureMLFFInput = z.infer<typeof OptimizeStructureMLFFSchema>;
/**
 * Schema for calculate_energy_mlff tool
 */
export declare const CalculateEnergyMLFFSchema: z.ZodObject<{
    structure: z.ZodUnion<[z.ZodString, z.ZodAny]>;
    mlff_model: z.ZodEnum<["chgnet", "m3gnet", "mace"]>;
    calculate_forces: z.ZodDefault<z.ZodBoolean>;
    calculate_stress: z.ZodDefault<z.ZodBoolean>;
}, "strip", z.ZodTypeAny, {
    mlff_model: "chgnet" | "m3gnet" | "mace";
    calculate_forces: boolean;
    calculate_stress: boolean;
    structure?: any;
}, {
    mlff_model: "chgnet" | "m3gnet" | "mace";
    structure?: any;
    calculate_forces?: boolean | undefined;
    calculate_stress?: boolean | undefined;
}>;
export type CalculateEnergyMLFFInput = z.infer<typeof CalculateEnergyMLFFSchema>;
/**
 * Schema for ground_state_search tool
 */
export declare const GroundStateSearchSchema: z.ZodObject<{
    composition: z.ZodArray<z.ZodString, "many">;
    space_groups: z.ZodOptional<z.ZodArray<z.ZodNumber, "many">>;
    num_structures_per_group: z.ZodDefault<z.ZodNumber>;
    mlff_model: z.ZodEnum<["chgnet", "m3gnet", "mace"]>;
    optimization_settings: z.ZodObject<{
        optimizer: z.ZodDefault<z.ZodEnum<["BFGS", "FIRE", "LBFGS"]>>;
        fmax: z.ZodDefault<z.ZodNumber>;
        steps: z.ZodDefault<z.ZodNumber>;
        constrain_symmetry: z.ZodDefault<z.ZodBoolean>;
    }, "strip", z.ZodTypeAny, {
        optimizer: "BFGS" | "FIRE" | "LBFGS";
        fmax: number;
        steps: number;
        constrain_symmetry: boolean;
    }, {
        optimizer?: "BFGS" | "FIRE" | "LBFGS" | undefined;
        fmax?: number | undefined;
        steps?: number | undefined;
        constrain_symmetry?: boolean | undefined;
    }>;
    parallel: z.ZodDefault<z.ZodBoolean>;
    n_workers: z.ZodOptional<z.ZodNumber>;
    output_directory: z.ZodOptional<z.ZodString>;
    save_trajectories: z.ZodDefault<z.ZodBoolean>;
}, "strip", z.ZodTypeAny, {
    composition: string[];
    parallel: boolean;
    mlff_model: "chgnet" | "m3gnet" | "mace";
    num_structures_per_group: number;
    optimization_settings: {
        optimizer: "BFGS" | "FIRE" | "LBFGS";
        fmax: number;
        steps: number;
        constrain_symmetry: boolean;
    };
    save_trajectories: boolean;
    space_groups?: number[] | undefined;
    output_directory?: string | undefined;
    n_workers?: number | undefined;
}, {
    composition: string[];
    mlff_model: "chgnet" | "m3gnet" | "mace";
    optimization_settings: {
        optimizer?: "BFGS" | "FIRE" | "LBFGS" | undefined;
        fmax?: number | undefined;
        steps?: number | undefined;
        constrain_symmetry?: boolean | undefined;
    };
    space_groups?: number[] | undefined;
    parallel?: boolean | undefined;
    output_directory?: string | undefined;
    num_structures_per_group?: number | undefined;
    n_workers?: number | undefined;
    save_trajectories?: boolean | undefined;
}>;
export type GroundStateSearchInput = z.infer<typeof GroundStateSearchSchema>;
/**
 * Schema for export_structure tool
 */
export declare const ExportStructureSchema: z.ZodObject<{
    structure: z.ZodUnion<[z.ZodString, z.ZodAny]>;
    formats: z.ZodArray<z.ZodEnum<["cif", "poscar", "contcar", "xyz", "extxyz", "json", "pdb", "pwscf", "castep", "lammps"]>, "many">;
    include_metadata: z.ZodDefault<z.ZodBoolean>;
    dft_software: z.ZodOptional<z.ZodEnum<["vasp", "qe", "cp2k", "castep"]>>;
    output_directory: z.ZodOptional<z.ZodString>;
}, "strip", z.ZodTypeAny, {
    formats: ("cif" | "poscar" | "contcar" | "xyz" | "extxyz" | "json" | "pdb" | "pwscf" | "castep" | "lammps")[];
    include_metadata: boolean;
    output_directory?: string | undefined;
    structure?: any;
    dft_software?: "castep" | "vasp" | "qe" | "cp2k" | undefined;
}, {
    formats: ("cif" | "poscar" | "contcar" | "xyz" | "extxyz" | "json" | "pdb" | "pwscf" | "castep" | "lammps")[];
    output_directory?: string | undefined;
    structure?: any;
    include_metadata?: boolean | undefined;
    dft_software?: "castep" | "vasp" | "qe" | "cp2k" | undefined;
}>;
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
export declare const VisualizationSchema: z.ZodObject<{
    structure: z.ZodAny;
    format: z.ZodDefault<z.ZodOptional<z.ZodEnum<["html", "png"]>>>;
    output_file: z.ZodOptional<z.ZodString>;
}, "strip", z.ZodTypeAny, {
    format: "html" | "png";
    structure?: any;
    output_file?: string | undefined;
}, {
    structure?: any;
    format?: "html" | "png" | undefined;
    output_file?: string | undefined;
}>;
export type VisualizationInput = z.infer<typeof VisualizationSchema>;
/**
 * Schema for create_defect tool
 */
export declare const CreateDefectSchema: z.ZodObject<{
    structure: z.ZodUnion<[z.ZodString, z.ZodAny]>;
    defect_type: z.ZodEnum<["vacancy", "substitution", "interstitial"]>;
    defect_site: z.ZodNumber;
    defect_species: z.ZodOptional<z.ZodString>;
    concentration: z.ZodDefault<z.ZodNumber>;
}, "strip", z.ZodTypeAny, {
    defect_type: "vacancy" | "substitution" | "interstitial";
    defect_site: number;
    concentration: number;
    structure?: any;
    defect_species?: string | undefined;
}, {
    defect_type: "vacancy" | "substitution" | "interstitial";
    defect_site: number;
    structure?: any;
    defect_species?: string | undefined;
    concentration?: number | undefined;
}>;
export type CreateDefectInput = z.infer<typeof CreateDefectSchema>;
/**
 * Schema for generate_molecular_crystal tool
 */
export declare const GenerateMolecularCrystalSchema: z.ZodObject<{
    molecules: z.ZodArray<z.ZodString, "many">;
    space_group: z.ZodUnion<[z.ZodNumber, z.ZodString]>;
    num_molecules: z.ZodOptional<z.ZodNumber>;
    lattice_params: z.ZodOptional<z.ZodObject<{
        a: z.ZodOptional<z.ZodNumber>;
        b: z.ZodOptional<z.ZodNumber>;
        c: z.ZodOptional<z.ZodNumber>;
        alpha: z.ZodOptional<z.ZodNumber>;
        beta: z.ZodOptional<z.ZodNumber>;
        gamma: z.ZodOptional<z.ZodNumber>;
    }, "strip", z.ZodTypeAny, {
        a?: number | undefined;
        b?: number | undefined;
        c?: number | undefined;
        alpha?: number | undefined;
        beta?: number | undefined;
        gamma?: number | undefined;
    }, {
        a?: number | undefined;
        b?: number | undefined;
        c?: number | undefined;
        alpha?: number | undefined;
        beta?: number | undefined;
        gamma?: number | undefined;
    }>>;
    volume_factor: z.ZodDefault<z.ZodNumber>;
    min_distance: z.ZodOptional<z.ZodRecord<z.ZodString, z.ZodNumber>>;
    seed: z.ZodOptional<z.ZodNumber>;
}, "strip", z.ZodTypeAny, {
    space_group: string | number;
    volume_factor: number;
    molecules: string[];
    lattice_params?: {
        a?: number | undefined;
        b?: number | undefined;
        c?: number | undefined;
        alpha?: number | undefined;
        beta?: number | undefined;
        gamma?: number | undefined;
    } | undefined;
    min_distance?: Record<string, number> | undefined;
    seed?: number | undefined;
    num_molecules?: number | undefined;
}, {
    space_group: string | number;
    molecules: string[];
    lattice_params?: {
        a?: number | undefined;
        b?: number | undefined;
        c?: number | undefined;
        alpha?: number | undefined;
        beta?: number | undefined;
        gamma?: number | undefined;
    } | undefined;
    volume_factor?: number | undefined;
    min_distance?: Record<string, number> | undefined;
    seed?: number | undefined;
    num_molecules?: number | undefined;
}>;
export type GenerateMolecularCrystalInput = z.infer<typeof GenerateMolecularCrystalSchema>;
/**
 * Schema for generate_nanostructure tool
 */
export declare const GenerateNanostructureSchema: z.ZodObject<{
    type: z.ZodEnum<["nanotube", "graphene", "nanoribbon", "fullerene", "mos2", "nanowire"]>;
    params: z.ZodObject<{
        n: z.ZodOptional<z.ZodNumber>;
        m: z.ZodOptional<z.ZodNumber>;
        length: z.ZodOptional<z.ZodNumber>;
        bond: z.ZodOptional<z.ZodNumber>;
        formula: z.ZodOptional<z.ZodString>;
        size: z.ZodOptional<z.ZodTuple<[z.ZodNumber, z.ZodNumber, z.ZodNumber], null>>;
        vacuum: z.ZodOptional<z.ZodNumber>;
        width: z.ZodOptional<z.ZodNumber>;
        ribbon_type: z.ZodOptional<z.ZodEnum<["armchair", "zigzag"]>>;
        saturated: z.ZodOptional<z.ZodBoolean>;
        name: z.ZodOptional<z.ZodString>;
        kind: z.ZodOptional<z.ZodString>;
        thickness: z.ZodOptional<z.ZodNumber>;
    }, "strip", z.ZodTypeAny, {
        length?: number | undefined;
        thickness?: number | undefined;
        vacuum?: number | undefined;
        n?: number | undefined;
        m?: number | undefined;
        bond?: number | undefined;
        formula?: string | undefined;
        size?: [number, number, number] | undefined;
        width?: number | undefined;
        ribbon_type?: "armchair" | "zigzag" | undefined;
        saturated?: boolean | undefined;
        name?: string | undefined;
        kind?: string | undefined;
    }, {
        length?: number | undefined;
        thickness?: number | undefined;
        vacuum?: number | undefined;
        n?: number | undefined;
        m?: number | undefined;
        bond?: number | undefined;
        formula?: string | undefined;
        size?: [number, number, number] | undefined;
        width?: number | undefined;
        ribbon_type?: "armchair" | "zigzag" | undefined;
        saturated?: boolean | undefined;
        name?: string | undefined;
        kind?: string | undefined;
    }>;
}, "strip", z.ZodTypeAny, {
    params: {
        length?: number | undefined;
        thickness?: number | undefined;
        vacuum?: number | undefined;
        n?: number | undefined;
        m?: number | undefined;
        bond?: number | undefined;
        formula?: string | undefined;
        size?: [number, number, number] | undefined;
        width?: number | undefined;
        ribbon_type?: "armchair" | "zigzag" | undefined;
        saturated?: boolean | undefined;
        name?: string | undefined;
        kind?: string | undefined;
    };
    type: "nanotube" | "graphene" | "nanoribbon" | "fullerene" | "mos2" | "nanowire";
}, {
    params: {
        length?: number | undefined;
        thickness?: number | undefined;
        vacuum?: number | undefined;
        n?: number | undefined;
        m?: number | undefined;
        bond?: number | undefined;
        formula?: string | undefined;
        size?: [number, number, number] | undefined;
        width?: number | undefined;
        ribbon_type?: "armchair" | "zigzag" | undefined;
        saturated?: boolean | undefined;
        name?: string | undefined;
        kind?: string | undefined;
    };
    type: "nanotube" | "graphene" | "nanoribbon" | "fullerene" | "mos2" | "nanowire";
}>;
export type GenerateNanostructureInput = z.infer<typeof GenerateNanostructureSchema>;
/**
 * Schema for explore_symmetry_relations tool
 */
export declare const ExploreSymmetryRelationsSchema: z.ZodEffects<z.ZodObject<{
    operation: z.ZodEnum<["get_subgroups", "get_path"]>;
    space_group: z.ZodOptional<z.ZodNumber>;
    start_group: z.ZodOptional<z.ZodNumber>;
    end_group: z.ZodOptional<z.ZodNumber>;
    max_depth: z.ZodOptional<z.ZodDefault<z.ZodNumber>>;
}, "strip", z.ZodTypeAny, {
    operation: "get_subgroups" | "get_path";
    space_group?: number | undefined;
    start_group?: number | undefined;
    end_group?: number | undefined;
    max_depth?: number | undefined;
}, {
    operation: "get_subgroups" | "get_path";
    space_group?: number | undefined;
    start_group?: number | undefined;
    end_group?: number | undefined;
    max_depth?: number | undefined;
}>, {
    operation: "get_subgroups" | "get_path";
    space_group?: number | undefined;
    start_group?: number | undefined;
    end_group?: number | undefined;
    max_depth?: number | undefined;
}, {
    operation: "get_subgroups" | "get_path";
    space_group?: number | undefined;
    start_group?: number | undefined;
    end_group?: number | undefined;
    max_depth?: number | undefined;
}>;
export type ExploreSymmetryRelationsInput = z.infer<typeof ExploreSymmetryRelationsSchema>;
export declare const BuildMoleculeSchema: z.ZodObject<{
    name: z.ZodString;
    vacuum: z.ZodOptional<z.ZodDefault<z.ZodNumber>>;
}, "strip", z.ZodTypeAny, {
    name: string;
    vacuum?: number | undefined;
}, {
    name: string;
    vacuum?: number | undefined;
}>;
export type BuildMoleculeInput = z.infer<typeof BuildMoleculeSchema>;
/**
 * Schema for create_alloy tool
 */
export declare const CreateAlloySchema: z.ZodObject<{
    structure: z.ZodUnion<[z.ZodString, z.ZodAny]>;
    substitutions: z.ZodRecord<z.ZodString, z.ZodObject<{
        element: z.ZodString;
        concentration: z.ZodNumber;
    }, "strip", z.ZodTypeAny, {
        element: string;
        concentration: number;
    }, {
        element: string;
        concentration: number;
    }>>;
    seed: z.ZodOptional<z.ZodNumber>;
}, "strip", z.ZodTypeAny, {
    substitutions: Record<string, {
        element: string;
        concentration: number;
    }>;
    seed?: number | undefined;
    structure?: any;
}, {
    substitutions: Record<string, {
        element: string;
        concentration: number;
    }>;
    seed?: number | undefined;
    structure?: any;
}>;
export type CreateAlloyInput = z.infer<typeof CreateAlloySchema>;
/**
 * Schema for create_heterostructure tool
 */
export declare const CreateHeterostructureSchema: z.ZodObject<{
    substrate: z.ZodUnion<[z.ZodString, z.ZodAny]>;
    overlayer: z.ZodUnion<[z.ZodString, z.ZodAny]>;
    interface_distance: z.ZodDefault<z.ZodNumber>;
    vacuum: z.ZodDefault<z.ZodNumber>;
}, "strip", z.ZodTypeAny, {
    vacuum: number;
    interface_distance: number;
    substrate?: any;
    overlayer?: any;
}, {
    vacuum?: number | undefined;
    substrate?: any;
    overlayer?: any;
    interface_distance?: number | undefined;
}>;
export type CreateHeterostructureInput = z.infer<typeof CreateHeterostructureSchema>;
/**
 * Schema for add_adsorbate tool
 */
export declare const AddAdsorbateSchema: z.ZodObject<{
    structure: z.ZodUnion<[z.ZodString, z.ZodAny]>;
    molecule: z.ZodUnion<[z.ZodString, z.ZodAny]>;
    site_index: z.ZodNumber;
    distance: z.ZodDefault<z.ZodNumber>;
}, "strip", z.ZodTypeAny, {
    site_index: number;
    distance: number;
    structure?: any;
    molecule?: any;
}, {
    site_index: number;
    structure?: any;
    molecule?: any;
    distance?: number | undefined;
}>;
export type AddAdsorbateInput = z.infer<typeof AddAdsorbateSchema>;
/**
 * Schema for apply_strain tool
 */
export declare const ApplyStrainSchema: z.ZodObject<{
    structure: z.ZodUnion<[z.ZodString, z.ZodAny]>;
    strain_tensor: z.ZodOptional<z.ZodArray<z.ZodNumber, "many">>;
    strain_type: z.ZodOptional<z.ZodEnum<["biaxial", "uniaxial", "hydrostatic"]>>;
    strain_value: z.ZodOptional<z.ZodNumber>;
}, "strip", z.ZodTypeAny, {
    structure?: any;
    strain_tensor?: number[] | undefined;
    strain_type?: "biaxial" | "uniaxial" | "hydrostatic" | undefined;
    strain_value?: number | undefined;
}, {
    structure?: any;
    strain_tensor?: number[] | undefined;
    strain_type?: "biaxial" | "uniaxial" | "hydrostatic" | undefined;
    strain_value?: number | undefined;
}>;
export type ApplyStrainInput = z.infer<typeof ApplyStrainSchema>;
/**
 * All tool definitions
 */
export declare const TOOL_DEFINITIONS: readonly ToolMetadata[];
//# sourceMappingURL=tools.d.ts.map