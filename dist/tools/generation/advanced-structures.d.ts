/**
 * Advanced Structure Generation Tools
 *
 * TypeScript handlers for advanced structure generation.
 */
import { z } from "zod";
import { Result } from "../../types/errors.js";
export declare const GeneratePrototypeSchema: z.ZodObject<{
    prototype: z.ZodEnum<["rocksalt", "zincblende", "wurtzite", "fluorite", "antifluorite", "perovskite", "spinel", "heusler", "rutile", "diamond"]>;
    elements: z.ZodRecord<z.ZodString, z.ZodString>;
    lattice_constant: z.ZodNumber;
    c_over_a: z.ZodDefault<z.ZodNumber>;
}, "strip", z.ZodTypeAny, {
    elements: Record<string, string>;
    prototype: "rocksalt" | "zincblende" | "wurtzite" | "fluorite" | "antifluorite" | "perovskite" | "spinel" | "heusler" | "rutile" | "diamond";
    lattice_constant: number;
    c_over_a: number;
}, {
    elements: Record<string, string>;
    prototype: "rocksalt" | "zincblende" | "wurtzite" | "fluorite" | "antifluorite" | "perovskite" | "spinel" | "heusler" | "rutile" | "diamond";
    lattice_constant: number;
    c_over_a?: number | undefined;
}>;
export declare const GenerateTwistedBilayerSchema: z.ZodObject<{
    material: z.ZodDefault<z.ZodEnum<["graphene", "MoS2", "WS2", "hBN"]>>;
    twist_angle: z.ZodNumber;
    layers: z.ZodDefault<z.ZodNumber>;
    stacking: z.ZodDefault<z.ZodEnum<["AA", "AB"]>>;
    interlayer_distance: z.ZodDefault<z.ZodNumber>;
    vacuum: z.ZodDefault<z.ZodNumber>;
}, "strip", z.ZodTypeAny, {
    material: "graphene" | "MoS2" | "WS2" | "hBN";
    twist_angle: number;
    vacuum: number;
    layers: number;
    stacking: "AA" | "AB";
    interlayer_distance: number;
}, {
    twist_angle: number;
    material?: "graphene" | "MoS2" | "WS2" | "hBN" | undefined;
    vacuum?: number | undefined;
    layers?: number | undefined;
    stacking?: "AA" | "AB" | undefined;
    interlayer_distance?: number | undefined;
}>;
export declare const GenerateHighEntropyAlloySchema: z.ZodObject<{
    elements: z.ZodArray<z.ZodString, "many">;
    concentrations: z.ZodOptional<z.ZodArray<z.ZodNumber, "many">>;
    structure_type: z.ZodDefault<z.ZodEnum<["fcc", "bcc", "hcp"]>>;
    supercell: z.ZodDefault<z.ZodTuple<[z.ZodNumber, z.ZodNumber, z.ZodNumber], null>>;
    lattice_constant: z.ZodOptional<z.ZodNumber>;
    seed: z.ZodOptional<z.ZodNumber>;
}, "strip", z.ZodTypeAny, {
    elements: string[];
    supercell: [number, number, number];
    structure_type: "fcc" | "bcc" | "hcp";
    seed?: number | undefined;
    lattice_constant?: number | undefined;
    concentrations?: number[] | undefined;
}, {
    elements: string[];
    seed?: number | undefined;
    supercell?: [number, number, number] | undefined;
    lattice_constant?: number | undefined;
    concentrations?: number[] | undefined;
    structure_type?: "fcc" | "bcc" | "hcp" | undefined;
}>;
export declare const Generate2DMaterialSchema: z.ZodObject<{
    material: z.ZodEnum<["hBN", "MoS2", "WS2", "MoSe2", "WSe2", "phosphorene", "silicene", "MXene"]>;
    size: z.ZodDefault<z.ZodTuple<[z.ZodNumber, z.ZodNumber, z.ZodNumber], null>>;
    vacuum: z.ZodDefault<z.ZodNumber>;
    extra_params: z.ZodOptional<z.ZodRecord<z.ZodString, z.ZodAny>>;
}, "strip", z.ZodTypeAny, {
    material: "MoS2" | "WS2" | "hBN" | "MoSe2" | "WSe2" | "phosphorene" | "silicene" | "MXene";
    vacuum: number;
    size: [number, number, number];
    extra_params?: Record<string, any> | undefined;
}, {
    material: "MoS2" | "WS2" | "hBN" | "MoSe2" | "WSe2" | "phosphorene" | "silicene" | "MXene";
    vacuum?: number | undefined;
    size?: [number, number, number] | undefined;
    extra_params?: Record<string, any> | undefined;
}>;
export declare const GenerateMOFSchema: z.ZodObject<{
    mof_type: z.ZodEnum<["MOF-5", "HKUST-1", "UiO-66", "ZIF-8"]>;
    functionalization: z.ZodOptional<z.ZodString>;
    size: z.ZodDefault<z.ZodTuple<[z.ZodNumber, z.ZodNumber, z.ZodNumber], null>>;
}, "strip", z.ZodTypeAny, {
    size: [number, number, number];
    mof_type: "MOF-5" | "HKUST-1" | "UiO-66" | "ZIF-8";
    functionalization?: string | undefined;
}, {
    mof_type: "MOF-5" | "HKUST-1" | "UiO-66" | "ZIF-8";
    size?: [number, number, number] | undefined;
    functionalization?: string | undefined;
}>;
export declare const GenerateCageSchema: z.ZodObject<{
    cage_type: z.ZodEnum<["C60", "C70", "C80", "clathrate_I", "clathrate_II"]>;
    guest: z.ZodOptional<z.ZodString>;
    vacuum: z.ZodDefault<z.ZodNumber>;
}, "strip", z.ZodTypeAny, {
    vacuum: number;
    cage_type: "C60" | "C70" | "C80" | "clathrate_I" | "clathrate_II";
    guest?: string | undefined;
}, {
    cage_type: "C60" | "C70" | "C80" | "clathrate_I" | "clathrate_II";
    vacuum?: number | undefined;
    guest?: string | undefined;
}>;
export declare function generatePrototype(input: unknown): Promise<Result<any>>;
export declare function handleGeneratePrototype(args: unknown): Promise<any>;
export declare function generateTwistedBilayer(input: unknown): Promise<Result<any>>;
export declare function handleGenerateTwistedBilayer(args: unknown): Promise<any>;
export declare function generateHighEntropyAlloy(input: unknown): Promise<Result<any>>;
export declare function handleGenerateHighEntropyAlloy(args: unknown): Promise<any>;
export declare function generate2DMaterial(input: unknown): Promise<Result<any>>;
export declare function handleGenerate2DMaterial(args: unknown): Promise<any>;
export declare function generateMOF(input: unknown): Promise<Result<any>>;
export declare function handleGenerateMOF(args: unknown): Promise<any>;
export declare function generateCage(input: unknown): Promise<Result<any>>;
export declare function handleGenerateCage(args: unknown): Promise<any>;
//# sourceMappingURL=advanced-structures.d.ts.map