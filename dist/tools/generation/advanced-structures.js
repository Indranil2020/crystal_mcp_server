/**
 * Advanced Structure Generation Tools
 *
 * TypeScript handlers for advanced structure generation.
 */
import { z } from "zod";
import { createSuccess, createFailure, createError, CrystalErrorCode } from "../../types/errors.js";
import { executePythonWithJSON } from "../../utils/python-bridge.js";
import { formatStructureOutput } from "../../utils/formatting.js";
// ============================================================================
// SCHEMAS
// ============================================================================
export const GeneratePrototypeSchema = z.object({
    prototype: z.enum([
        "rocksalt", "zincblende", "wurtzite", "fluorite", "antifluorite",
        "perovskite", "spinel", "heusler", "rutile", "diamond"
    ]).describe("Prototype structure type"),
    elements: z.record(z.string(), z.string())
        .describe("Mapping of site labels to elements, e.g., {A: 'Ca', B: 'Ti', X: 'O'}"),
    lattice_constant: z.number().positive()
        .describe("Lattice constant 'a' in Angstroms"),
    c_over_a: z.number().positive().default(1.0)
        .describe("c/a ratio for non-cubic systems")
});
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
export const GenerateHighEntropyAlloySchema = z.object({
    elements: z.array(z.string()).min(4)
        .describe("List of 4+ elements for HEA"),
    concentrations: z.array(z.number()).optional()
        .describe("Concentrations (default: equimolar)"),
    structure_type: z.enum(["fcc", "bcc", "hcp"]).default("fcc")
        .describe("Base crystal structure"),
    supercell: z.tuple([z.number().int(), z.number().int(), z.number().int()]).default([3, 3, 3])
        .describe("Supercell size"),
    lattice_constant: z.number().positive().optional()
        .describe("Lattice constant (estimated if not provided)"),
    seed: z.number().int().optional()
        .describe("Random seed for reproducibility")
});
export const Generate2DMaterialSchema = z.object({
    material: z.enum(["hBN", "MoS2", "WS2", "MoSe2", "WSe2", "phosphorene", "silicene", "MXene"])
        .describe("2D material type"),
    size: z.tuple([z.number().int(), z.number().int(), z.number().int()]).default([1, 1, 1])
        .describe("Supercell size"),
    vacuum: z.number().positive().default(15.0)
        .describe("Vacuum padding"),
    extra_params: z.record(z.string(), z.any()).optional()
        .describe("Material-specific parameters")
});
export const GenerateMOFSchema = z.object({
    mof_type: z.enum(["MOF-5", "HKUST-1", "UiO-66", "ZIF-8"])
        .describe("MOF type"),
    functionalization: z.string().optional()
        .describe("Linker functionalization"),
    size: z.tuple([z.number().int(), z.number().int(), z.number().int()]).default([1, 1, 1])
        .describe("Supercell size")
});
export const GenerateCageSchema = z.object({
    cage_type: z.enum(["C60", "C70", "C80", "clathrate_I", "clathrate_II"])
        .describe("Cage structure type"),
    guest: z.string().optional()
        .describe("Guest atom for endohedral structures"),
    vacuum: z.number().positive().default(5.0)
});
// ============================================================================
// HANDLER FUNCTIONS
// ============================================================================
export async function generatePrototype(input) {
    const parsed = GeneratePrototypeSchema.safeParse(input);
    if (!parsed.success) {
        return createFailure(createError(CrystalErrorCode.INVALID_INPUT, "Invalid prototype parameters", { zodErrors: parsed.error.errors }));
    }
    const params = { ...parsed.data, operation: "prototype" };
    const result = await executePythonWithJSON("advanced_structures.py", params, { timeout: 60000 });
    if (!result.success)
        return createFailure(result.error);
    const pythonResult = result.data.data;
    if (!pythonResult.success) {
        return createFailure(createError(CrystalErrorCode.GENERATION_FAILED, pythonResult.error?.message || "Prototype generation failed"));
    }
    return createSuccess(pythonResult);
}
export async function handleGeneratePrototype(args) {
    const result = await generatePrototype(args);
    if (!result.success) {
        return {
            content: [{ type: "text", text: `❌ **Prototype Generation Failed**\n\n${result.error.message}` }],
            isError: true
        };
    }
    const data = result.data;
    const outputText = formatStructureOutput(data.structure, undefined);
    return {
        content: [{
                type: "text",
                text: `## 🏛️ ${data.prototype.toUpperCase()} Prototype Generated\n\n**Description:** ${data.description}\n**Space Group:** ${data.space_group}\n\n${outputText}`
            }]
    };
}
export async function generateTwistedBilayer(input) {
    const parsed = GenerateTwistedBilayerSchema.safeParse(input);
    if (!parsed.success) {
        return createFailure(createError(CrystalErrorCode.INVALID_INPUT, "Invalid twisted bilayer parameters", { zodErrors: parsed.error.errors }));
    }
    const params = { ...parsed.data, operation: "twisted_bilayer" };
    const result = await executePythonWithJSON("advanced_structures.py", params, { timeout: 120000 });
    if (!result.success)
        return createFailure(result.error);
    const pythonResult = result.data.data;
    if (!pythonResult.success) {
        return createFailure(createError(CrystalErrorCode.GENERATION_FAILED, pythonResult.error?.message || "Twisted bilayer generation failed"));
    }
    return createSuccess(pythonResult);
}
export async function handleGenerateTwistedBilayer(args) {
    const result = await generateTwistedBilayer(args);
    if (!result.success) {
        return {
            content: [{ type: "text", text: `❌ **Twisted Bilayer Generation Failed**\n\n${result.error.message}` }],
            isError: true
        };
    }
    const data = result.data;
    const warnings = data.warnings?.length > 0 ? `\n> ⚠️ ${data.warnings.join('\n> ')}\n` : "";
    const outputText = formatStructureOutput(data.structure, undefined);
    return {
        content: [{
                type: "text",
                text: `## 🔄 Twisted ${data.material} Bilayer\n\n**Twist Angle:** ${data.twist_angle}°\n**Commensurate Indices:** (${data.commensurate_indices[0]}, ${data.commensurate_indices[1]})\n**Atoms:** ${data.n_atoms}${warnings}\n\n${outputText}`
            }]
    };
}
export async function generateHighEntropyAlloy(input) {
    const parsed = GenerateHighEntropyAlloySchema.safeParse(input);
    if (!parsed.success) {
        return createFailure(createError(CrystalErrorCode.INVALID_INPUT, "Invalid HEA parameters", { zodErrors: parsed.error.errors }));
    }
    const params = { ...parsed.data, operation: "high_entropy_alloy" };
    const result = await executePythonWithJSON("advanced_structures.py", params, { timeout: 60000 });
    if (!result.success)
        return createFailure(result.error);
    const pythonResult = result.data.data;
    if (!pythonResult.success) {
        return createFailure(createError(CrystalErrorCode.GENERATION_FAILED, pythonResult.error?.message || "HEA generation failed"));
    }
    return createSuccess(pythonResult);
}
export async function handleGenerateHighEntropyAlloy(args) {
    const result = await generateHighEntropyAlloy(args);
    if (!result.success) {
        return {
            content: [{ type: "text", text: `❌ **HEA Generation Failed**\n\n${result.error.message}` }],
            isError: true
        };
    }
    const data = result.data;
    const comp = Object.entries(data.actual_composition).map(([e, n]) => `${e}:${n}`).join(", ");
    const outputText = formatStructureOutput(data.structure, undefined);
    return {
        content: [{
                type: "text",
                text: `## 🌈 High-Entropy Alloy Generated\n\n**Elements:** ${data.elements.join(", ")}\n**Structure:** ${data.structure_type.toUpperCase()}\n**Composition:** ${comp}\n**Atoms:** ${data.n_atoms}\n\n${outputText}`
            }]
    };
}
export async function generate2DMaterial(input) {
    const parsed = Generate2DMaterialSchema.safeParse(input);
    if (!parsed.success) {
        return createFailure(createError(CrystalErrorCode.INVALID_INPUT, "Invalid 2D material parameters", { zodErrors: parsed.error.errors }));
    }
    const params = { ...parsed.data, operation: "2d_material" };
    const result = await executePythonWithJSON("advanced_structures.py", params, { timeout: 60000 });
    if (!result.success)
        return createFailure(result.error);
    const pythonResult = result.data.data;
    if (!pythonResult.success) {
        return createFailure(createError(CrystalErrorCode.GENERATION_FAILED, pythonResult.error?.message || "2D material generation failed"));
    }
    return createSuccess(pythonResult);
}
export async function handleGenerate2DMaterial(args) {
    const result = await generate2DMaterial(args);
    if (!result.success) {
        return {
            content: [{ type: "text", text: `❌ **2D Material Generation Failed**\n\n${result.error.message}` }],
            isError: true
        };
    }
    const data = result.data;
    const outputText = data.structure ? formatStructureOutput(data.structure, undefined) : "";
    return {
        content: [{
                type: "text",
                text: `## 📄 2D Material: ${data.material}\n\n${outputText}`
            }]
    };
}
export async function generateMOF(input) {
    const parsed = GenerateMOFSchema.safeParse(input);
    if (!parsed.success) {
        return createFailure(createError(CrystalErrorCode.INVALID_INPUT, "Invalid MOF parameters", { zodErrors: parsed.error.errors }));
    }
    const params = { ...parsed.data, operation: "mof" };
    const result = await executePythonWithJSON("advanced_structures.py", params, { timeout: 60000 });
    if (!result.success)
        return createFailure(result.error);
    const pythonResult = result.data.data;
    if (!pythonResult.success) {
        return createFailure(createError(CrystalErrorCode.GENERATION_FAILED, pythonResult.error?.message || "MOF generation failed"));
    }
    return createSuccess(pythonResult);
}
export async function handleGenerateMOF(args) {
    const result = await generateMOF(args);
    if (!result.success) {
        return {
            content: [{ type: "text", text: `❌ **MOF Generation Failed**\n\n${result.error.message}` }],
            isError: true
        };
    }
    const data = result.data;
    const warnings = data.warnings?.length > 0 ? `\n> ⚠️ ${data.warnings.join('\n> ')}\n` : "";
    const outputText = formatStructureOutput(data.structure, undefined);
    return {
        content: [{
                type: "text",
                text: `## 🏗️ ${data.mof_type} Generated\n\n**Metal Node:** ${data.metal_node}\n**Linker:** ${data.linker}\n**Topology:** ${data.topology}\n**Space Group:** ${data.space_group}\n**Atoms:** ${data.n_atoms}${warnings}\n\n${outputText}`
            }]
    };
}
export async function generateCage(input) {
    const parsed = GenerateCageSchema.safeParse(input);
    if (!parsed.success) {
        return createFailure(createError(CrystalErrorCode.INVALID_INPUT, "Invalid cage parameters", { zodErrors: parsed.error.errors }));
    }
    const params = { ...parsed.data, operation: "cage", extra_params: { vacuum: parsed.data.vacuum } };
    const result = await executePythonWithJSON("advanced_structures.py", params, { timeout: 60000 });
    if (!result.success)
        return createFailure(result.error);
    const pythonResult = result.data.data;
    if (!pythonResult.success) {
        return createFailure(createError(CrystalErrorCode.GENERATION_FAILED, pythonResult.error?.message || "Cage generation failed"));
    }
    return createSuccess(pythonResult);
}
export async function handleGenerateCage(args) {
    const result = await generateCage(args);
    if (!result.success) {
        return {
            content: [{ type: "text", text: `❌ **Cage Generation Failed**\n\n${result.error.message}` }],
            isError: true
        };
    }
    const data = result.data;
    const guest = data.guest ? `\n**Guest:** ${data.guest}` : "";
    const outputText = formatStructureOutput(data.structure, undefined);
    return {
        content: [{
                type: "text",
                text: `## ⚽ Cage Structure: ${data.cage_type}${guest}\n\n**Atoms:** ${data.n_atoms}\n\n${outputText}`
            }]
    };
}
//# sourceMappingURL=advanced-structures.js.map