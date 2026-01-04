/**
 * Advanced Structure Generation Tools
 *
 * TypeScript handlers for advanced structure generation.
 */
import { createSuccess, createFailure, createError, CrystalErrorCode } from "../../types/errors.js";
import { executePythonWithJSON } from "../../utils/python-bridge.js";
import { formatStructureOutput } from "../../utils/formatting.js";
import { Generate2DMaterialSchema, GenerateCageSchema, GenerateHighEntropyAlloySchema, GenerateMOFSchema, GeneratePrototypeSchema, GenerateTwistedBilayerSchema } from "../../types/tools.js";
// ============================================================================
// SCHEMAS
// ============================================================================
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
    const pythonResult = result.data;
    if (!pythonResult.success) {
        return createFailure(createError(CrystalErrorCode.GENERATION_FAILED, pythonResult.error?.message || "Prototype generation failed"));
    }
    return createSuccess(pythonResult);
}
export async function handleGeneratePrototype(args) {
    const result = await generatePrototype(args);
    if (!result.success) {
        return {
            content: [{ type: "text", text: `**Prototype Generation Failed**\n\n${result.error.message}` }],
            isError: true
        };
    }
    const data = result.data;
    const outputText = formatStructureOutput(data.structure, undefined);
    const jsonData = JSON.stringify({
        success: true,
        structure: data.structure,
        prototype: data.prototype,
        description: data.description,
        space_group: data.space_group
    });
    return {
        content: [
            {
                type: "text",
                text: `## ${data.prototype.toUpperCase()} Prototype Generated\n\n**Description:** ${data.description}\n**Space Group:** ${data.space_group}\n\n${outputText}`
            },
            {
                type: "text",
                text: `\n\n<json-data>\n${jsonData}\n</json-data>`
            }
        ]
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
    const pythonResult = result.data;
    if (!pythonResult.success) {
        return createFailure(createError(CrystalErrorCode.GENERATION_FAILED, pythonResult.error?.message || "Twisted bilayer generation failed"));
    }
    return createSuccess(pythonResult);
}
export async function handleGenerateTwistedBilayer(args) {
    const result = await generateTwistedBilayer(args);
    if (!result.success) {
        return {
            content: [{ type: "text", text: `**Twisted Bilayer Generation Failed**\n\n${result.error.message}` }],
            isError: true
        };
    }
    const data = result.data;
    const warnings = data.warnings?.length > 0 ? `\n> ⚠️ ${data.warnings.join('\n> ')}\n` : "";
    const outputText = formatStructureOutput(data.structure, undefined);
    // Include raw JSON data for the frontend viewer
    const jsonData = JSON.stringify({
        success: true,
        structure: data.structure,
        material: data.material,
        twist_angle: data.twist_angle,
        commensurate_indices: data.commensurate_indices,
        n_atoms: data.n_atoms
    });
    return {
        content: [
            {
                type: "text",
                text: `## 🔄 Twisted ${data.material} Bilayer\n\n**Twist Angle:** ${data.twist_angle}°\n**Commensurate Indices:** (${data.commensurate_indices[0]}, ${data.commensurate_indices[1]})\n**Atoms:** ${data.n_atoms}${warnings}\n\n${outputText}`
            },
            {
                type: "text",
                text: `\n\n<json-data>\n${jsonData}\n</json-data>`
            }
        ]
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
    const pythonResult = result.data;
    if (!pythonResult.success) {
        return createFailure(createError(CrystalErrorCode.GENERATION_FAILED, pythonResult.error?.message || "HEA generation failed"));
    }
    return createSuccess(pythonResult);
}
export async function handleGenerateHighEntropyAlloy(args) {
    const result = await generateHighEntropyAlloy(args);
    if (!result.success) {
        return {
            content: [{ type: "text", text: `**HEA Generation Failed**\n\n${result.error.message}` }],
            isError: true
        };
    }
    const data = result.data;
    const comp = Object.entries(data.actual_composition).map(([e, n]) => `${e}:${n}`).join(", ");
    const outputText = formatStructureOutput(data.structure, undefined);
    // Include raw JSON data for the frontend viewer
    const jsonData = JSON.stringify({
        success: true,
        structure: data.structure,
        elements: data.elements,
        composition: data.actual_composition,
        structure_type: data.structure_type
    });
    return {
        content: [
            {
                type: "text",
                text: `## 🌈 High-Entropy Alloy Generated\n\n**Elements:** ${data.elements.join(", ")}\n**Structure:** ${data.structure_type.toUpperCase()}\n**Composition:** ${comp}\n**Atoms:** ${data.n_atoms}\n\n${outputText}`
            },
            {
                type: "text",
                text: `\n\n<json-data>\n${jsonData}\n</json-data>`
            }
        ]
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
    const pythonResult = result.data;
    if (!pythonResult.success) {
        return createFailure(createError(CrystalErrorCode.GENERATION_FAILED, pythonResult.error?.message || "2D material generation failed"));
    }
    return createSuccess(pythonResult);
}
export async function handleGenerate2DMaterial(args) {
    const result = await generate2DMaterial(args);
    if (!result.success) {
        return {
            content: [{ type: "text", text: `**2D Material Generation Failed**\n\n${result.error.message}` }],
            isError: true
        };
    }
    const data = result.data;
    const outputText = data.structure ? formatStructureOutput(data.structure, undefined) : "";
    // Include raw JSON data for the frontend viewer
    const jsonData = JSON.stringify({
        success: true,
        structure: data.structure,
        material: data.material
    });
    return {
        content: [
            {
                type: "text",
                text: `## 📄 2D Material: ${data.material}\n\n${outputText}`
            },
            {
                type: "text",
                text: `\n\n<json-data>\n${jsonData}\n</json-data>`
            }
        ]
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
    const pythonResult = result.data;
    if (!pythonResult.success) {
        return createFailure(createError(CrystalErrorCode.GENERATION_FAILED, pythonResult.error?.message || "MOF generation failed"));
    }
    return createSuccess(pythonResult);
}
export async function handleGenerateMOF(args) {
    const result = await generateMOF(args);
    if (!result.success) {
        return {
            content: [{ type: "text", text: `**MOF Generation Failed**\n\n${result.error.message}` }],
            isError: true
        };
    }
    const data = result.data;
    const warnings = data.warnings?.length > 0 ? `\n> ⚠️ ${data.warnings.join('\n> ')}\n` : "";
    const outputText = formatStructureOutput(data.structure, undefined);
    // Include raw JSON data for the frontend viewer
    const jsonData = JSON.stringify({
        success: true,
        structure: data.structure,
        mof_type: data.mof_type,
        metal_node: data.metal_node,
        linker: data.linker,
        topology: data.topology
    });
    return {
        content: [
            {
                type: "text",
                text: `## 🏗️ ${data.mof_type} Generated\n\n**Metal Node:** ${data.metal_node}\n**Linker:** ${data.linker}\n**Topology:** ${data.topology}\n**Space Group:** ${data.space_group}\n**Atoms:** ${data.n_atoms}${warnings}\n\n${outputText}`
            },
            {
                type: "text",
                text: `\n\n<json-data>\n${jsonData}\n</json-data>`
            }
        ]
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
    const pythonResult = result.data;
    if (!pythonResult.success) {
        return createFailure(createError(CrystalErrorCode.GENERATION_FAILED, pythonResult.error?.message || "Cage generation failed"));
    }
    return createSuccess(pythonResult);
}
export async function handleGenerateCage(args) {
    const result = await generateCage(args);
    if (!result.success) {
        return {
            content: [{ type: "text", text: `**Cage Generation Failed**\n\n${result.error.message}` }],
            isError: true
        };
    }
    const data = result.data;
    const guest = data.guest ? `\n**Guest:** ${data.guest}` : "";
    const outputText = formatStructureOutput(data.structure, undefined);
    // Include raw JSON data for the frontend viewer
    const jsonData = JSON.stringify({
        success: true,
        structure: data.structure,
        cage_type: data.cage_type,
        guest: data.guest
    });
    return {
        content: [
            {
                type: "text",
                text: `## ⚽ Cage Structure: ${data.cage_type}${guest}\n\n**Atoms:** ${data.n_atoms}\n\n${outputText}`
            },
            {
                type: "text",
                text: `\n\n<json-data>\n${jsonData}\n</json-data>`
            }
        ]
    };
}
//# sourceMappingURL=advanced-structures.js.map