/**
 * Molecular Cluster Builder Tool
 * 
 * Tool for generating molecular clusters (dimers, trimers, n-mers)
 * with various stacking arrangements for quantum chemistry.
 */

import { BuildMolecularClusterSchema } from "../../types/tools.js";
import { Result, createSuccess, createFailure, createError, CrystalErrorCode } from "../../types/errors.js";
import { executePythonWithJSON } from "../../utils/python-bridge.js";

export async function buildMolecularCluster(input: unknown): Promise<Result<any>> {
    // Sanitize input: remove null values as Zod optional() doesn't accept null
    const sanitizedInput = input && typeof input === 'object'
        ? Object.fromEntries(
            Object.entries(input as Record<string, unknown>)
                .filter(([_, v]) => v !== null)
        )
        : input;

    const parsed = BuildMolecularClusterSchema.safeParse(sanitizedInput);
    if (!parsed.success) {
        return createFailure(createError(
            CrystalErrorCode.INVALID_INPUT,
            "Invalid input parameters for molecular cluster",
            { zodErrors: parsed.error.errors },
            [
                "Ensure 'molecules' is an array with at least one molecule",
                "Each molecule needs an 'identifier' field",
                "Example: {molecules: [{identifier: 'benzene', count: 2}]}"
            ],
            true
        ));
    }

    // Normalize input: Handle array-based rotations (LLM quirk)
    // If rotation_x/y/z is an array, move it to 'rotations' structure
    const data = { ...parsed.data } as any;
    const axes = ['x', 'y', 'z'];

    // Initialize rotations array if needed
    if (!data.rotations) data.rotations = [];

    axes.forEach(axis => {
        const key = `rotation_${axis}`;
        const val = data[key];
        if (Array.isArray(val)) {
            // Merge into rotations list
            val.forEach((rotVal: number, i: number) => {
                if (!data.rotations[i]) data.rotations[i] = { x: 0, y: 0, z: 0 };
                data.rotations[i][axis] = rotVal;
            });
            // Clear the array field to avoid confusing backend which expects scalar
            delete data[key];
        }
    });

    const result = await executePythonWithJSON<any, any>(
        "molecular_cluster_generator.py",
        data,
        { timeout: 600000 }  // Longer timeout for cluster generation
    );

    if (!result.success) {
        return createFailure(result.error);
    }

    const pythonResult = result.data;

    if (!pythonResult.success) {
        return createFailure(createError(
            CrystalErrorCode.GENERATION_FAILED,
            pythonResult.error?.message || "Failed to generate molecular cluster",
            pythonResult.error?.details,
            [pythonResult.error?.suggestion || "Check molecule names and try again"]
        ));
    }

    return createSuccess(pythonResult);
}

export async function handleBuildMolecularCluster(args: unknown): Promise<any> {
    const result = await buildMolecularCluster(args);

    if (!result.success) {
        return {
            content: [{
                type: "text",
                text: `**Molecular Cluster Generation Failed**\n\n${result.error.message}`
            }],
            isError: true
        };
    }

    const data = result.data;
    const structure = data.structure;
    const metadata = structure?.metadata || {};

    // Format formulas
    const formulas = data.formulas?.join(" + ") || metadata.formula || "Unknown";

    // Format warnings
    const warnings = metadata.warnings || [];
    const warningText = warnings.length > 0
        ? `\n\n**⚠️ Warnings:**\n${warnings.map((w: string) => `- ${w}`).join("\n")}`
        : "";

    // Include raw JSON data for the frontend viewer
    const jsonData = JSON.stringify({
        success: true,
        structure: structure
    });

    return {
        content: [{
            type: "text",
            text: `### Molecular Cluster Built: ${formulas}\n\n` +
                `- **Molecules**: ${data.n_molecules}\n` +
                `- **Total Atoms**: ${data.n_atoms}\n` +
                `- **Stacking**: ${data.stacking_type || "auto"}\n` +
                `- **Distance**: ${data.intermolecular_distance?.toFixed(2) || "auto"} Å\n` +
                `- **Box Size**: ${structure.lattice.a.toFixed(2)} × ${structure.lattice.b.toFixed(2)} × ${structure.lattice.c.toFixed(2)} Å\n` +
                warningText + "\n\n" +
                `*Cluster structure data is available in the response.*`
        }, {
            type: "text",
            text: `\n\n**To visualize this structure:**\n` +
                `Use the \`generate_visualization\` tool with:\n` +
                `\`structure_file: "${data.auto_saved_path || data.output_file || 'STRUCTURE_FILE_PATH'}"\``
        }, {
            type: "text",
            text: `\n\n<json-data>\n${jsonData}\n</json-data>`
        }]
    };
}
