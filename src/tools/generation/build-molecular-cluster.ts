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
    const parsed = BuildMolecularClusterSchema.safeParse(input);
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

    const result = await executePythonWithJSON<typeof parsed.data, any>(
        "molecular_cluster_generator.py",
        parsed.data,
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
                `- **Box Size**: ${structure.lattice.a.toFixed(2)} × ${structure.lattice.b.toFixed(2)} × ${structure.lattice.c.toFixed(2)} Å\n\n` +
                `*Cluster structure data is available in the response.*`
        }, {
            type: "text",
            text: `\n\n<json-data>\n${jsonData}\n</json-data>`
        }]
    };
}
