/**
 * Molecule Builder Tool
 *
 * Tool for generating isolated molecular structures.
 */
import { BuildMoleculeSchema } from "../../types/tools.js";
import { createSuccess, createFailure, createError, CrystalErrorCode } from "../../types/errors.js";
import { executePythonWithJSON } from "../../utils/python-bridge.js";
export async function buildMolecule(input) {
    const parsed = BuildMoleculeSchema.safeParse(input);
    if (!parsed.success) {
        return createFailure(createError(CrystalErrorCode.INVALID_INPUT, "Invalid input parameters", { zodErrors: parsed.error.errors }, ["Check molecule name"], true));
    }
    const result = await executePythonWithJSON("molecule_generator.py", parsed.data, { timeout: 30000 });
    if (!result.success) {
        return createFailure(result.error);
    }
    const pythonResult = result.data;
    if (!pythonResult.success) {
        return createFailure(createError(CrystalErrorCode.GENERATION_FAILED, pythonResult.error.message || "Failed to generate molecule", pythonResult.error.details));
    }
    return createSuccess(pythonResult);
}
export async function handleBuildMolecule(args) {
    const result = await buildMolecule(args);
    if (!result.success) {
        return {
            content: [{
                    type: "text",
                    text: `❌ **Molecule Building Failed**\n\n${result.error.message}`
                }],
            isError: true
        };
    }
    const structure = result.data.structure;
    const metadata = structure.metadata;
    return {
        content: [{
                type: "text",
                text: `### 🧪 Molecule Built: ${metadata.formula}\n\n` +
                    `- **Atoms**: ${metadata.natoms}\n` +
                    `- **Box Size**: ${structure.lattice.a.toFixed(2)} x ${structure.lattice.b.toFixed(2)} x ${structure.lattice.c.toFixed(2)} Å\n\n` +
                    `*Structure data is available in the response.*`
            }, {
                type: "text",
                text: JSON.stringify(structure, null, 2)
            }]
    };
}
//# sourceMappingURL=build-molecule.js.map