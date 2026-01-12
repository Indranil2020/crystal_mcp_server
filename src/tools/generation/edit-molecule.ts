/**
 * Molecule Editor Tool
 *
 * Tool for modifying molecular structures through natural language
 * or structured operations. Enables iterative molecular design workflows.
 */

import { EditMoleculeSchema } from "../../types/tools.js";
import { Result, createSuccess, createFailure, createError, CrystalErrorCode } from "../../types/errors.js";
import { executePythonWithJSON } from "../../utils/python-bridge.js";

interface EditMoleculeResult {
    success: boolean;
    structure?: {
        lattice: Record<string, unknown>;
        atoms: Array<Record<string, unknown>>;
        metadata: {
            formula: string;
            natoms: number;
            smiles: string;
            canonical_smiles?: string;
            molecular_weight?: number;
        };
    };
    annotations?: {
        stereocenters: Array<Record<string, unknown>>;
        rings: Array<Record<string, unknown>>;
        functional_groups: Array<Record<string, unknown>>;
        atom_labels: Record<string, string>;
        edit_history: Array<Record<string, unknown>>;
    };
    edit_result?: {
        operations_applied: Array<Record<string, unknown>>;
        successful_count: number;
        failed_count: number;
        validation?: {
            valid: boolean;
            issues: Array<Record<string, unknown>>;
            warnings_count: number;
            errors_count: number;
        };
    };
    warnings?: string[];
    error?: {
        code: string;
        message: string;
        details?: Record<string, unknown>;
    };
}

export async function editMolecule(input: unknown): Promise<Result<EditMoleculeResult>> {
    const parsed = EditMoleculeSchema.safeParse(input);
    if (!parsed.success) {
        return createFailure(createError(
            CrystalErrorCode.INVALID_INPUT,
            "Invalid input parameters for edit_molecule",
            { zodErrors: parsed.error.errors },
            ["Ensure 'molecule' and 'operations' are provided", "Check operation format"],
            true
        ));
    }

    console.error(`[MCP DEBUG] 🔄 Executing molecule_editor.py for editing molecule`);

    const result = await executePythonWithJSON<typeof parsed.data, EditMoleculeResult>(
        "molecule_editor.py",
        parsed.data,
        { timeout: 60000 }  // 60s for complex edits with optimization
    );

    if (!result.success) {
        console.error(`[MCP DEBUG] ❌ Python execution failed: ${result.error.message}`);
        return createFailure(result.error);
    }

    const pythonResult = result.data;
    console.error(`[MCP DEBUG] ✅ Python execution successful. Result success: ${pythonResult.success}`);

    if (!pythonResult.success) {
        console.error(`[MCP DEBUG] ❌ Python logical error: ${pythonResult.error?.message}`);
        return createFailure(createError(
            CrystalErrorCode.GENERATION_FAILED,
            pythonResult.error?.message || "Failed to edit molecule",
            pythonResult.error?.details
        ));
    }

    return createSuccess(pythonResult);
}

export async function handleEditMolecule(args: unknown): Promise<{
    content: Array<{ type: string; text: string }>;
    isError?: boolean;
}> {
    const result = await editMolecule(args);

    if (!result.success) {
        return {
            content: [{
                type: "text",
                text: `**Molecule Editing Failed**\n\n${result.error.message}`
            }],
            isError: true
        };
    }

    const data = result.data;
    const structure = data.structure;
    const metadata = structure?.metadata;
    const editResult = data.edit_result;

    // Build summary text
    const formula = metadata?.formula ?? 'Unknown';
    const natoms = metadata?.natoms ?? 'N/A';
    const smiles = metadata?.smiles ?? 'N/A';

    let summaryText = `### Molecule Edited: ${formula}\n\n`;
    summaryText += `- **Atoms**: ${natoms}\n`;
    summaryText += `- **SMILES**: \`${smiles}\`\n`;

    if (editResult) {
        summaryText += `- **Operations Applied**: ${editResult.successful_count} successful`;
        if (editResult.failed_count > 0) {
            summaryText += `, ${editResult.failed_count} failed`;
        }
        summaryText += `\n`;

        if (editResult.validation) {
            const validation = editResult.validation;
            if (validation.valid) {
                summaryText += `- **Validation**: ✓ Valid\n`;
            } else {
                summaryText += `- **Validation**: ✗ ${validation.errors_count} errors, ${validation.warnings_count} warnings\n`;
            }
        }
    }

    if (structure?.lattice) {
        const lattice = structure.lattice as { a?: number; b?: number; c?: number };
        if (lattice.a && lattice.b && lattice.c) {
            summaryText += `- **Box Size**: ${lattice.a.toFixed(2)} x ${lattice.b.toFixed(2)} x ${lattice.c.toFixed(2)} Å\n`;
        }
    }

    // Add warnings if any
    if (data.warnings && data.warnings.length > 0) {
        summaryText += `\n**Warnings:**\n`;
        for (const warning of data.warnings) {
            summaryText += `- ${warning}\n`;
        }
    }

    summaryText += `\n*Structure data and annotations are available in the response.*`;

    // Include raw JSON data for the frontend viewer
    const jsonData = JSON.stringify({
        success: true,
        structure: structure,
        annotations: data.annotations,
        edit_result: editResult
    });

    return {
        content: [
            {
                type: "text",
                text: summaryText
            },
            {
                type: "text",
                text: `\n\n<json-data>\n${jsonData}\n</json-data>`
            }
        ]
    };
}
