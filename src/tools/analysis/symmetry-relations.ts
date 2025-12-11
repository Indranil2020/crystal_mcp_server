/**
 * Symmetry Relations Analysis Tool
 * 
 * Tool for exploring group-subgroup relationships and symmetry paths.
 */

import { ExploreSymmetryRelationsSchema } from "../../types/tools.js";
import { Result, createSuccess, createFailure, createError, CrystalErrorCode } from "../../types/errors.js";
import { executePythonWithJSON } from "../../utils/python-bridge.js";

function formatSymmetryRelationsOutput(result: any): string {
    let output = `## 🕸️ Symmetry Relations Analysis\n\n`;

    if (result.subgroups) {
        output += `### Maximal Subgroups of Group ${result.space_group.number} (${result.space_group.symbol})\n\n`;

        // Group by crystal system for better readability
        const bySystem: Record<string, any[]> = {};
        result.subgroups.forEach((sg: any) => {
            const sys = sg.crystal_system || "Unknown";
            if (!bySystem[sys]) bySystem[sys] = [];
            bySystem[sys].push(sg);
        });

        output += `| Crystal System | Space Group | Symbol | Point Group |\n`;
        output += `|----------------|-------------|--------|-------------|\n`;

        Object.keys(bySystem).sort().forEach(sys => {
            bySystem[sys]?.forEach(sg => {
                output += `| ${sys} | ${sg.number} | ${sg.symbol} | ${sg.point_group} |\n`;
            });
        });

        output += `\n**Total Unique Maximal Subgroups:** ${result.subgroups.length}\n`;
    }

    if (result.path) {
        output += `### Symmetry Path\n\n`;
        const pathStr = result.path.map((p: any) => `**${p.number}** (${p.symbol})`).join(" → ");
        output += `${pathStr}\n\n`;
        output += `**Path Length:** ${result.length} steps\n`;
    }

    return output;
}

export async function exploreSymmetryRelations(input: unknown): Promise<Result<any>> {
    const parsed = ExploreSymmetryRelationsSchema.safeParse(input);
    if (!parsed.success) {
        return createFailure(createError(
            CrystalErrorCode.INVALID_INPUT,
            "Invalid input parameters",
            { zodErrors: parsed.error.errors },
            ["Check space group numbers (1-230)"],
            true
        ));
    }

    const result = await executePythonWithJSON<typeof parsed.data, any>(
        "symmetry_transformations.py",
        parsed.data,
        { timeout: 30000 }
    );

    if (!result.success) {
        return createFailure(result.error);
    }

    const pythonResult = result.data.data;

    if (!pythonResult.success) {
        return createFailure(createError(
            CrystalErrorCode.PYTHON_EXECUTION_FAILED,
            pythonResult.error.message,
            pythonResult.error.details
        ));
    }

    return createSuccess(pythonResult);
}

export async function handleExploreSymmetryRelations(args: unknown): Promise<any> {
    const result = await exploreSymmetryRelations(args);

    if (!result.success) {
        return {
            content: [{
                type: "text",
                text: `❌ **Symmetry Analysis Failed**\n\n${result.error.message}`
            }],
            isError: true
        };
    }

    const outputText = formatSymmetryRelationsOutput(result.data);

    return {
        content: [{
            type: "text",
            text: outputText
        }]
    };
}
