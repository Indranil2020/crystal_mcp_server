/**
 * Calculate Energy with MLFF Tool
 *
 * Tool for calculating energy, forces, and stress using Machine Learning Force Fields.
 */
import { CalculateEnergyMLFFSchema } from "../../types/tools.js";
import { createSuccess, createFailure, createError, CrystalErrorCode } from "../../types/errors.js";
import { executePythonWithJSON } from "../../utils/python-bridge.js";
export async function calculateEnergyMLFF(input) {
    const parsed = CalculateEnergyMLFFSchema.safeParse(input);
    if (!parsed.success) {
        return createFailure(createError(CrystalErrorCode.INVALID_INPUT, "Invalid input parameters", { zodErrors: parsed.error.errors }, ["Check MLFF model name (chgnet, m3gnet, mace)", "Verify structure input"], true));
    }
    const params = { ...parsed.data, operation: "calculate_energy", structure_dict: parsed.data.structure };
    const result = await executePythonWithJSON("mlff_calculator.py", params, { timeout: 60000 } // 1 minute timeout for energy calculation
    );
    if (!result.success) {
        return createFailure(result.error);
    }
    const pythonResult = result.data.data;
    if (!pythonResult.success) {
        return createFailure(createError(pythonResult.error.code, pythonResult.error.message, pythonResult.error.details, ["Check MLFF model installation", "Verify structure format", "Consider using CPU instead of GPU"], false));
    }
    return createSuccess(pythonResult);
}
export async function handleCalculateEnergyMLFF(args) {
    const result = await calculateEnergyMLFF(args);
    if (!result.success) {
        return {
            content: [{
                    type: "text",
                    text: `❌ **MLFF Energy Calculation Failed**\n\n${result.error.message}\n\n**Suggestions:**\n${result.error.suggestions.map(s => `- ${s}`).join('\n')}`
                }],
            isError: true
        };
    }
    const data = result.data.data;
    let outputText = `✅ **MLFF Energy Calculation Complete**\n\n`;
    outputText += `**Model:** ${data.model}\n`;
    outputText += `**Energy:** ${data.energy.toFixed(6)} eV\n`;
    if (data.energy_per_atom !== undefined) {
        outputText += `**Energy per atom:** ${data.energy_per_atom.toFixed(6)} eV/atom\n`;
    }
    if (data.forces && data.forces.length > 0) {
        const maxForce = Math.max(...data.forces.map((f) => Math.sqrt((f[0] ?? 0) ** 2 + (f[1] ?? 0) ** 2 + (f[2] ?? 0) ** 2)));
        outputText += `**Maximum force:** ${maxForce.toFixed(6)} eV/Å\n`;
    }
    if (data.stress) {
        outputText += `**Stress tensor available:** Yes\n`;
    }
    return {
        content: [{
                type: "text",
                text: outputText
            }]
    };
}
//# sourceMappingURL=calculate-energy-mlff.js.map