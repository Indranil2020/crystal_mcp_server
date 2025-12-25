/**
 * Calculate Energy with MLFF Tool
 *
 * Tool for calculating energy, forces, and stress using Machine Learning Force Fields.
 */

import { CalculateEnergyMLFFSchema } from "../../types/tools.js";
import { Result, createSuccess, createFailure, createError, CrystalErrorCode, ERROR_MESSAGES } from "../../types/errors.js";
import { executePythonWithJSON } from "../../utils/python-bridge.js";

export async function calculateEnergyMLFF(input: unknown): Promise<Result<any>> {
  const parsed = CalculateEnergyMLFFSchema.safeParse(input);
  if (!parsed.success) {
    return createFailure(createError(
      CrystalErrorCode.INVALID_INPUT,
      "Invalid input parameters",
      { zodErrors: parsed.error.errors },
      ["Check MLFF model name (chgnet, m3gnet, mace)", "Verify structure input"],
      true
    ));
  }

  const params = { ...parsed.data, operation: "calculate_energy", structure_dict: parsed.data.structure };

  const result = await executePythonWithJSON<typeof params, any>(
    "mlff_calculator.py",
    params,
    { timeout: 60000 } // 1 minute timeout for energy calculation
  );

  if (!result.success) {
    return createFailure(result.error);
  }

  const pythonResult = result.data;

  if (!pythonResult.success) {
    if (pythonResult.error?.code === "MODEL_NOT_AVAILABLE") {
      const model = pythonResult.error.details?.mlff_model ?? parsed.data.mlff_model;
      const messageInfo = ERROR_MESSAGES.MODEL_NOT_AVAILABLE(String(model));
      return createFailure(createError(
        CrystalErrorCode.MODEL_NOT_AVAILABLE,
        messageInfo.message,
        pythonResult.error.details ?? {},
        [...messageInfo.suggestions],
        false
      ));
    }

    return createFailure(createError(
      pythonResult.error.code as CrystalErrorCode,
      pythonResult.error.message,
      pythonResult.error.details,
      ["Check MLFF model installation", "Verify structure format", "Consider using CPU instead of GPU"],
      false
    ));
  }

  return createSuccess(pythonResult);
}

export async function handleCalculateEnergyMLFF(args: unknown): Promise<any> {
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

  const data = result.data;

  let outputText = `✅ **MLFF Energy Calculation Complete**\n\n`;
  outputText += `**Model:** ${data.model}\n`;
  outputText += `**Energy:** ${data.energy.toFixed(6)} eV\n`;

  if (data.energy_per_atom !== undefined) {
    outputText += `**Energy per atom:** ${data.energy_per_atom.toFixed(6)} eV/atom\n`;
  }

  if (data.forces && data.forces.length > 0) {
    const maxForce = Math.max(...data.forces.map((f: number[]) =>
      Math.sqrt((f[0] ?? 0)**2 + (f[1] ?? 0)**2 + (f[2] ?? 0)**2)
    ));
    outputText += `**Maximum force:** ${maxForce.toFixed(6)} eV/Angstrom\n`;
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
