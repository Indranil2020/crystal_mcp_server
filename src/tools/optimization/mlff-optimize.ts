/**
 * Optimize Structure with MLFF Tool
 * 
 * Tool for optimizing crystal structures using Machine Learning Force Fields.
 */

import { OptimizeStructureMLFFSchema } from "../../types/tools.js";
import { Result, createSuccess, createFailure, createError, CrystalErrorCode, ERROR_MESSAGES } from "../../types/errors.js";
import { executePythonWithJSON } from "../../utils/python-bridge.js";
import { formatOptimizationOutput } from "../../utils/formatting.js";

export async function optimizeStructureMLFF(input: unknown): Promise<Result<any>> {
  const parsed = OptimizeStructureMLFFSchema.safeParse(input);
  if (!parsed.success) {
    return createFailure(createError(
      CrystalErrorCode.INVALID_INPUT,
      "Invalid input parameters",
      { zodErrors: parsed.error.errors },
      ["Check MLFF model name (chgnet, m3gnet, mace)", "Verify optimizer settings"],
      true
    ));
  }

  const params = { ...parsed.data, operation: "optimize", structure_dict: parsed.data.structure };

  const result = await executePythonWithJSON<typeof params, any>(
    "mlff_calculator.py",
    params,
    { timeout: 300000 } // 5 minute timeout for optimization
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
      ["Check MLFF model installation", "Try reducing number of optimization steps", "Consider using CPU instead of GPU"],
      false
    ));
  }

  return createSuccess(pythonResult);
}

export async function handleOptimizeStructureMLFF(args: unknown): Promise<any> {
  const result = await optimizeStructureMLFF(args);
  
  if (!result.success) {
    return {
      content: [{
        type: "text",
        text: `❌ **MLFF Optimization Failed**\n\n${result.error.message}\n\n**Suggestions:**\n${result.error.suggestions.map(s => `- ${s}`).join('\n')}`
      }],
      isError: true
    };
  }

  const outputText = formatOptimizationOutput(result.data);
  
  return {
    content: [{
      type: "text",
      text: outputText
    }]
  };
}
