/**
 * Validate Structure Tool
 * 
 * Tool for validating crystal structures.
 */

import { ValidateStructureSchema } from "../../types/tools.js";
import { Result, createSuccess, createFailure, createError, CrystalErrorCode } from "../../types/errors.js";
import { executePythonWithJSON } from "../../utils/python-bridge.js";
import { formatValidationOutput } from "../../utils/formatting.js";

export async function validateStructure(input: unknown): Promise<Result<any>> {
  const parsed = ValidateStructureSchema.safeParse(input);
  if (!parsed.success) {
    return createFailure(createError(
      CrystalErrorCode.INVALID_INPUT,
      "Invalid input parameters",
      { zodErrors: parsed.error.errors },
      ["Check structure and validation parameters"],
      true
    ));
  }

  const params = { structure_dict: parsed.data.structure, ...parsed.data };

  const result = await executePythonWithJSON<typeof params, any>(
    "validators.py",
    params,
    { timeout: 30000 }
  );

  if (!result.success) {
    return createFailure(result.error);
  }

  const pythonResult = result.data.data;
  
  if (!pythonResult.success) {
    return createFailure(createError(
      pythonResult.error.code as CrystalErrorCode,
      pythonResult.error.message,
      pythonResult.error.details
    ));
  }

  return createSuccess(pythonResult);
}

export async function handleValidateStructure(args: unknown): Promise<any> {
  const result = await validateStructure(args);
  
  if (!result.success) {
    return {
      content: [{
        type: "text",
        text: `❌ **Validation Failed**\n\n${result.error.message}`
      }],
      isError: true
    };
  }

  const outputText = formatValidationOutput(result.data);
  
  return {
    content: [{
      type: "text",
      text: outputText
    }]
  };
}
