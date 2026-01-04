/**
 * Analyze Symmetry Tool
 * 
 * Tool for detecting and analyzing crystal symmetry.
 */

import { AnalyzeSymmetrySchema } from "../../types/tools.js";
import { Result, createSuccess, createFailure, createError, CrystalErrorCode } from "../../types/errors.js";
import { executePythonWithJSON } from "../../utils/python-bridge.js";
import { formatSymmetryOutput } from "../../utils/formatting.js";

export async function analyzeSymmetry(input: unknown): Promise<Result<any>> {
  const parsed = AnalyzeSymmetrySchema.safeParse(input);
  if (!parsed.success) {
    return createFailure(createError(
      CrystalErrorCode.INVALID_INPUT,
      "Invalid input parameters",
      { zodErrors: parsed.error.errors },
      ["Check symprec and angle_tolerance values"],
      true
    ));
  }

  const params = parsed.data;

  const result = await executePythonWithJSON<typeof params, any>(
    "symmetry_analyzer.py",
    params,
    { timeout: 60000 }
  );

  if (!result.success) {
    return createFailure(result.error);
  }

  const pythonResult = result.data;
  
  if (!pythonResult.success) {
    return createFailure(createError(
      pythonResult.error.code as CrystalErrorCode,
      pythonResult.error.message,
      pythonResult.error.details
    ));
  }

  return createSuccess(pythonResult);
}

export async function handleAnalyzeSymmetry(args: unknown): Promise<any> {
  const result = await analyzeSymmetry(args);
  
  if (!result.success) {
    return {
      content: [{
        type: "text",
        text: `**Symmetry Analysis Failed**\n\n${result.error.message}`
      }],
      isError: true
    };
  }

  const outputText = formatSymmetryOutput(result.data);
  
  // Include raw JSON data for the frontend viewer
  // Prefer conventional cell, then primitive cell, then null
  const structure = result.data.conventional_cell || result.data.primitive_cell;

  const jsonData = JSON.stringify({
      success: true,
      structure: structure,
      space_group: result.data.space_group,
      symmetry_operations: result.data.symmetry_operations
  });

  return {
    content: [
      {
        type: "text",
        text: outputText
      },
      {
        type: "text",
        text: `\n\n<json-data>\n${jsonData}\n</json-data>`
      }
    ]
  };
}
