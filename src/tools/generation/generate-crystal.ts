/**
 * Generate Crystal Tool
 *
 * Main tool for generating crystal structures with PyXtal.
 * Follows defensive programming with Result<T> pattern.
 */

import { GenerateCrystalSchema } from "../../types/tools.js";
import { CrystalStructure, StructureGenerationResult } from "../../types/crystal.js";
import { Result, createSuccess, createFailure, createError, CrystalErrorCode } from "../../types/errors.js";
import { executePythonWithJSON } from "../../utils/python-bridge.js";
import { formatStructureOutput } from "../../utils/formatting.js";
import { generateCIF, generatePOSCAR, generateXYZ, generateJSON } from "../../utils/structure-formats.js";

/**
 * Generate a crystal structure with specified composition and space group.
 * 
 * @param input - Generation parameters
 * @returns Result containing structure data or error
 */
export async function generateCrystal(input: unknown): Promise<Result<StructureGenerationResult>> {
  // Validate input schema
  const parsed = GenerateCrystalSchema.safeParse(input);
  if (!parsed.success) {
    return createFailure(createError(
      CrystalErrorCode.INVALID_INPUT,
      "Invalid input parameters",
      {
        zodErrors: parsed.error.errors,
        receivedInput: input
      },
      ["Check input parameters match the required schema", "See documentation for valid parameter types"],
      true
    ));
  }

  const params = parsed.data;

  // Execute Python crystal generator
  const result = await executePythonWithJSON<typeof params, any>(
    "crystal_generator.py",
    params,
    {
      timeout: 60000 // 60 second timeout
    }
  );

  if (!result.success) {
    return createFailure(result.error);
  }

  const pythonResult = result.data;

  // Check if Python execution succeeded
  if (!pythonResult.success) {
    const error = pythonResult.error;
    return createFailure(createError(
      error.code as CrystalErrorCode || CrystalErrorCode.GENERATION_FAILED,
      error.message,
      error.details,
      error.suggestions || [],
      error.recoverable !== false
    ));
  }

  // Extract structure data
  const structure = pythonResult.structure;
  const validation = pythonResult.validation;
  const metadata = pythonResult.metadata;

  // Create result
  const generationResult: StructureGenerationResult = {
    structure,
    validation: {
      valid: validation.valid,
      issues: validation.errors || [],
      warnings: validation.warnings || []
    },
    generation_metadata: {
      attempts: metadata?.attempts || 1,
      seed: params.seed,
      generation_time_ms: metadata?.generation_time_ms || 0
    },
    files: {
      cif: generateCIF(structure),
      poscar: generatePOSCAR(structure),
      xyz: generateXYZ(structure),
      json: generateJSON(structure)
    }
  };

  return createSuccess(generationResult);
}

/**
 * MCP tool handler for generate_crystal.
 */
export async function handleGenerateCrystal(args: unknown): Promise<any> {
  const result = await generateCrystal(args);
  
  if (!result.success) {
    return {
      content: [{
        type: "text",
        text: `❌ **Crystal Generation Failed**\n\n**Error:** ${result.error.message}\n\n**Code:** ${result.error.code}\n\n**Suggestions:**\n${result.error.suggestions.map(s => `- ${s}`).join('\n')}`
      }],
      isError: true
    };
  }
  
  const data = result.data;
  const outputText = formatStructureOutput(data.structure, data.validation);
  
  return {
    content: [{
      type: "text",
      text: outputText
    }]
  };
}
