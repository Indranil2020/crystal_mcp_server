/**
 * Generate Slab Tool
 * 
 * Tool for creating surface slabs for DFT calculations.
 */

import { GenerateSlabSchema } from "../../types/tools.js";
import { Result, createSuccess, createFailure, createError, CrystalErrorCode } from "../../types/errors.js";
import { executePythonWithJSON } from "../../utils/python-bridge.js";
import { formatStructureOutput } from "../../utils/formatting.js";

export async function generateSlab(input: unknown): Promise<Result<any>> {
  const parsed = GenerateSlabSchema.safeParse(input);
  if (!parsed.success) {
    return createFailure(createError(
      CrystalErrorCode.INVALID_INPUT,
      "Invalid input parameters",
      { zodErrors: parsed.error.errors },
      ["Check miller_indices format: [h, k, l]", "Ensure thickness >= 1 and vacuum >= 0"],
      true
    ));
  }

  // Rename structure -> structure_dict for Python compatibility
  // Note: fix_atoms is filtered out as Python doesn't support it yet
  const { structure, fix_atoms, ...rest } = parsed.data;
  const params = { structure_dict: structure, ...rest, operation: "slab" };

  const result = await executePythonWithJSON<typeof params, any>(
    "structure_tools.py",
    params,
    { timeout: 120000 }
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

export async function handleGenerateSlab(args: unknown): Promise<any> {
  const result = await generateSlab(args);
  
  if (!result.success) {
    return {
      content: [{
        type: "text",
        text: `❌ **Slab Generation Failed**\n\n${result.error.message}`
      }],
      isError: true
    };
  }

  const data = result.data;
  let output = `## 📐 Surface Slab Generated\n\n`;
  output += `**Miller Indices:** (${data.miller_indices.join(' ')})\n`;
  output += `**Slab Thickness:** ${data.slab_thickness.toFixed(3)} Å\n`;
  output += `**Vacuum Thickness:** ${data.vacuum_thickness.toFixed(3)} Å\n`;
  output += `**Surface Area:** ${data.surface_area.toFixed(3)} ų\n`;
  if (data.fixed_atoms && data.fixed_atoms.length > 0) {
    output += `**Fixed Atoms:** ${data.fixed_atoms.length} atoms (bottom layers)\n`;
  }
  output += `\n${formatStructureOutput(data.slab, undefined)}`;
  
  return {
    content: [{
      type: "text",
      text: output
    }]
  };
}
