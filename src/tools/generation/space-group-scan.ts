/**
 * Space Group Scan Tool
 * 
 * Tool for scanning multiple space groups to find stable structures.
 */

import { SpaceGroupScanSchema } from "../../types/tools.js";
import { Result, createSuccess, createFailure, createError, CrystalErrorCode } from "../../types/errors.js";
import { executePythonWithJSON } from "../../utils/python-bridge.js";
import { formatSpaceGroupScanOutput } from "../../utils/formatting.js";

export async function spaceGroupScan(input: unknown): Promise<Result<any>> {
  const parsed = SpaceGroupScanSchema.safeParse(input);
  if (!parsed.success) {
    return createFailure(createError(
      CrystalErrorCode.INVALID_INPUT,
      "Invalid input parameters",
      { zodErrors: parsed.error.errors },
      ["Check space_groups array or space_group_range parameters"],
      true
    ));
  }

  const params = parsed.data;

  // If no space groups specified, use all 230
  if (!params.space_groups && !params.space_group_range) {
    params.space_groups = Array.from({ length: 230 }, (_, i) => i + 1);
  }

  const result = await executePythonWithJSON<typeof params, any>(
    "space_group_scanner.py",
    params,
    { timeout: 600000 } // 10 minute timeout for large scans
  );

  if (!result.success) {
    return createFailure(result.error);
  }

  const scanResults = result.data;
  return createSuccess(scanResults);
}

export async function handleSpaceGroupScan(args: unknown): Promise<any> {
  const result = await spaceGroupScan(args);
  
  if (!result.success) {
    return {
      content: [{
        type: "text",
        text: `❌ **Space Group Scan Failed**\n\n${result.error.message}`
      }],
      isError: true
    };
  }

  const outputText = formatSpaceGroupScanOutput(result.data.generated_structures);
  
  return {
    content: [{
      type: "text",
      text: outputText
    }]
  };
}
