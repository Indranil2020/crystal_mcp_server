/**
 * Export Structure Tool
 *
 * Tool for exporting crystal structures to various file formats.
 */

import { ExportStructureSchema } from "../../types/tools.js";
import { Result, createSuccess, createFailure, createError, CrystalErrorCode } from "../../types/errors.js";
import { writeFileSync } from "../../utils/file-io.js";
import { generateCIF, generatePOSCAR, generateXYZ, generateJSON } from "../../utils/structure-formats.js";

export async function exportStructure(input: unknown): Promise<Result<any>> {
  const parsed = ExportStructureSchema.safeParse(input);
  if (!parsed.success) {
    return createFailure(createError(
      CrystalErrorCode.INVALID_INPUT,
      "Invalid input parameters",
      { zodErrors: parsed.error.errors },
      ["Check formats array and structure object"],
      true
    ));
  }

  const { structure, formats, output_directory } = parsed.data;

  const exportedFiles: Record<string, string> = {};
  const filePaths: Record<string, string> = {};

  for (const format of formats) {
    let content = '';
    let filename = `structure.${format}`;

    if (format === 'cif') {
      content = generateCIF(structure);
    } else if (format === 'poscar') {
      content = generatePOSCAR(structure);
      filename = 'POSCAR';
    } else if (format === 'xyz') {
      content = generateXYZ(structure);
    } else if (format === 'json') {
      content = generateJSON(structure);
    } else {
      continue;
    }

    exportedFiles[format] = content;

    if (output_directory) {
      const filepath = `${output_directory}/${filename}`;
      const writeResult = await writeFileSync(filepath, content);
      if (writeResult.success) {
        filePaths[format] = filepath;
      }
    }
  }

  return createSuccess({
    files: exportedFiles,
    file_paths: filePaths
  });
}

export async function handleExportStructure(args: unknown): Promise<any> {
  const result = await exportStructure(args);
  
  if (!result.success) {
    return {
      content: [{
        type: "text",
        text: `❌ **Export Failed**\n\n${result.error.message}`
      }],
      isError: true
    };
  }

  const data = result.data;
  let output = `## 📁 Structure Exported\n\n`;
  output += `**Formats:** ${Object.keys(data.files).join(', ')}\n\n`;
  
  if (Object.keys(data.file_paths).length > 0) {
    output += `### Files Written\n\n`;
    Object.entries(data.file_paths).forEach(([format, path]) => {
      output += `- **${format.toUpperCase()}:** \`${path}\`\n`;
    });
  }
  
  return {
    content: [{
      type: "text",
      text: output
    }]
  };
}
