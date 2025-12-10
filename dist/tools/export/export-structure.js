/**
 * Export Structure Tool
 *
 * Tool for exporting crystal structures to various file formats.
 */
import { ExportStructureSchema } from "../../types/tools.js";
import { createSuccess, createFailure, createError, CrystalErrorCode } from "../../types/errors.js";
import { writeFileSync } from "../../utils/file-io.js";
export async function exportStructure(input) {
    const parsed = ExportStructureSchema.safeParse(input);
    if (!parsed.success) {
        return createFailure(createError(CrystalErrorCode.INVALID_INPUT, "Invalid input parameters", { zodErrors: parsed.error.errors }, ["Check formats array and structure object"], true));
    }
    const { structure, formats, output_directory } = parsed.data;
    const exportedFiles = {};
    const filePaths = {};
    for (const format of formats) {
        let content = '';
        let filename = `structure.${format}`;
        if (format === 'cif') {
            content = generateCIF(structure);
        }
        else if (format === 'poscar') {
            content = generatePOSCAR(structure);
            filename = 'POSCAR';
        }
        else if (format === 'xyz') {
            content = generateXYZ(structure);
        }
        else if (format === 'json') {
            content = JSON.stringify(structure, null, 2);
        }
        else {
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
function generateCIF(structure) {
    const { lattice, space_group, atoms } = structure;
    let cif = `data_crystal\n`;
    cif += `_symmetry_space_group_name_H-M '${space_group.symbol}'\n`;
    cif += `_space_group_IT_number ${space_group.number}\n`;
    cif += `_cell_length_a ${lattice.a.toFixed(6)}\n`;
    cif += `_cell_length_b ${lattice.b.toFixed(6)}\n`;
    cif += `_cell_length_c ${lattice.c.toFixed(6)}\n`;
    cif += `_cell_angle_alpha ${lattice.alpha.toFixed(4)}\n`;
    cif += `_cell_angle_beta ${lattice.beta.toFixed(4)}\n`;
    cif += `_cell_angle_gamma ${lattice.gamma.toFixed(4)}\n`;
    cif += `_cell_volume ${lattice.volume.toFixed(4)}\n`;
    cif += `\nloop_\n`;
    cif += `_atom_site_label\n`;
    cif += `_atom_site_type_symbol\n`;
    cif += `_atom_site_fract_x\n`;
    cif += `_atom_site_fract_y\n`;
    cif += `_atom_site_fract_z\n`;
    atoms.forEach((atom, i) => {
        cif += `${atom.element}${i + 1} ${atom.element} ${atom.coords[0].toFixed(6)} ${atom.coords[1].toFixed(6)} ${atom.coords[2].toFixed(6)}\n`;
    });
    return cif;
}
function generatePOSCAR(structure) {
    const { lattice, atoms, metadata } = structure;
    let poscar = `${metadata.formula}\n`;
    poscar += `1.0\n`;
    lattice.matrix.forEach((vec) => {
        if (vec && vec[0] !== undefined && vec[1] !== undefined && vec[2] !== undefined) {
            poscar += `  ${vec[0].toFixed(12)}  ${vec[1].toFixed(12)}  ${vec[2].toFixed(12)}\n`;
        }
    });
    const elementCounts = {};
    const elementOrder = [];
    atoms.forEach((atom) => {
        const element = atom.element;
        if (element) {
            if (!elementCounts[element]) {
                elementCounts[element] = 0;
                elementOrder.push(element);
            }
            elementCounts[element]++;
        }
    });
    poscar += elementOrder.join(' ') + '\n';
    poscar += elementOrder.map(el => elementCounts[el]).join(' ') + '\n';
    poscar += `Direct\n`;
    elementOrder.forEach(element => {
        atoms.filter((a) => a.element === element).forEach((atom) => {
            poscar += `  ${atom.coords[0].toFixed(12)}  ${atom.coords[1].toFixed(12)}  ${atom.coords[2].toFixed(12)}\n`;
        });
    });
    return poscar;
}
function generateXYZ(structure) {
    const { atoms, metadata } = structure;
    let xyz = `${atoms.length}\n`;
    xyz += `${metadata.formula}\n`;
    atoms.forEach((atom) => {
        xyz += `${atom.element}  ${atom.cartesian[0].toFixed(8)}  ${atom.cartesian[1].toFixed(8)}  ${atom.cartesian[2].toFixed(8)}\n`;
    });
    return xyz;
}
export async function handleExportStructure(args) {
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
//# sourceMappingURL=export-structure.js.map