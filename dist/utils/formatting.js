/**
 * Formatting Utilities
 *
 * Functions for formatting crystal structure data into human-readable output.
 */
/**
 * Format crystal structure as human-readable text.
 */
export function formatStructureOutput(structure, validation) {
    // Handle missing fields gracefully
    const lattice = structure?.lattice || {};
    const space_group = structure?.space_group || {};
    const atoms = structure?.atoms || structure?.sites || [];
    const metadata = structure?.metadata || {};
    let output = `## Generated Crystal Structure\n\n`;
    // Basic information - only show if available
    if (metadata.formula) {
        output += `**Formula:** ${metadata.formula}\n`;
    }
    if (space_group.number !== undefined || space_group.symbol) {
        const sgInfo = space_group.number !== undefined
            ? `${space_group.number}${space_group.symbol ? ` (${space_group.symbol})` : ''}`
            : space_group.symbol || 'Unknown';
        output += `**Space Group:** ${sgInfo}\n`;
    }
    if (space_group.crystal_system) {
        output += `**Crystal System:** ${space_group.crystal_system}\n`;
    }
    if (space_group.point_group) {
        output += `**Point Group:** ${space_group.point_group}\n`;
    }
    if (space_group.hall_symbol) {
        output += `**Hall Symbol:** ${space_group.hall_symbol}\n`;
    }
    output += '\n';
    // Lattice parameters - only show if available
    if (lattice.a !== undefined) {
        output += `### Lattice Parameters\n\n`;
        output += `| Parameter | Value |\n`;
        output += `|-----------|-------|\n`;
        output += `| a | ${(lattice.a ?? 0).toFixed(4)} A |\n`;
        output += `| b | ${(lattice.b ?? 0).toFixed(4)} A |\n`;
        output += `| c | ${(lattice.c ?? 0).toFixed(4)} A |\n`;
        output += `| alpha | ${(lattice.alpha ?? 0).toFixed(2)} deg |\n`;
        output += `| beta | ${(lattice.beta ?? 0).toFixed(2)} deg |\n`;
        output += `| gamma | ${(lattice.gamma ?? 0).toFixed(2)} deg |\n`;
        if (lattice.volume !== undefined) {
            output += `| Volume | ${lattice.volume.toFixed(3)} A^3 |\n`;
        }
        output += `\n`;
    }
    // Physical properties - only show if available
    const hasMetadata = metadata.natoms !== undefined || metadata.density !== undefined || metadata.n_atoms !== undefined;
    if (hasMetadata) {
        output += `### Physical Properties\n\n`;
        output += `| Property | Value |\n`;
        output += `|----------|-------|\n`;
        const atomCount = metadata.natoms ?? metadata.n_atoms ?? atoms.length;
        if (atomCount !== undefined) {
            output += `| Number of Atoms | ${atomCount} |\n`;
        }
        if (metadata.density !== undefined) {
            output += `| Density | ${metadata.density.toFixed(3)} g/cm^3 |\n`;
        }
        if (metadata.packing_fraction !== undefined) {
            output += `| Packing Fraction | ${metadata.packing_fraction.toFixed(3)} |\n`;
        }
        output += `\n`;
    }
    // Validation results
    if (validation) {
        if (!validation.valid || validation.issues.length > 0) {
            output += `### Validation Issues\n\n`;
            validation.issues.forEach(issue => {
                output += `- **${issue}**\n`;
            });
            output += `\n`;
        }
        if (validation.warnings.length > 0) {
            output += `### Warnings\n\n`;
            validation.warnings.forEach(warning => {
                output += `- ${warning}\n`;
            });
            output += `\n`;
        }
    }
    // Atomic positions - only show if atoms are available
    if (atoms && atoms.length > 0) {
        const hasCartesian = atoms.some((a) => a.cartesian);
        output += `### Atomic Positions\n\n`;
        if (hasCartesian) {
            output += `| # | Element | Wyckoff | Fractional Coordinates | Cartesian (A) |\n`;
            output += `|---|---------|---------|------------------------|---------------|\n`;
        }
        else {
            output += `| # | Element | Wyckoff | Fractional Coordinates |\n`;
            output += `|---|---------|---------|------------------------|\n`;
        }
        atoms.forEach((atom, i) => {
            const coords = atom.coords || atom.abc || [0, 0, 0];
            const frac = coords.map((x) => (x ?? 0).toFixed(4)).join(', ');
            const element = atom.element || atom.species?.[0]?.element || 'X';
            const wyckoff = atom.wyckoff || '-';
            if (hasCartesian) {
                const cart = atom.cartesian
                    ? atom.cartesian.map((x) => (x ?? 0).toFixed(3)).join(', ')
                    : '-';
                output += `| ${i + 1} | ${element} | ${wyckoff} | (${frac}) | (${cart}) |\n`;
            }
            else {
                output += `| ${i + 1} | ${element} | ${wyckoff} | (${frac}) |\n`;
            }
        });
        output += `\n`;
    }
    // Lattice vectors - only show if matrix is available
    if (lattice.matrix && Array.isArray(lattice.matrix)) {
        output += `### Lattice Vectors\n\n`;
        output += `\`\`\`\n`;
        lattice.matrix.forEach((vec, i) => {
            const label = ['a', 'b', 'c'][i];
            output += `${label} = [${vec.map((x) => (x ?? 0).toFixed(6)).join(', ')}]\n`;
        });
        output += `\`\`\`\n\n`;
    }
    output += `Structure generated successfully!\n`;
    return output;
}
/**
 * Format space group scan results.
 */
export function formatSpaceGroupScanOutput(results) {
    let output = `## Space Group Scan Results\n\n`;
    const successful = results.filter(r => r.success);
    const failed = results.filter(r => !r.success);
    output += `**Total:** ${results.length} space groups\n`;
    output += `**Successful:** ${successful.length}\n`;
    output += `**Failed:** ${failed.length}\n\n`;
    if (successful.length > 0) {
        output += `### Successfully Generated\n\n`;
        output += `| Space Group | Symbol | Crystal System | Atoms | Volume (A^3) |\n`;
        output += `|-------------|--------|----------------|-------|------------|\n`;
        successful.forEach(result => {
            // Handle both nested (structure.space_group) and flat (spacegroup_symbol) formats
            const spgNumber = result.space_group || result.structure?.space_group?.number || '?';
            const spgSymbol = result.spacegroup_symbol || result.structure?.space_group?.symbol || '?';
            const crystalSystem = result.crystal_system || result.structure?.space_group?.crystal_system || '?';
            const nAtoms = result.structure?.metadata?.natoms || result.structure?.metadata?.n_atoms || result.structure?.atoms?.length || '?';
            const volume = result.structure?.lattice?.volume;
            const volumeStr = volume !== undefined ? volume.toFixed(2) : '?';
            output += `| ${spgNumber} | ${spgSymbol} | ${crystalSystem} | ${nAtoms} | ${volumeStr} |\n`;
        });
        output += `\n`;
    }
    if (failed.length > 0) {
        output += `### Failed to Generate\n\n`;
        output += `| Space Group | Error |\n`;
        output += `|-------------|-------|\n`;
        failed.forEach(result => {
            const sg = result.space_group || '?';
            const errMsg = result.error_message || result.error || 'Unknown error';
            output += `| ${sg} | ${errMsg} |\n`;
        });
        output += `\n`;
    }
    return output;
}
/**
 * Format optimization results.
 */
export function formatOptimizationOutput(result) {
    let output = `## MLFF Optimization Results\n\n`;
    output += `### Energy\n\n`;
    output += `| Property | Value |\n`;
    output += `|----------|-------|\n`;
    output += `| Initial Energy | ${result.initial_energy.toFixed(4)} eV |\n`;
    output += `| Final Energy | ${result.final_energy.toFixed(4)} eV |\n`;
    output += `| Energy Change | ${result.energy_change.toFixed(4)} eV |\n`;
    output += `| Energy per Atom | ${(result.final_energy / result.optimized_structure.metadata.natoms).toFixed(4)} eV |\n\n`;
    output += `### Forces\n\n`;
    output += `| Property | Value |\n`;
    output += `|----------|-------|\n`;
    output += `| Initial Max Force | ${result.max_force_initial.toFixed(4)} eV/Å |\n`;
    output += `| Final Max Force | ${result.max_force_final.toFixed(4)} eV/Å |\n\n`;
    output += `### Convergence\n\n`;
    output += `| Property | Value |\n`;
    output += `|----------|-------|\n`;
    output += `| Steps | ${result.n_steps} |\n`;
    output += `| Converged | ${result.converged ? 'Yes' : 'No'} |\n`;
    if (result.preserved_symmetry !== undefined) {
        output += `| Symmetry Preserved | ${result.preserved_symmetry ? 'Yes' : 'No'} |\n`;
    }
    output += `\n`;
    if (result.timing) {
        output += `### Timing\n\n`;
        output += `| Phase | Time |\n`;
        output += `|-------|------|\n`;
        output += `| Optimization | ${(result.timing.optimization_time_ms / 1000).toFixed(2)} s |\n`;
        output += `| Total | ${(result.timing.total_time_ms / 1000).toFixed(2)} s |\n\n`;
    }
    return output;
}
/**
 * Format ground state search results.
 */
export function formatGroundStateSearchOutput(result) {
    let output = `## 🏆 Ground State Search Results\n\n`;
    const gs = result.ground_state;
    output += `### Ground State Structure\n\n`;
    output += `**Space Group:** ${gs.space_group} (${gs.structure.space_group.symbol})\n`;
    output += `**Energy:** ${gs.energy.toFixed(4)} eV\n`;
    output += `**Energy per Atom:** ${(gs.energy / gs.structure.metadata.natoms).toFixed(4)} eV\n`;
    output += `**Formula:** ${gs.structure.metadata.formula}\n\n`;
    output += `### Energy Ranking (Top 10)\n\n`;
    output += `| Rank | Space Group | Symbol | Energy (eV) | ΔE (meV/atom) |\n`;
    output += `|------|-------------|--------|-------------|---------------|\n`;
    const topTen = result.energy_ranking.slice(0, 10);
    topTen.forEach((entry) => {
        const deltaE = entry.energy_above_ground_state * 1000; // Convert to meV
        output += `| ${entry.rank} | ${entry.space_group} | ${entry.symbol || '-'} | ${entry.energy.toFixed(4)} | ${deltaE.toFixed(2)} |\n`;
    });
    output += `\n`;
    output += `### Statistics\n\n`;
    output += `| Property | Value |\n`;
    output += `|----------|-------|\n`;
    output += `| Total Structures | ${result.statistics.total_structures} |\n`;
    output += `| Converged | ${result.statistics.converged_structures} |\n`;
    output += `| Unique Space Groups | ${result.statistics.unique_space_groups} |\n`;
    output += `| Energy Range | ${result.statistics.energy_range[0].toFixed(4)} to ${result.statistics.energy_range[1].toFixed(4)} eV |\n`;
    output += `| Most Stable System | ${result.statistics.most_stable_crystal_system} |\n\n`;
    return output;
}
/**
 * Format symmetry analysis results.
 */
export function formatSymmetryOutput(result) {
    let output = `## 🔷 Symmetry Analysis Results\n\n`;
    const spg = result.space_group;
    const symops = result.symmetry_operations;
    output += `### Space Group\n\n`;
    output += `| Property | Value |\n`;
    output += `|----------|-------|\n`;
    output += `| Number | ${spg.number} |\n`;
    output += `| Symbol | ${spg.symbol} |\n`;
    output += `| Hall Symbol | ${spg.hall_symbol} |\n`;
    output += `| Point Group | ${spg.point_group} |\n`;
    output += `| Crystal System | ${spg.crystal_system} |\n\n`;
    output += `### Symmetry Operations\n\n`;
    output += `**Total Operations:** ${symops.n_operations}\n\n`;
    if (result.wyckoff_positions && result.wyckoff_positions.length > 0) {
        output += `### Wyckoff Positions\n\n`;
        output += `| Element | Wyckoff | Multiplicity | Coordinates |\n`;
        output += `|---------|---------|--------------|-------------|\n`;
        result.wyckoff_positions.forEach((wp) => {
            const coords = wp.coords.map((x) => x.toFixed(4)).join(', ');
            output += `| ${wp.element} | ${wp.wyckoff} | ${wp.multiplicity} | (${coords}) |\n`;
        });
        output += `\n`;
    }
    if (result.primitive_cell) {
        output += `### Primitive Cell Found\n\n`;
        output += `**Atoms:** ${result.primitive_cell.metadata.natoms}\n\n`;
    }
    if (result.conventional_cell) {
        output += `### Conventional Cell\n\n`;
        output += `**Atoms:** ${result.conventional_cell.metadata.natoms}\n`;
        output += `**Standardized:** ${result.is_standardized ? 'Yes' : 'No'}\n\n`;
    }
    return output;
}
/**
 * Format validation results.
 */
export function formatValidationOutput(result) {
    let output = `## Structure Validation Results\n\n`;
    output += `**Status:** ${result.valid ? '✅ Valid' : '❌ Invalid'}\n\n`;
    if (result.errors && result.errors.length > 0) {
        output += `### ❌ Errors\n\n`;
        result.errors.forEach((error) => {
            output += `- **[${error.severity.toUpperCase()}]** ${error.message}\n`;
            if (error.details && Object.keys(error.details).length > 0) {
                output += `  - Details: ${JSON.stringify(error.details)}\n`;
            }
        });
        output += `\n`;
    }
    if (result.warnings && result.warnings.length > 0) {
        output += `### ⚠️ Warnings\n\n`;
        result.warnings.forEach((warning) => {
            output += `- ${warning.message || warning}\n`;
        });
        output += `\n`;
    }
    if (result.metrics && Object.keys(result.metrics).length > 0) {
        output += `### 📊 Metrics\n\n`;
        output += `| Metric | Value |\n`;
        output += `|--------|-------|\n`;
        for (const [key, value] of Object.entries(result.metrics)) {
            const formattedValue = typeof value === 'number' ? value.toFixed(4) : String(value);
            const formattedKey = key.replace(/_/g, ' ').replace(/\b\w/g, l => l.toUpperCase());
            output += `| ${formattedKey} | ${formattedValue} |\n`;
        }
        output += `\n`;
    }
    if (result.suggestions && result.suggestions.length > 0) {
        output += `### 💡 Suggestions\n\n`;
        result.suggestions.forEach((suggestion) => {
            output += `- ${suggestion}\n`;
        });
        output += `\n`;
    }
    return output;
}
/**
 * Format error message with suggestions.
 */
export function formatError(error) {
    let output = `## ❌ Error\n\n`;
    output += `**Code:** \`${error.code}\`\n`;
    output += `**Message:** ${error.message}\n\n`;
    if (error.details && Object.keys(error.details).length > 0) {
        output += `### Details\n\n`;
        output += `\`\`\`json\n${JSON.stringify(error.details, null, 2)}\n\`\`\`\n\n`;
    }
    if (error.suggestions && error.suggestions.length > 0) {
        output += `### 💡 Suggestions\n\n`;
        error.suggestions.forEach((suggestion) => {
            output += `- ${suggestion}\n`;
        });
        output += `\n`;
    }
    return output;
}
//# sourceMappingURL=formatting.js.map