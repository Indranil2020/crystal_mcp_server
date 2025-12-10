/**
 * Formatting Utilities
 * 
 * Functions for formatting crystal structure data into human-readable output.
 */

import { CrystalStructure, StructureValidation } from "../types/crystal.js";

/**
 * Format crystal structure as human-readable text.
 */
export function formatStructureOutput(
  structure: CrystalStructure,
  validation?: StructureValidation
): string {
  const { lattice, space_group, atoms, metadata } = structure;
  
  let output = `## 🔬 Generated Crystal Structure\n\n`;
  
  // Basic information
  output += `**Formula:** ${metadata.formula}\n`;
  output += `**Space Group:** ${space_group.number} (${space_group.symbol})\n`;
  output += `**Crystal System:** ${space_group.crystal_system}\n`;
  output += `**Point Group:** ${space_group.point_group}\n`;
  output += `**Hall Symbol:** ${space_group.hall_symbol}\n\n`;
  
  // Lattice parameters
  output += `### Lattice Parameters\n\n`;
  output += `| Parameter | Value |\n`;
  output += `|-----------|-------|\n`;
  output += `| a | ${lattice.a.toFixed(4)} Å |\n`;
  output += `| b | ${lattice.b.toFixed(4)} Å |\n`;
  output += `| c | ${lattice.c.toFixed(4)} Å |\n`;
  output += `| α | ${lattice.alpha.toFixed(2)}° |\n`;
  output += `| β | ${lattice.beta.toFixed(2)}° |\n`;
  output += `| γ | ${lattice.gamma.toFixed(2)}° |\n`;
  output += `| Volume | ${lattice.volume.toFixed(3)} ų |\n\n`;
  
  // Physical properties
  output += `### Physical Properties\n\n`;
  output += `| Property | Value |\n`;
  output += `|----------|-------|\n`;
  output += `| Number of Atoms | ${metadata.natoms} |\n`;
  output += `| Density | ${metadata.density.toFixed(3)} g/cm³ |\n`;
  if (metadata.packing_fraction !== undefined) {
    output += `| Packing Fraction | ${metadata.packing_fraction.toFixed(3)} |\n`;
  }
  output += `\n`;
  
  // Validation results
  if (validation) {
    if (!validation.valid || validation.issues.length > 0) {
      output += `### ⚠️ Validation Issues\n\n`;
      validation.issues.forEach(issue => {
        output += `- **${issue}**\n`;
      });
      output += `\n`;
    }
    
    if (validation.warnings.length > 0) {
      output += `### ℹ️ Warnings\n\n`;
      validation.warnings.forEach(warning => {
        output += `- ${warning}\n`;
      });
      output += `\n`;
    }
  }
  
  // Atomic positions
  output += `### Atomic Positions\n\n`;
  output += `| # | Element | Wyckoff | Fractional Coordinates | Cartesian (Å) |\n`;
  output += `|---|---------|---------|------------------------|---------------|\n`;
  
  atoms.forEach((atom, i) => {
    const frac = atom.coords.map(x => x.toFixed(4)).join(', ');
    const cart = atom.cartesian.map(x => x.toFixed(3)).join(', ');
    const wyckoff = atom.wyckoff || '-';
    output += `| ${i+1} | ${atom.element} | ${wyckoff} | (${frac}) | (${cart}) |\n`;
  });
  
  output += `\n`;
  
  // Lattice vectors
  output += `### Lattice Vectors\n\n`;
  output += `\`\`\`\n`;
  lattice.matrix.forEach((vec, i) => {
    const label = ['a', 'b', 'c'][i];
    output += `${label} = [${vec.map(x => x.toFixed(6)).join(', ')}]\n`;
  });
  output += `\`\`\`\n\n`;
  
  output += `✅ Structure generated successfully!\n`;
  
  return output;
}

/**
 * Format space group scan results.
 */
export function formatSpaceGroupScanOutput(results: any[]): string {
  let output = `## 🔍 Space Group Scan Results\n\n`;
  
  const successful = results.filter(r => r.success);
  const failed = results.filter(r => !r.success);
  
  output += `**Total:** ${results.length} space groups\n`;
  output += `**Successful:** ${successful.length}\n`;
  output += `**Failed:** ${failed.length}\n\n`;
  
  if (successful.length > 0) {
    output += `### ✅ Successfully Generated\n\n`;
    output += `| Space Group | Symbol | Crystal System | Atoms | Volume (ų) |\n`;
    output += `|-------------|--------|----------------|-------|------------|\n`;
    
    successful.forEach(result => {
      const spg = result.structure.space_group;
      const meta = result.structure.metadata;
      output += `| ${spg.number} | ${spg.symbol} | ${spg.crystal_system} | ${meta.natoms} | ${result.structure.lattice.volume.toFixed(2)} |\n`;
    });
    
    output += `\n`;
  }
  
  if (failed.length > 0) {
    output += `### ❌ Failed to Generate\n\n`;
    output += `| Space Group | Error |\n`;
    output += `|-------------|-------|\n`;
    
    failed.forEach(result => {
      output += `| ${result.space_group} | ${result.error_message} |\n`;
    });
    
    output += `\n`;
  }
  
  return output;
}

/**
 * Format optimization results.
 */
export function formatOptimizationOutput(result: any): string {
  let output = `## ⚡ MLFF Optimization Results\n\n`;
  
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
  output += `| Converged | ${result.converged ? 'Yes ✅' : 'No ❌'} |\n`;
  if (result.preserved_symmetry !== undefined) {
    output += `| Symmetry Preserved | ${result.preserved_symmetry ? 'Yes ✅' : 'No ❌'} |\n`;
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
export function formatGroundStateSearchOutput(result: any): string {
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
  topTen.forEach((entry: any) => {
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
export function formatSymmetryOutput(result: any): string {
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
    
    result.wyckoff_positions.forEach((wp: any) => {
      const coords = wp.coords.map((x: number) => x.toFixed(4)).join(', ');
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
    output += `**Standardized:** ${result.is_standardized ? 'Yes ✅' : 'No'}\n\n`;
  }
  
  return output;
}

/**
 * Format validation results.
 */
export function formatValidationOutput(result: any): string {
  let output = `## ✓ Structure Validation Results\n\n`;
  
  output += `**Status:** ${result.valid ? '✅ Valid' : '❌ Invalid'}\n\n`;
  
  if (result.errors && result.errors.length > 0) {
    output += `### ❌ Errors\n\n`;
    result.errors.forEach((error: any) => {
      output += `- **[${error.severity.toUpperCase()}]** ${error.message}\n`;
      if (error.details && Object.keys(error.details).length > 0) {
        output += `  - Details: ${JSON.stringify(error.details)}\n`;
      }
    });
    output += `\n`;
  }
  
  if (result.warnings && result.warnings.length > 0) {
    output += `### ⚠️ Warnings\n\n`;
    result.warnings.forEach((warning: any) => {
      output += `- ${warning.message || warning}\n`;
    });
    output += `\n`;
  }
  
  if (result.metrics && Object.keys(result.metrics).length > 0) {
    output += `### 📊 Metrics\n\n`;
    output += `| Metric | Value |\n`;
    output += `|--------|-------|\n`;
    
    for (const [key, value] of Object.entries(result.metrics)) {
      const formattedValue = typeof value === 'number' ? (value as number).toFixed(4) : String(value);
      const formattedKey = key.replace(/_/g, ' ').replace(/\b\w/g, l => l.toUpperCase());
      output += `| ${formattedKey} | ${formattedValue} |\n`;
    }
    
    output += `\n`;
  }
  
  if (result.suggestions && result.suggestions.length > 0) {
    output += `### 💡 Suggestions\n\n`;
    result.suggestions.forEach((suggestion: string) => {
      output += `- ${suggestion}\n`;
    });
    output += `\n`;
  }
  
  return output;
}

/**
 * Format error message with suggestions.
 */
export function formatError(error: any): string {
  let output = `## ❌ Error\n\n`;
  
  output += `**Code:** \`${error.code}\`\n`;
  output += `**Message:** ${error.message}\n\n`;
  
  if (error.details && Object.keys(error.details).length > 0) {
    output += `### Details\n\n`;
    output += `\`\`\`json\n${JSON.stringify(error.details, null, 2)}\n\`\`\`\n\n`;
  }
  
  if (error.suggestions && error.suggestions.length > 0) {
    output += `### 💡 Suggestions\n\n`;
    error.suggestions.forEach((suggestion: string) => {
      output += `- ${suggestion}\n`;
    });
    output += `\n`;
  }
  
  return output;
}
