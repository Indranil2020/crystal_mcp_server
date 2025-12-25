/**
 * Formatting Utilities
 *
 * Functions for formatting crystal structure data into human-readable output.
 */
import { CrystalStructure, StructureValidation } from "../types/crystal.js";
/**
 * Format crystal structure as human-readable text.
 */
export declare function formatStructureOutput(structure: CrystalStructure | any, validation?: StructureValidation): string;
/**
 * Format space group scan results.
 */
export declare function formatSpaceGroupScanOutput(results: any[]): string;
/**
 * Format optimization results.
 */
export declare function formatOptimizationOutput(result: any): string;
/**
 * Format ground state search results.
 */
export declare function formatGroundStateSearchOutput(result: any): string;
/**
 * Format symmetry analysis results.
 */
export declare function formatSymmetryOutput(result: any): string;
/**
 * Format validation results.
 */
export declare function formatValidationOutput(result: any): string;
/**
 * Format error message with suggestions.
 */
export declare function formatError(error: any): string;
//# sourceMappingURL=formatting.d.ts.map