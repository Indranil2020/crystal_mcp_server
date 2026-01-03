/**
 * Generate Crystal Tool
 *
 * Main tool for generating crystal structures with PyXtal.
 * Follows defensive programming with Result<T> pattern.
 */
import { StructureGenerationResult } from "../../types/crystal.js";
import { Result } from "../../types/errors.js";
/**
 * Generate a crystal structure with specified composition and space group.
 *
 * @param input - Generation parameters
 * @returns Result containing structure data or error
 */
export declare function generateCrystal(input: unknown): Promise<Result<StructureGenerationResult>>;
/**
 * MCP tool handler for generate_crystal.
 *
 * Returns two content items:
 * 1. Human-readable markdown for LLM consumption
 * 2. JSON data block for programmatic access (web GUI, CLI, etc.)
 */
export declare function handleGenerateCrystal(args: unknown): Promise<any>;
//# sourceMappingURL=generate-crystal.d.ts.map