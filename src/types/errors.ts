/**
 * Error Type Definitions for Crystal MCP Server
 * 
 * This module defines error codes, error structures, and result types
 * for defensive programming without try/catch blocks.
 */

/**
 * Error codes categorized by domain
 */
export enum CrystalErrorCode {
  // Input validation errors (E1xxx)
  INVALID_INPUT = "E1000",
  INVALID_SPACE_GROUP = "E1001",
  INVALID_COMPOSITION = "E1002",
  INVALID_LATTICE_PARAMS = "E1003",
  INVALID_WYCKOFF_POSITION = "E1004",
  INVALID_VOLUME_FACTOR = "E1005",
  INVALID_MIN_DISTANCE = "E1006",
  INVALID_NUM_ATOMS = "E1007",
  INVALID_PARAMETER = "E1008",
  INVALID_OPTIMIZER = "E1009",
  INVALID_STRUCTURE = "E1010",
  INVALID_USAGE = "E1011",

  // Generation failures (E2xxx)
  GENERATION_TIMEOUT = "E2001",
  MAX_ATTEMPTS_EXCEEDED = "E2002",
  INCOMPATIBLE_COMPOSITION = "E2003",
  INVALID_WYCKOFF_ASSIGNMENT = "E2004",
  DISTANCE_CONSTRAINT_VIOLATION = "E2005",
  SYMMETRY_CONSTRAINT_VIOLATION = "E2006",
  GENERATION_FAILED = "E2007",

  // Structure validation errors (E3xxx)
  ATOMS_TOO_CLOSE = "E3001",
  INVALID_SYMMETRY = "E3002",
  UNREALISTIC_DENSITY = "E3003",
  INVALID_CELL_PARAMETERS = "E3004",
  OVERLAPPING_ATOMS = "E3005",
  CONVERSION_FAILED = "E3006",
  SYMMETRY_DETECTION_FAILED = "E3007",
  REFINEMENT_FAILED = "E3008",

  // MLFF errors (E4xxx)
  MODEL_NOT_AVAILABLE = "E4001",
  OPTIMIZATION_FAILED = "E4002",
  ENERGY_CALCULATION_FAILED = "E4003",
  CONVERGENCE_FAILURE = "E4004",
  MODEL_LOAD_FAILED = "E4005",

  // File I/O errors (E5xxx)
  FILE_READ_ERROR = "E5001",
  FILE_WRITE_ERROR = "E5002",
  INVALID_FORMAT = "E5003",
  FILE_NOT_FOUND = "E5004",
  CORRUPTED_DATA = "E5005",

  // Python bridge errors (E6xxx)
  PYTHON_EXECUTION_FAILED = "E6001",
  PYTHON_NOT_FOUND = "E6002",
  SCRIPT_NOT_FOUND = "E6003",
  INVALID_JSON_RESPONSE = "E6004",
  TIMEOUT_ERROR = "E6005",

  // General errors (E9xxx)
  EXECUTION_ERROR = "E9000",
  UNKNOWN_ERROR = "E9001",
  NOT_IMPLEMENTED = "E9002",
  INVALID_OPERATION = "E9003",
}

/**
 * Detailed error information structure
 */
export interface CrystalError {
  readonly code: CrystalErrorCode;
  readonly message: string;
  readonly details: Record<string, unknown>;
  readonly suggestions: readonly string[];
  readonly recoverable: boolean;
  readonly timestamp: string;
}

/**
 * Factory function to create CrystalError instances
 */
export function createError(
  code: CrystalErrorCode,
  message: string,
  details: Record<string, unknown> = {},
  suggestions: string[] = [],
  recoverable: boolean = true
): CrystalError {
  return {
    code,
    message,
    details,
    suggestions,
    recoverable,
    timestamp: new Date().toISOString()
  };
}

/**
 * Result type for operations that may fail
 * Uses discriminated union for type-safe error handling
 */
export type Result<T> =
  | { success: true; data: T }
  | { success: false; error: CrystalError };

/**
 * Helper function to create a success result
 */
export function createSuccess<T>(data: T): Result<T> {
  return { success: true, data };
}

/**
 * Helper function to create a failure result
 */
export function createFailure<T>(error: CrystalError): Result<T> {
  return { success: false, error };
}

/**
 * Optional type for operations that may return no value
 */
export type Optional<T> = T | null;

/**
 * Validation result for input validation
 */
export interface ValidationResult {
  readonly valid: boolean;
  readonly errors: readonly CrystalError[];
  readonly warnings: readonly string[];
}

/**
 * Create a successful validation result
 */
export function createValidationSuccess(): ValidationResult {
  return {
    valid: true,
    errors: [],
    warnings: []
  };
}

/**
 * Create a failed validation result
 */
export function createValidationFailure(
  errors: CrystalError[],
  warnings: string[] = []
): ValidationResult {
  return {
    valid: false,
    errors,
    warnings
  };
}

/**
 * Common error messages with actionable suggestions
 */
export const ERROR_MESSAGES = {
  INVALID_SPACE_GROUP: (spg: number | string) => ({
    message: `Invalid space group: ${spg}. Must be between 1 and 230.`,
    suggestions: [
      "Check the space group number is correct",
      "Refer to International Tables for Crystallography",
      "Use Hermann-Mauguin notation if providing a symbol"
    ]
  }),

  INCOMPATIBLE_COMPOSITION: (composition: string[], spg: number) => ({
    message: `Composition ${composition.join(", ")} is incompatible with space group ${spg}.`,
    suggestions: [
      "Adjust the number of atoms to match Wyckoff position multiplicities",
      "Try a different space group with compatible symmetry",
      "Consult the International Tables for compatible compositions",
      "Increase volume_factor to 1.2-1.5 to allow more space"
    ]
  }),

  MAX_ATTEMPTS_EXCEEDED: (attempts: number) => ({
    message: `Structure generation failed after ${attempts} attempts.`,
    suggestions: [
      "Increase volume_factor (try 1.2-1.5)",
      "Relax minimum distance constraints",
      "Try a different space group with similar symmetry",
      "Reduce the number of atoms if specified",
      "Check that composition is compatible with space group"
    ]
  }),

  DISTANCE_CONSTRAINT_VIOLATION: (pair: string, distance: number, minDist: number) => ({
    message: `Distance between ${pair} is ${distance.toFixed(3)} Å, less than minimum ${minDist.toFixed(3)} Å.`,
    suggestions: [
      `Reduce min_distance for ${pair} to ${(distance * 0.9).toFixed(3)} Å`,
      "Increase volume_factor to allow more space",
      "Try a different space group with lower symmetry"
    ]
  }),

  MODEL_NOT_AVAILABLE: (model: string) => ({
    message: `MLFF model '${model}' is not available or failed to load.`,
    suggestions: [
      `Install the required package: pip install ${model}`,
      "Check that the model is properly installed",
      "Try a different MLFF model (chgnet, m3gnet, mace)",
      "Verify Python environment has required dependencies"
    ]
  })
} as const;
