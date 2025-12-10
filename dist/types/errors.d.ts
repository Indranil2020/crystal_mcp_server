/**
 * Error Type Definitions for Crystal MCP Server
 *
 * This module defines error codes, error structures, and result types
 * for defensive programming without try/catch blocks.
 */
/**
 * Error codes categorized by domain
 */
export declare enum CrystalErrorCode {
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
    GENERATION_TIMEOUT = "E2001",
    MAX_ATTEMPTS_EXCEEDED = "E2002",
    INCOMPATIBLE_COMPOSITION = "E2003",
    INVALID_WYCKOFF_ASSIGNMENT = "E2004",
    DISTANCE_CONSTRAINT_VIOLATION = "E2005",
    SYMMETRY_CONSTRAINT_VIOLATION = "E2006",
    GENERATION_FAILED = "E2007",
    ATOMS_TOO_CLOSE = "E3001",
    INVALID_SYMMETRY = "E3002",
    UNREALISTIC_DENSITY = "E3003",
    INVALID_CELL_PARAMETERS = "E3004",
    OVERLAPPING_ATOMS = "E3005",
    CONVERSION_FAILED = "E3006",
    SYMMETRY_DETECTION_FAILED = "E3007",
    REFINEMENT_FAILED = "E3008",
    MODEL_NOT_AVAILABLE = "E4001",
    OPTIMIZATION_FAILED = "E4002",
    ENERGY_CALCULATION_FAILED = "E4003",
    CONVERGENCE_FAILURE = "E4004",
    MODEL_LOAD_FAILED = "E4005",
    FILE_READ_ERROR = "E5001",
    FILE_WRITE_ERROR = "E5002",
    INVALID_FORMAT = "E5003",
    FILE_NOT_FOUND = "E5004",
    CORRUPTED_DATA = "E5005",
    PYTHON_EXECUTION_FAILED = "E6001",
    PYTHON_NOT_FOUND = "E6002",
    SCRIPT_NOT_FOUND = "E6003",
    INVALID_JSON_RESPONSE = "E6004",
    TIMEOUT_ERROR = "E6005",
    UNKNOWN_ERROR = "E9001",
    NOT_IMPLEMENTED = "E9002",
    INVALID_OPERATION = "E9003"
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
export declare function createError(code: CrystalErrorCode, message: string, details?: Record<string, unknown>, suggestions?: string[], recoverable?: boolean): CrystalError;
/**
 * Result type for operations that may fail
 * Uses discriminated union for type-safe error handling
 */
export type Result<T> = {
    success: true;
    data: T;
} | {
    success: false;
    error: CrystalError;
};
/**
 * Helper function to create a success result
 */
export declare function createSuccess<T>(data: T): Result<T>;
/**
 * Helper function to create a failure result
 */
export declare function createFailure<T>(error: CrystalError): Result<T>;
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
export declare function createValidationSuccess(): ValidationResult;
/**
 * Create a failed validation result
 */
export declare function createValidationFailure(errors: CrystalError[], warnings?: string[]): ValidationResult;
/**
 * Common error messages with actionable suggestions
 */
export declare const ERROR_MESSAGES: {
    readonly INVALID_SPACE_GROUP: (spg: number | string) => {
        message: string;
        suggestions: string[];
    };
    readonly INCOMPATIBLE_COMPOSITION: (composition: string[], spg: number) => {
        message: string;
        suggestions: string[];
    };
    readonly MAX_ATTEMPTS_EXCEEDED: (attempts: number) => {
        message: string;
        suggestions: string[];
    };
    readonly DISTANCE_CONSTRAINT_VIOLATION: (pair: string, distance: number, minDist: number) => {
        message: string;
        suggestions: string[];
    };
    readonly MODEL_NOT_AVAILABLE: (model: string) => {
        message: string;
        suggestions: string[];
    };
};
//# sourceMappingURL=errors.d.ts.map