/**
 * Error Type Definitions for Crystal MCP Server
 *
 * This module defines error codes, error structures, and result types
 * for defensive programming without try/catch blocks.
 */
/**
 * Error codes categorized by domain
 */
export var CrystalErrorCode;
(function (CrystalErrorCode) {
    // Input validation errors (E1xxx)
    CrystalErrorCode["INVALID_INPUT"] = "E1000";
    CrystalErrorCode["INVALID_SPACE_GROUP"] = "E1001";
    CrystalErrorCode["INVALID_COMPOSITION"] = "E1002";
    CrystalErrorCode["INVALID_LATTICE_PARAMS"] = "E1003";
    CrystalErrorCode["INVALID_WYCKOFF_POSITION"] = "E1004";
    CrystalErrorCode["INVALID_VOLUME_FACTOR"] = "E1005";
    CrystalErrorCode["INVALID_MIN_DISTANCE"] = "E1006";
    CrystalErrorCode["INVALID_NUM_ATOMS"] = "E1007";
    CrystalErrorCode["INVALID_PARAMETER"] = "E1008";
    CrystalErrorCode["INVALID_OPTIMIZER"] = "E1009";
    CrystalErrorCode["INVALID_STRUCTURE"] = "E1010";
    CrystalErrorCode["INVALID_USAGE"] = "E1011";
    // Generation failures (E2xxx)
    CrystalErrorCode["GENERATION_TIMEOUT"] = "E2001";
    CrystalErrorCode["MAX_ATTEMPTS_EXCEEDED"] = "E2002";
    CrystalErrorCode["INCOMPATIBLE_COMPOSITION"] = "E2003";
    CrystalErrorCode["INVALID_WYCKOFF_ASSIGNMENT"] = "E2004";
    CrystalErrorCode["DISTANCE_CONSTRAINT_VIOLATION"] = "E2005";
    CrystalErrorCode["SYMMETRY_CONSTRAINT_VIOLATION"] = "E2006";
    CrystalErrorCode["GENERATION_FAILED"] = "E2007";
    // Structure validation errors (E3xxx)
    CrystalErrorCode["ATOMS_TOO_CLOSE"] = "E3001";
    CrystalErrorCode["INVALID_SYMMETRY"] = "E3002";
    CrystalErrorCode["UNREALISTIC_DENSITY"] = "E3003";
    CrystalErrorCode["INVALID_CELL_PARAMETERS"] = "E3004";
    CrystalErrorCode["OVERLAPPING_ATOMS"] = "E3005";
    CrystalErrorCode["CONVERSION_FAILED"] = "E3006";
    CrystalErrorCode["SYMMETRY_DETECTION_FAILED"] = "E3007";
    CrystalErrorCode["REFINEMENT_FAILED"] = "E3008";
    // MLFF errors (E4xxx)
    CrystalErrorCode["MODEL_NOT_AVAILABLE"] = "E4001";
    CrystalErrorCode["OPTIMIZATION_FAILED"] = "E4002";
    CrystalErrorCode["ENERGY_CALCULATION_FAILED"] = "E4003";
    CrystalErrorCode["CONVERGENCE_FAILURE"] = "E4004";
    CrystalErrorCode["MODEL_LOAD_FAILED"] = "E4005";
    // File I/O errors (E5xxx)
    CrystalErrorCode["FILE_READ_ERROR"] = "E5001";
    CrystalErrorCode["FILE_WRITE_ERROR"] = "E5002";
    CrystalErrorCode["INVALID_FORMAT"] = "E5003";
    CrystalErrorCode["FILE_NOT_FOUND"] = "E5004";
    CrystalErrorCode["CORRUPTED_DATA"] = "E5005";
    // Python bridge errors (E6xxx)
    CrystalErrorCode["PYTHON_EXECUTION_FAILED"] = "E6001";
    CrystalErrorCode["PYTHON_NOT_FOUND"] = "E6002";
    CrystalErrorCode["SCRIPT_NOT_FOUND"] = "E6003";
    CrystalErrorCode["INVALID_JSON_RESPONSE"] = "E6004";
    CrystalErrorCode["TIMEOUT_ERROR"] = "E6005";
    // General errors (E9xxx)
    CrystalErrorCode["EXECUTION_ERROR"] = "E9000";
    CrystalErrorCode["UNKNOWN_ERROR"] = "E9001";
    CrystalErrorCode["NOT_IMPLEMENTED"] = "E9002";
    CrystalErrorCode["INVALID_OPERATION"] = "E9003";
})(CrystalErrorCode || (CrystalErrorCode = {}));
/**
 * Factory function to create CrystalError instances
 */
export function createError(code, message, details = {}, suggestions = [], recoverable = true) {
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
 * Helper function to create a success result
 */
export function createSuccess(data) {
    return { success: true, data };
}
/**
 * Helper function to create a failure result
 */
export function createFailure(error) {
    return { success: false, error };
}
/**
 * Create a successful validation result
 */
export function createValidationSuccess() {
    return {
        valid: true,
        errors: [],
        warnings: []
    };
}
/**
 * Create a failed validation result
 */
export function createValidationFailure(errors, warnings = []) {
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
    INVALID_SPACE_GROUP: (spg) => ({
        message: `Invalid space group: ${spg}. Must be between 1 and 230.`,
        suggestions: [
            "Check the space group number is correct",
            "Refer to International Tables for Crystallography",
            "Use Hermann-Mauguin notation if providing a symbol"
        ]
    }),
    INCOMPATIBLE_COMPOSITION: (composition, spg) => ({
        message: `Composition ${composition.join(", ")} is incompatible with space group ${spg}.`,
        suggestions: [
            "Adjust the number of atoms to match Wyckoff position multiplicities",
            "Try a different space group with compatible symmetry",
            "Consult the International Tables for compatible compositions",
            "Increase volume_factor to 1.2-1.5 to allow more space"
        ]
    }),
    MAX_ATTEMPTS_EXCEEDED: (attempts) => ({
        message: `Structure generation failed after ${attempts} attempts.`,
        suggestions: [
            "Increase volume_factor (try 1.2-1.5)",
            "Relax minimum distance constraints",
            "Try a different space group with similar symmetry",
            "Reduce the number of atoms if specified",
            "Check that composition is compatible with space group"
        ]
    }),
    DISTANCE_CONSTRAINT_VIOLATION: (pair, distance, minDist) => ({
        message: `Distance between ${pair} is ${distance.toFixed(3)} Å, less than minimum ${minDist.toFixed(3)} Å.`,
        suggestions: [
            `Reduce min_distance for ${pair} to ${(distance * 0.9).toFixed(3)} Å`,
            "Increase volume_factor to allow more space",
            "Try a different space group with lower symmetry"
        ]
    }),
    MODEL_NOT_AVAILABLE: (model) => ({
        message: `MLFF model '${model}' is not available or failed to load.`,
        suggestions: [
            `Install the required package: pip install ${model}`,
            "Check that the model is properly installed",
            "Try a different MLFF model (chgnet, m3gnet, mace)",
            "Verify Python environment has required dependencies"
        ]
    })
};
//# sourceMappingURL=errors.js.map