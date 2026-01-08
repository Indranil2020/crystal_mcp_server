/**
 * Debug Configuration
 * 
 * Enable debug logging throughout the application to trace request flow
 */

// Set to true to enable verbose debug logging
export const DEBUG_ENABLED = true;

// Debug log levels - ALL enabled for thorough debugging
export const DEBUG_LEVELS = {
    MCP_CLIENT: true,         // MCP bridge communication
    LLM_CLIENT: true,         // Ollama API calls
    TOOL_ORCHESTRATOR: true,  // Tool call orchestration (renamed from TOOL_ORCHESTRATOR)
    REDUX_STORE: true,        // State changes
    CONVERTERS: true,         // Structure conversion
    VIEWERS: true,            // MolStar/Kekule loading
    CHAT_PANEL: true,         // Chat UI interactions
};

// Timestamp helper
function timestamp(): string {
    return new Date().toISOString().split('T')[1].slice(0, 12);
}

// Debug logger with module tag
export function debug(module: keyof typeof DEBUG_LEVELS, ...args: unknown[]): void {
    if (!DEBUG_ENABLED) return;
    if (!DEBUG_LEVELS[module]) return;

    const prefix = `[${timestamp()}][${module}]`;
    console.log(prefix, ...args);
}

// Error logger (always enabled)
export function debugError(module: keyof typeof DEBUG_LEVELS, error: unknown, context?: string): void {
    const prefix = `[${timestamp()}][${module}][ERROR]`;
    console.error(prefix, context || '', error);
}

// Performance timer
export function debugTimer(module: keyof typeof DEBUG_LEVELS, label: string): () => void {
    if (!DEBUG_ENABLED || !DEBUG_LEVELS[module]) return () => { };

    const start = performance.now();
    debug(module, `⏱️ START: ${label}`);

    return () => {
        const duration = (performance.now() - start).toFixed(2);
        debug(module, `⏱️ END: ${label} (${duration}ms)`);
    };
}

// Object inspector (truncated)
export function debugInspect(obj: unknown, maxLength = 500): string {
    try {
        const str = JSON.stringify(obj, null, 2);
        if (str.length > maxLength) {
            return str.slice(0, maxLength) + '... [TRUNCATED]';
        }
        return str;
    } catch {
        return String(obj);
    }
}
