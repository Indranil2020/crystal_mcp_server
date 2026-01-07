/**
 * MCP Types - TypeScript interfaces for MCP communication
 */

/** MCP Tool definition */
export interface McpTool {
    name: string;
    description?: string;
    inputSchema?: object;
}

/** Content item in MCP response */
export interface McpContentItem {
    type: 'text' | 'image' | 'resource';
    text?: string;
    data?: string;
    mimeType?: string;
}

/** MCP Tool call result */
export interface McpToolResult {
    content: McpContentItem[];
    isError: boolean;
}

/** JSON-RPC request format */
export interface JsonRpcRequest {
    jsonrpc: '2.0';
    id: number | string;
    method: string;
    params?: object;
}

/** JSON-RPC response format */
export interface JsonRpcResponse {
    jsonrpc: '2.0';
    id: number | string;
    result?: unknown;
    error?: {
        code: number;
        message: string;
        data?: unknown;
    };
}

/** MCP Server initialization params */
export interface McpInitializeParams {
    protocolVersion: string;
    capabilities: object;
    clientInfo: {
        name: string;
        version: string;
    };
}

/** Parsed structure from MCP response */
export interface McpStructureResponse {
    success: boolean;
    structure?: StructureData;
    source?: string;
    error?: string;
}

/** Structure data from Python backend */
export interface StructureData {
    lattice: LatticeData;
    atoms: AtomData[];
    space_group?: SpaceGroupData;
    metadata?: StructureMetadata;
}

export interface LatticeData {
    a: number;
    b: number;
    c: number;
    alpha: number;
    beta: number;
    gamma: number;
    matrix: number[][];
    volume?: number;
}

export interface AtomData {
    element: string;
    coords: [number, number, number];
    cartesian?: [number, number, number];
    species?: { element: string; occupation: number }[];
    label?: string;
}

export interface SpaceGroupData {
    number: number;
    symbol: string;
    international?: string;
}

export interface StructureMetadata {
    formula: string;
    natoms: number;
    n_molecules?: number;
    stacking_type?: string;
    pbc?: [boolean, boolean, boolean];
}
