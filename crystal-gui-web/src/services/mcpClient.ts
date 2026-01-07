/**
 * MCP Client Service - HTTP bridge communication
 * 
 * Ported from: crystal-gui/src/mcp_client.rs
 * Communicates with the FastAPI bridge server
 * 
 * DEBUG: Comprehensive logging enabled
 */

import type { McpTool, McpToolResult, McpStructureResponse, StructureData } from '../types';
import { debug, debugError, debugTimer, debugInspect } from '../debug';

const DEFAULT_BRIDGE_URL = 'http://localhost:8080';

interface McpCallParams {
    name: string;
    arguments: Record<string, unknown>;
}

/**
 * MCP Client for communicating with the bridge server
 */
export class McpClient {
    private bridgeUrl: string;
    private initialized: boolean = false;

    constructor(bridgeUrl: string = DEFAULT_BRIDGE_URL) {
        this.bridgeUrl = bridgeUrl;
        debug('MCP_CLIENT', `Created MCP client with bridge URL: ${bridgeUrl}`);
    }

    /**
     * Check connection to the bridge server
     */
    async checkConnection(): Promise<boolean> {
        const endTimer = debugTimer('MCP_CLIENT', 'checkConnection');
        debug('MCP_CLIENT', `Checking connection to: ${this.bridgeUrl}/health`);

        try {
            const response = await fetch(`${this.bridgeUrl}/health`, {
                method: 'GET',
                headers: { 'Content-Type': 'application/json' },
            });

            const ok = response.ok;
            debug('MCP_CLIENT', `Health check result: ${ok ? 'OK' : 'FAILED'} (status: ${response.status})`);

            if (ok) {
                const data = await response.json();
                debug('MCP_CLIENT', `Health response:`, debugInspect(data));
            }

            endTimer();
            return ok;
        } catch (error) {
            debugError('MCP_CLIENT', error, 'Health check failed');
            endTimer();
            return false;
        }
    }

    /**
     * Initialize the MCP connection through the bridge
     */
    async initialize(): Promise<void> {
        const endTimer = debugTimer('MCP_CLIENT', 'initialize');
        debug('MCP_CLIENT', 'Initializing MCP connection...');

        try {
            const body = {
                protocolVersion: '2024-11-05',
                capabilities: {},
                clientInfo: {
                    name: 'crystal-gui-web',
                    version: '1.0.0',
                },
            };
            debug('MCP_CLIENT', `Sending initialize request:`, debugInspect(body));

            const response = await fetch(`${this.bridgeUrl}/mcp/initialize`, {
                method: 'POST',
                headers: { 'Content-Type': 'application/json' },
                body: JSON.stringify(body),
            });

            debug('MCP_CLIENT', `Initialize response status: ${response.status} ${response.statusText}`);

            if (!response.ok) {
                const errorText = await response.text();
                debugError('MCP_CLIENT', errorText, 'Initialize failed');
                throw new Error(`Failed to initialize MCP: ${response.statusText} - ${errorText}`);
            }

            const result = await response.json();
            debug('MCP_CLIENT', `Initialize result:`, debugInspect(result));

            this.initialized = true;
            debug('MCP_CLIENT', '✅ MCP initialized successfully');
            endTimer();
        } catch (error) {
            debugError('MCP_CLIENT', error, 'Initialize exception');
            endTimer();
            throw error;
        }
    }

    /**
     * List available tools from the MCP server
     */
    async listTools(): Promise<McpTool[]> {
        const endTimer = debugTimer('MCP_CLIENT', 'listTools');

        if (!this.initialized) {
            debug('MCP_CLIENT', 'Not initialized, calling initialize first...');
            await this.initialize();
        }

        debug('MCP_CLIENT', 'Fetching tool list...');

        try {
            const response = await fetch(`${this.bridgeUrl}/mcp/tools`, {
                method: 'GET',
                headers: { 'Content-Type': 'application/json' },
            });

            debug('MCP_CLIENT', `Tools response status: ${response.status}`);

            if (!response.ok) {
                const errorText = await response.text();
                debugError('MCP_CLIENT', errorText, 'List tools failed');
                throw new Error(`Failed to list tools: ${response.statusText}`);
            }

            const data = await response.json();
            const tools = data.tools || [];
            debug('MCP_CLIENT', `✅ Received ${tools.length} tools:`, tools.map((t: McpTool) => t.name));

            endTimer();
            return tools;
        } catch (error) {
            debugError('MCP_CLIENT', error, 'List tools exception');
            endTimer();
            throw error;
        }
    }

    /**
     * Call a tool through the MCP bridge
     */
    async callTool(params: McpCallParams): Promise<McpToolResult> {
        const endTimer = debugTimer('MCP_CLIENT', `callTool(${params.name})`);

        if (!this.initialized) {
            debug('MCP_CLIENT', 'Not initialized, calling initialize first...');
            await this.initialize();
        }

        debug('MCP_CLIENT', `📞 Calling tool: ${params.name}`);
        debug('MCP_CLIENT', `Arguments:`, debugInspect(params.arguments));

        try {
            const response = await fetch(`${this.bridgeUrl}/mcp/call`, {
                method: 'POST',
                headers: { 'Content-Type': 'application/json' },
                body: JSON.stringify({
                    name: params.name,
                    arguments: params.arguments,
                }),
            });

            debug('MCP_CLIENT', `Tool call response status: ${response.status}`);

            if (!response.ok) {
                const errorText = await response.text();
                debugError('MCP_CLIENT', errorText, `Tool ${params.name} call failed`);
                throw new Error(`Tool call failed: ${response.statusText} - ${errorText}`);
            }

            const result: McpToolResult = await response.json();
            debug('MCP_CLIENT', `✅ Tool ${params.name} result:`, debugInspect(result, 1000));

            endTimer();
            return result;
        } catch (error) {
            debugError('MCP_CLIENT', error, `Tool ${params.name} exception`);
            endTimer();
            throw error;
        }
    }

    /**
     * Parse structure from MCP tool result
     * Extracts <json-data> tag and parses structure
     */
    parseStructureFromResult(result: McpToolResult): McpStructureResponse {
        debug('MCP_CLIENT', 'Parsing structure from result...');
        debug('MCP_CLIENT', `Result isError: ${result.isError}`);
        debug('MCP_CLIENT', `Result content items: ${result.content.length}`);

        if (result.isError) {
            debug('MCP_CLIENT', '❌ Tool result has error flag');
            return { success: false, error: 'Tool execution failed' };
        }

        for (let i = 0; i < result.content.length; i++) {
            const item = result.content[i];
            debug('MCP_CLIENT', `Content item ${i}: type=${item.type}`);

            if (item.type === 'text' && item.text) {
                debug('MCP_CLIENT', `Text content length: ${item.text.length}`);

                // Look for <json-data> tag
                const jsonDataMatch = item.text.match(/<json-data>([\s\S]*?)<\/json-data>/);
                if (jsonDataMatch) {
                    debug('MCP_CLIENT', `Found <json-data> tag, content length: ${jsonDataMatch[1].length}`);
                    try {
                        const parsed = JSON.parse(jsonDataMatch[1]);
                        const structure: StructureData = parsed.structure || parsed;
                        debug('MCP_CLIENT', `✅ Parsed structure: ${structure.atoms?.length} atoms`);
                        return {
                            success: true,
                            structure,
                            source: parsed.source || 'mcp',
                        };
                    } catch (e) {
                        debugError('MCP_CLIENT', e, 'Failed to parse structure JSON');
                    }
                } else {
                    debug('MCP_CLIENT', 'No <json-data> tag found in this content item');
                    debug('MCP_CLIENT', `Text preview: ${item.text.slice(0, 200)}...`);
                }
            }
        }

        debug('MCP_CLIENT', '❌ No structure data found in response');
        return { success: false, error: 'No structure data found in response' };
    }

    /**
     * Get bridge URL
     */
    getBridgeUrl(): string {
        return this.bridgeUrl;
    }

    /**
     * Set bridge URL
     */
    setBridgeUrl(url: string): void {
        debug('MCP_CLIENT', `Changing bridge URL from ${this.bridgeUrl} to ${url}`);
        this.bridgeUrl = url;
        this.initialized = false;
    }

    /**
     * Check if initialized
     */
    isInitialized(): boolean {
        return this.initialized;
    }
}

// Singleton instance
export const mcpClient = new McpClient();
