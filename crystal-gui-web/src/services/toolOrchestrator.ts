/**
 * Tool Orchestrator - Coordinates LLM tool calls with MCP execution
 * 
 * Ported from: crystal-gui/src/app.rs (send_chat_message, parse_tool_call, call_tool)
 */

import { v4 as uuidv4 } from 'uuid';
import { mcpClient } from './mcpClient';
import { llmClient } from './llmClient';
import type { McpTool, ChatResult, ToolCall, Structure, StructureData } from '../types';

export interface ToolOrchestrationResult {
    success: boolean;
    toolName?: string;
    arguments?: Record<string, unknown>;
    structure?: Structure;
    textResponse?: string;
    error?: string;
}

/**
 * Tool Orchestrator Service
 * 
 * Handles the complete flow:
 * 1. User message → LLM with tool schemas
 * 2. Parse tool_calls from LLM response
 * 3. Execute via MCP client
 * 4. Parse structure from result
 */
export class ToolOrchestrator {
    private tools: McpTool[] = [];

    /**
     * Initialize by fetching available tools
     */
    async initialize(): Promise<McpTool[]> {
        this.tools = await mcpClient.listTools();
        console.log('[ToolOrchestrator] Initialized with', this.tools.length, 'tools');
        return this.tools;
    }

    /**
     * Set tools directly (for external injection from Redux state)
     */
    setTools(tools: McpTool[]): void {
        this.tools = tools;
        console.log('[ToolOrchestrator] Tools set externally:', this.tools.length, 'tools');
    }

    /**
     * Get cached tools
     */
    getTools(): McpTool[] {
        return this.tools;
    }

    /**
     * Process a user message through the LLM → MCP pipeline
     */
    async processMessage(
        userMessage: string,
        chatHistory: Array<{ role: string; content: string }>
    ): Promise<ToolOrchestrationResult> {
        // Build messages array
        const messages = [
            ...chatHistory,
            { role: 'user', content: userMessage },
        ];

        // Send to LLM with tool schemas
        const llmResult: ChatResult = await llmClient.chatWithTools(messages, this.tools);

        if (llmResult.type === 'text') {
            // LLM responded with text, no tool call
            return {
                success: true,
                textResponse: llmResult.content,
            };
        }

        // Process tool calls
        const toolCall = llmResult.toolCalls[0]; // Take first tool call
        if (!toolCall) {
            return {
                success: false,
                error: 'LLM returned empty tool calls array',
            };
        }

        return await this.executeTool(toolCall);
    }

    /**
     * Execute a single tool call via MCP
     */
    async executeTool(toolCall: ToolCall): Promise<ToolOrchestrationResult> {
        const { name, arguments: args } = toolCall.function;

        console.debug('[ToolOrchestrator] Executing tool:', name, args);

        try {
            // Call MCP tool
            const result = await mcpClient.callTool({
                name,
                arguments: args as Record<string, unknown>,
            });

            // Parse structure from result
            const structureResponse = mcpClient.parseStructureFromResult(result);

            if (structureResponse.success && structureResponse.structure) {
                // Create Structure object for state
                const structure = this.createStructureFromData(
                    structureResponse.structure,
                    name,
                    structureResponse.source
                );

                return {
                    success: true,
                    toolName: name,
                    arguments: args as Record<string, unknown>,
                    structure,
                };
            }

            // Tool succeeded but no structure (e.g., analysis tool)
            const textContent = result.content
                .filter(c => c.type === 'text')
                .map(c => c.text)
                .join('\n');

            return {
                success: true,
                toolName: name,
                arguments: args as Record<string, unknown>,
                textResponse: textContent,
            };
        } catch (error) {
            return {
                success: false,
                toolName: name,
                arguments: args as Record<string, unknown>,
                error: error instanceof Error ? error.message : 'Unknown error',
            };
        }
    }

    /**
     * Create a Structure object from raw structure data
     */
    private createStructureFromData(
        data: StructureData,
        toolName: string,
        _source?: string
    ): Structure {
        const formula = data.metadata?.formula || 'Unknown';

        return {
            id: uuidv4(),
            name: `${formula} (${toolName})`,
            data,
            source: 'mcp',
            createdAt: Date.now(),
            modifiedAt: Date.now(),
            visible: true,
        };
    }
}

// Singleton instance
export const toolOrchestrator = new ToolOrchestrator();
