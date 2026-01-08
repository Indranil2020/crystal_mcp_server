/**
 * Tool Orchestrator - Coordinates LLM tool calls with MCP execution
 *
 * Ported from: crystal-gui/src/app.rs (send_chat_message, parse_tool_call, call_tool)
 *
 * DEBUG: Comprehensive logging at every step of the orchestration pipeline.
 */

import { v4 as uuidv4 } from 'uuid';
import { mcpClient } from './mcpClient';
import { llmClient } from './llmClient';
import { debug, debugError } from '../debug';
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
 *
 * DEBUG: Full logging enabled for pipeline tracing.
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

        debug('TOOL_ORCHESTRATOR', '═'.repeat(50));
        debug('TOOL_ORCHESTRATOR', '🔧 EXECUTING TOOL');
        debug('TOOL_ORCHESTRATOR', `  Tool: ${name}`);
        debug('TOOL_ORCHESTRATOR', `  Arguments: ${JSON.stringify(args)}`);

        // Call MCP tool
        debug('TOOL_ORCHESTRATOR', '  Calling MCP bridge...');
        const result = await mcpClient.callTool({
            name,
            arguments: args as Record<string, unknown>,
        });

        debug('TOOL_ORCHESTRATOR', `  MCP result received: isError=${result.isError}, content items=${result.content.length}`);

        // Parse structure from result
        debug('TOOL_ORCHESTRATOR', '  Parsing structure from result...');
        const structureResponse = mcpClient.parseStructureFromResult(result);

        debug('TOOL_ORCHESTRATOR', `  Structure parsing: success=${structureResponse.success}`);
        if (structureResponse.error) {
            debug('TOOL_ORCHESTRATOR', `  Parse error: ${structureResponse.error}`);
        }

        if (structureResponse.success && structureResponse.structure) {
            debug('TOOL_ORCHESTRATOR', '  ✓ Structure data extracted');
            debug('TOOL_ORCHESTRATOR', `    Atoms: ${structureResponse.structure.atoms?.length || 0}`);
            debug('TOOL_ORCHESTRATOR', `    Formula: ${structureResponse.structure.metadata?.formula || 'unknown'}`);

            // Create Structure object for state
            const structure = this.createStructureFromData(
                structureResponse.structure,
                name,
                structureResponse.source
            );

            debug('TOOL_ORCHESTRATOR', '  ✅ STRUCTURE CREATED FOR REDUX');
            debug('TOOL_ORCHESTRATOR', `    ID: ${structure.id}`);
            debug('TOOL_ORCHESTRATOR', `    Name: ${structure.name}`);
            debug('TOOL_ORCHESTRATOR', `    Atoms: ${structure.data.atoms.length}`);
            debug('TOOL_ORCHESTRATOR', `    Lattice: a=${structure.data.lattice.a.toFixed(2)}`);
            debug('TOOL_ORCHESTRATOR', '═'.repeat(50));

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

        debug('TOOL_ORCHESTRATOR', '  No structure in result, returning text response');
        debug('TOOL_ORCHESTRATOR', `  Text length: ${textContent.length}`);
        debug('TOOL_ORCHESTRATOR', '═'.repeat(50));

        return {
            success: true,
            toolName: name,
            arguments: args as Record<string, unknown>,
            textResponse: textContent,
        };
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
        const structureId = uuidv4();

        debug('TOOL_ORCHESTRATOR', `  Creating Structure object:`);
        debug('TOOL_ORCHESTRATOR', `    ID: ${structureId}`);
        debug('TOOL_ORCHESTRATOR', `    Formula: ${formula}`);
        debug('TOOL_ORCHESTRATOR', `    Tool: ${toolName}`);

        // Validate critical data
        if (!data.lattice) {
            debugError('TOOL_ORCHESTRATOR', 'Missing lattice data!', 'createStructureFromData');
        }
        if (!data.atoms || data.atoms.length === 0) {
            debugError('TOOL_ORCHESTRATOR', 'Missing or empty atoms array!', 'createStructureFromData');
        }

        // Log first atom for verification
        if (data.atoms && data.atoms.length > 0) {
            const first = data.atoms[0];
            debug('TOOL_ORCHESTRATOR', `    First atom: ${first.element} coords=${JSON.stringify(first.coords)} cartesian=${JSON.stringify(first.cartesian)}`);
        }

        return {
            id: structureId,
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
