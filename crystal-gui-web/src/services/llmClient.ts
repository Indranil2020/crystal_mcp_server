/**
 * LLM Client Service - Ollama API communication
 * 
 * Ported from: crystal-gui/src/llm_client.rs
 * Communicates with local Ollama server for LLM inference
 * 
 * DEBUG: Comprehensive logging enabled
 */

import type { ChatResult, McpTool, ToolCall } from '../types';
import { debug, debugError, debugTimer, debugInspect } from '../debug';

const DEFAULT_OLLAMA_URL = 'http://localhost:11434';
const DEFAULT_MODEL = 'qwen2.5:7b';

interface OllamaToolDef {
    type: 'function';
    function: {
        name: string;
        description: string;
        parameters: object;
    };
}

interface OllamaChatRequest {
    model: string;
    messages: Array<{ role: string; content: string }>;
    stream: boolean;
    tools?: OllamaToolDef[];
    format?: string;
    options?: {
        temperature?: number;
        top_p?: number;
        seed?: number;
    };
}

interface OllamaChatResponse {
    model: string;
    message: {
        role: string;
        content: string;
        tool_calls?: ToolCall[];
    };
    done: boolean;
}

/**
 * Convert MCP tools to Ollama format
 */
function mcpToolsToOllamaFormat(tools: McpTool[]): OllamaToolDef[] {
    debug('LLM_CLIENT', `Converting ${tools.length} MCP tools to Ollama format`);
    return tools.map(tool => ({
        type: 'function',
        function: {
            name: tool.name,
            description: tool.description || '',
            parameters: tool.inputSchema || { type: 'object', properties: {} },
        },
    }));
}

/**
 * LLM Client for communicating with Ollama
 */
export class LlmClient {
    private baseUrl: string;
    private model: string;
    private availableModels: string[] = [];
    private _temperature: number = 0; // Default 0 for deterministic tool calling

    constructor(baseUrl: string = DEFAULT_OLLAMA_URL, model: string = DEFAULT_MODEL) {
        this.baseUrl = baseUrl;
        this.model = model;
        debug('LLM_CLIENT', `Created LLM client: url=${baseUrl}, model=${model}, temperature=${this._temperature}`);
    }

    /**
     * Check connection to Ollama and fetch available models
     */
    async checkConnection(): Promise<boolean> {
        const endTimer = debugTimer('LLM_CLIENT', 'checkConnection');
        debug('LLM_CLIENT', `Checking connection to: ${this.baseUrl}/api/tags`);

        try {
            const response = await fetch(`${this.baseUrl}/api/tags`);
            debug('LLM_CLIENT', `Response status: ${response.status}`);

            if (response.ok) {
                const data = await response.json();
                this.availableModels = (data.models || []).map((m: { name: string }) => m.name);
                debug('LLM_CLIENT', `✅ Connected. Available models: ${this.availableModels.join(', ')}`);
                endTimer();
                return true;
            }
            debug('LLM_CLIENT', '❌ Connection failed: not OK response');
            endTimer();
            return false;
        } catch (error) {
            debugError('LLM_CLIENT', error, 'Connection check failed');
            endTimer();
            return false;
        }
    }

    /**
     * Get available models
     */
    getAvailableModels(): string[] {
        return this.availableModels;
    }

    /**
     * Set model to use
     */
    setModel(model: string): void {
        debug('LLM_CLIENT', `Changing model from ${this.model} to ${model}`);
        this.model = model;
    }

    /**
     * Get current model
     */
    getModel(): string {
        return this.model;
    }

    /**
     * Set temperature (0-1). Lower = more deterministic, higher = more creative.
     * Default 0 for precise tool calling.
     */
    setTemperature(temp: number): void {
        this._temperature = Math.max(0, Math.min(1, temp));
        debug('LLM_CLIENT', `Temperature set to ${this._temperature}`);
    }

    /**
     * Get current temperature
     */
    getTemperature(): number {
        return this._temperature;
    }

    /**
     * Chat with native Ollama tool calling
     * 
     * @param messages - Chat message history
     * @param tools - MCP tools to make available
     * @returns ChatResult - either tool calls or text response
     */
    async chatWithTools(
        messages: Array<{ role: string; content: string }>,
        tools: McpTool[]
    ): Promise<ChatResult> {
        const endTimer = debugTimer('LLM_CLIENT', 'chatWithTools');

        debug('LLM_CLIENT', '='.repeat(60));
        debug('LLM_CLIENT', '📤 STARTING LLM CHAT WITH TOOLS');
        debug('LLM_CLIENT', `Model: ${this.model}`);
        debug('LLM_CLIENT', `Message count: ${messages.length}`);
        debug('LLM_CLIENT', `Tool count: ${tools.length}`);
        debug('LLM_CLIENT', `Last message: "${messages[messages.length - 1]?.content?.slice(0, 100)}..."`);

        const ollamaTools = mcpToolsToOllamaFormat(tools);
        debug('LLM_CLIENT', `Converted tools: ${ollamaTools.map(t => t.function.name).join(', ')}`);

        const request: OllamaChatRequest = {
            model: this.model,
            messages: messages,
            stream: false,
            tools: ollamaTools.length > 0 ? ollamaTools : undefined,
            // DO NOT use format: 'json' - it DISABLES native tool calling!
            format: undefined,
            options: {
                temperature: this._temperature,
            },
        };

        debug('LLM_CLIENT', `Full request:`, debugInspect(request, 2000));

        try {
            debug('LLM_CLIENT', `⏳ Sending POST to ${this.baseUrl}/api/chat...`);
            const startTime = performance.now();

            const response = await fetch(`${this.baseUrl}/api/chat`, {
                method: 'POST',
                headers: { 'Content-Type': 'application/json' },
                body: JSON.stringify(request),
            });

            const elapsed = ((performance.now() - startTime) / 1000).toFixed(2);
            debug('LLM_CLIENT', `📥 Response received in ${elapsed}s, status: ${response.status}`);

            if (!response.ok) {
                const errorText = await response.text();
                debugError('LLM_CLIENT', errorText, 'LLM request failed');
                throw new Error(`LLM request failed: ${response.statusText} - ${errorText}`);
            }

            const data: OllamaChatResponse = await response.json();
            debug('LLM_CLIENT', `Response model: ${data.model}`);
            debug('LLM_CLIENT', `Response done: ${data.done}`);
            debug('LLM_CLIENT', `Message role: ${data.message.role}`);
            debug('LLM_CLIENT', `Message content length: ${data.message.content?.length || 0}`);
            debug('LLM_CLIENT', `Tool calls present: ${!!data.message.tool_calls}`);

            if (data.message.tool_calls) {
                debug('LLM_CLIENT', `Tool calls count: ${data.message.tool_calls.length}`);
            }

            // Check for tool calls first
            if (data.message.tool_calls && data.message.tool_calls.length > 0) {
                debug('LLM_CLIENT', '✅ LLM returned TOOL CALLS:');
                data.message.tool_calls.forEach((tc, i) => {
                    debug('LLM_CLIENT', `  [${i}] ${tc.function.name}(${debugInspect(tc.function.arguments)})`);
                });
                endTimer();
                return {
                    type: 'tool_calls',
                    toolCalls: data.message.tool_calls,
                };
            }

            // Fall back to text content
            debug('LLM_CLIENT', '📝 LLM returned TEXT response:');
            debug('LLM_CLIENT', `Content preview: ${data.message.content?.slice(0, 200)}...`);
            endTimer();
            return {
                type: 'text',
                content: data.message.content,
            };
        } catch (error) {
            debugError('LLM_CLIENT', error, 'chatWithTools exception');
            debug('LLM_CLIENT', '='.repeat(60));
            endTimer();
            throw error;
        }
    }

    /**
     * Plain chat without tools (legacy method)
     */
    async chat(
        messages: Array<{ role: string; content: string }>,
        format?: string
    ): Promise<string> {
        const endTimer = debugTimer('LLM_CLIENT', 'chat (plain)');
        debug('LLM_CLIENT', `Plain chat with ${messages.length} messages`);

        const request: OllamaChatRequest = {
            model: this.model,
            messages,
            stream: false,
            format,
            options: {
                temperature: this._temperature,
            },
        };

        try {
            const response = await fetch(`${this.baseUrl}/api/chat`, {
                method: 'POST',
                headers: { 'Content-Type': 'application/json' },
                body: JSON.stringify(request),
            });

            debug('LLM_CLIENT', `Response status: ${response.status}`);

            if (!response.ok) {
                const errorText = await response.text();
                debugError('LLM_CLIENT', errorText, 'Plain chat failed');
                throw new Error(`LLM request failed: ${response.statusText}`);
            }

            const data: OllamaChatResponse = await response.json();
            debug('LLM_CLIENT', `Response content length: ${data.message.content?.length}`);
            endTimer();
            return data.message.content;
        } catch (error) {
            debugError('LLM_CLIENT', error, 'Plain chat exception');
            endTimer();
            throw error;
        }
    }
}

// Singleton instance
export const llmClient = new LlmClient();
