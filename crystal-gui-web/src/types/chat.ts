/**
 * Chat Types - TypeScript interfaces for LLM chat communication
 */

/** Chat message role */
export type MessageRole = 'user' | 'assistant' | 'system' | 'tool';

/** Tool call from LLM response */
export interface ToolCall {
    function: {
        name: string;
        arguments: Record<string, unknown>;
    };
}

/** Chat message */
export interface ChatMessage {
    id: string;
    role: MessageRole;
    content: string;
    timestamp: number;
    toolCalls?: ToolCall[];
    toolResult?: {
        toolName: string;
        success: boolean;
        structureId?: string;
    };
}

/** Ollama chat request format */
export interface OllamaChatRequest {
    model: string;
    messages: Array<{
        role: string;
        content: string;
        tool_calls?: ToolCall[];
    }>;
    stream: boolean;
    tools?: OllamaToolDef[];
    format?: string;
}

/** Ollama tool definition */
export interface OllamaToolDef {
    type: 'function';
    function: {
        name: string;
        description: string;
        parameters: object;
    };
}

/** Ollama chat response */
export interface OllamaChatResponse {
    model: string;
    message: {
        role: string;
        content: string;
        tool_calls?: ToolCall[];
    };
    done: boolean;
}

/** Chat result - either tool calls or text */
export type ChatResult =
    | { type: 'tool_calls'; toolCalls: ToolCall[] }
    | { type: 'text'; content: string };

/** Tool execution state */
export type ToolExecutionStatus = 'pending' | 'running' | 'success' | 'error';

/** Tool execution tracker */
export interface ToolExecution {
    id: string;
    toolName: string;
    arguments: Record<string, unknown>;
    status: ToolExecutionStatus;
    startTime: number;
    endTime?: number;
    result?: string;
    error?: string;
    structureId?: string;
}
