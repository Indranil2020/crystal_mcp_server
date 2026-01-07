/**
 * ChatPanel - LLM Chat Interface with Tool Visualization
 * 
 * Features:
 * - Message history display
 * - Tool call visualization (loading, success, error states)
 * - Natural language input
 * - Auto-scroll with manual override
 */

import { useState, useRef, useEffect } from 'react';
import { v4 as uuidv4 } from 'uuid';
import { useAppDispatch, useAppSelector } from '../../store/hooks';
import {
    addMessage,
    setProcessing,
    startToolExecution,
    updateToolExecution,
} from '../../store/chatSlice';
import { addStructure } from '../../store/structureSlice';
import { toolOrchestrator } from '../../services';
import type { ChatMessage } from '../../types';
import ChatMessageItem from './ChatMessageItem';
import ToolCallCard from './ToolCallCard';

export default function ChatPanel() {
    const dispatch = useAppDispatch();
    const { messages, toolExecutions, isProcessing } = useAppSelector(state => state.chat);
    const { tools, connectionStatus } = useAppSelector(state => state.mcp);

    const [input, setInput] = useState('');
    const messagesEndRef = useRef<HTMLDivElement>(null);
    const inputRef = useRef<HTMLTextAreaElement>(null);

    // Auto-scroll to bottom on new messages
    useEffect(() => {
        messagesEndRef.current?.scrollIntoView({ behavior: 'smooth' });
    }, [messages, toolExecutions]);

    // Handle message submission
    async function handleSubmit(e: React.FormEvent) {
        e.preventDefault();
        if (!input.trim() || isProcessing) return;

        const userMessage: ChatMessage = {
            id: uuidv4(),
            role: 'user',
            content: input.trim(),
            timestamp: Date.now(),
        };

        dispatch(addMessage(userMessage));
        setInput('');
        dispatch(setProcessing(true));

        try {
            // Initialize tools if needed
            if (tools.length === 0) {
                await toolOrchestrator.initialize();
            }

            // Build chat history for LLM
            const history = messages.map(m => ({
                role: m.role,
                content: m.content,
            }));

            // Process through LLM → MCP pipeline
            const result = await toolOrchestrator.processMessage(userMessage.content, history);

            if (result.toolName) {
                // Tool was called - track execution
                const executionId = uuidv4();
                dispatch(startToolExecution({
                    id: executionId,
                    toolName: result.toolName,
                    arguments: result.arguments || {},
                }));

                if (result.success && result.structure) {
                    // Add structure to state
                    dispatch(addStructure(result.structure));
                    dispatch(updateToolExecution({
                        id: executionId,
                        status: 'success',
                        structureId: result.structure.id,
                        result: `Generated ${result.structure.name}`,
                    }));

                    // Add assistant message
                    dispatch(addMessage({
                        id: uuidv4(),
                        role: 'assistant',
                        content: `Successfully generated structure: ${result.structure.name}`,
                        timestamp: Date.now(),
                        toolResult: {
                            toolName: result.toolName,
                            success: true,
                            structureId: result.structure.id,
                        },
                    }));
                } else if (result.success && result.textResponse) {
                    dispatch(updateToolExecution({
                        id: executionId,
                        status: 'success',
                        result: result.textResponse.slice(0, 100),
                    }));

                    dispatch(addMessage({
                        id: uuidv4(),
                        role: 'assistant',
                        content: result.textResponse,
                        timestamp: Date.now(),
                    }));
                } else {
                    dispatch(updateToolExecution({
                        id: executionId,
                        status: 'error',
                        error: result.error || 'Unknown error',
                    }));
                }
            } else if (result.textResponse) {
                // Plain text response (no tool call)
                dispatch(addMessage({
                    id: uuidv4(),
                    role: 'assistant',
                    content: result.textResponse,
                    timestamp: Date.now(),
                }));
            }
        } catch (error) {
            console.error('[ChatPanel] Error:', error);
            dispatch(addMessage({
                id: uuidv4(),
                role: 'assistant',
                content: `Error: ${error instanceof Error ? error.message : 'Unknown error'}`,
                timestamp: Date.now(),
            }));
        } finally {
            dispatch(setProcessing(false));
        }
    }

    // Handle Enter key (Shift+Enter for newline)
    function handleKeyDown(e: React.KeyboardEvent) {
        if (e.key === 'Enter' && !e.shiftKey) {
            e.preventDefault();
            handleSubmit(e);
        }
    }

    const isConnected = connectionStatus === 'connected';

    return (
        <div className="flex flex-col h-full bg-slate-800">
            {/* Header */}
            <div className="px-4 py-3 border-b border-slate-700 flex items-center justify-between">
                <h2 className="text-sm font-semibold text-slate-200">Crystal Assistant</h2>
                <div className={`w-2 h-2 rounded-full ${isConnected ? 'bg-green-500' : 'bg-red-500'}`} />
            </div>

            {/* Messages */}
            <div className="flex-1 overflow-y-auto p-4 space-y-4">
                {messages.length === 0 ? (
                    <div className="text-center text-slate-500 mt-8">
                        <p className="text-lg mb-2">👋 Welcome to Crystal GUI</p>
                        <p className="text-sm">Ask me to generate molecules or crystal structures.</p>
                        <p className="text-sm text-slate-600 mt-4">Examples:</p>
                        <p className="text-sm text-blue-400">"Generate benzene"</p>
                        <p className="text-sm text-blue-400">"Create silicon with spacegroup 227"</p>
                        <p className="text-sm text-blue-400">"Build caffeine-aspirin dimer"</p>
                    </div>
                ) : (
                    messages.map(message => (
                        <ChatMessageItem key={message.id} message={message} />
                    ))
                )}

                {/* Active Tool Executions */}
                {toolExecutions
                    .filter(e => e.status === 'running')
                    .map(execution => (
                        <ToolCallCard key={execution.id} execution={execution} />
                    ))}

                {/* Processing indicator */}
                {isProcessing && (
                    <div className="flex items-center gap-2 text-slate-400">
                        <div className="w-2 h-2 bg-blue-500 rounded-full animate-pulse" />
                        <span className="text-sm">Thinking...</span>
                    </div>
                )}

                <div ref={messagesEndRef} />
            </div>

            {/* Input */}
            <form onSubmit={handleSubmit} className="p-4 border-t border-slate-700">
                <div className="flex gap-2">
                    <textarea
                        ref={inputRef}
                        value={input}
                        onChange={e => setInput(e.target.value)}
                        onKeyDown={handleKeyDown}
                        placeholder={isConnected ? "Ask me to generate a structure..." : "Connecting..."}
                        disabled={!isConnected || isProcessing}
                        rows={1}
                        className="flex-1 bg-slate-700 text-slate-100 rounded-lg px-4 py-2 
                       placeholder-slate-500 resize-none focus:outline-none 
                       focus:ring-2 focus:ring-blue-500 disabled:opacity-50"
                    />
                    <button
                        type="submit"
                        disabled={!isConnected || isProcessing || !input.trim()}
                        className="px-4 py-2 bg-blue-600 text-white rounded-lg font-medium
                       hover:bg-blue-500 disabled:opacity-50 disabled:cursor-not-allowed
                       transition-colors"
                    >
                        Send
                    </button>
                </div>
            </form>
        </div>
    );
}
