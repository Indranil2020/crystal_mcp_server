/**
 * Chat Slice - Redux state for chat messages and tool executions
 */

import { createSlice, type PayloadAction } from '@reduxjs/toolkit';
import type { ChatMessage, ToolExecution, ToolExecutionStatus } from '../types';

interface ChatState {
    messages: ChatMessage[];
    toolExecutions: ToolExecution[];
    isProcessing: boolean;
    currentInput: string;
    error: string | null;
}

const initialState: ChatState = {
    messages: [],
    toolExecutions: [],
    isProcessing: false,
    currentInput: '',
    error: null,
};

export const chatSlice = createSlice({
    name: 'chat',
    initialState,
    reducers: {
        addMessage: (state, action: PayloadAction<ChatMessage>) => {
            state.messages.push(action.payload);
        },
        updateMessage: (state, action: PayloadAction<{ id: string; updates: Partial<ChatMessage> }>) => {
            const index = state.messages.findIndex(m => m.id === action.payload.id);
            if (index !== -1) {
                state.messages[index] = { ...state.messages[index], ...action.payload.updates };
            }
        },
        clearMessages: (state) => {
            state.messages = [];
            state.toolExecutions = [];
        },
        setProcessing: (state, action: PayloadAction<boolean>) => {
            state.isProcessing = action.payload;
        },
        setCurrentInput: (state, action: PayloadAction<string>) => {
            state.currentInput = action.payload;
        },
        setError: (state, action: PayloadAction<string | null>) => {
            state.error = action.payload;
        },
        // Tool execution tracking
        startToolExecution: (state, action: PayloadAction<{ id: string; toolName: string; arguments: Record<string, unknown> }>) => {
            state.toolExecutions.push({
                id: action.payload.id,
                toolName: action.payload.toolName,
                arguments: action.payload.arguments,
                status: 'running',
                startTime: Date.now(),
            });
        },
        updateToolExecution: (state, action: PayloadAction<{ id: string; status: ToolExecutionStatus; result?: string; error?: string; structureId?: string }>) => {
            const execution = state.toolExecutions.find(e => e.id === action.payload.id);
            if (execution) {
                execution.status = action.payload.status;
                execution.endTime = Date.now();
                if (action.payload.result) execution.result = action.payload.result;
                if (action.payload.error) execution.error = action.payload.error;
                if (action.payload.structureId) execution.structureId = action.payload.structureId;
            }
        },
    },
});

export const {
    addMessage,
    updateMessage,
    clearMessages,
    setProcessing,
    setCurrentInput,
    setError,
    startToolExecution,
    updateToolExecution,
} = chatSlice.actions;

export default chatSlice.reducer;
