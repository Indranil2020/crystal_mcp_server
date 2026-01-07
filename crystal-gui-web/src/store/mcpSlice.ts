/**
 * MCP Slice - Redux state for MCP connection and tools
 */

import { createSlice, type PayloadAction } from '@reduxjs/toolkit';
import type { McpTool } from '../types';

type ConnectionStatus = 'disconnected' | 'connecting' | 'connected' | 'error';

interface McpState {
    connectionStatus: ConnectionStatus;
    tools: McpTool[];
    error: string | null;
    bridgeUrl: string;
    lastPingTime: number | null;
}

const initialState: McpState = {
    connectionStatus: 'disconnected',
    tools: [],
    error: null,
    bridgeUrl: 'http://localhost:8080',
    lastPingTime: null,
};

export const mcpSlice = createSlice({
    name: 'mcp',
    initialState,
    reducers: {
        setConnectionStatus: (state, action: PayloadAction<ConnectionStatus>) => {
            state.connectionStatus = action.payload;
        },
        setTools: (state, action: PayloadAction<McpTool[]>) => {
            state.tools = action.payload;
        },
        setError: (state, action: PayloadAction<string | null>) => {
            state.error = action.payload;
            if (action.payload) {
                state.connectionStatus = 'error';
            }
        },
        setBridgeUrl: (state, action: PayloadAction<string>) => {
            state.bridgeUrl = action.payload;
        },
        updatePingTime: (state) => {
            state.lastPingTime = Date.now();
        },
        disconnect: (state) => {
            state.connectionStatus = 'disconnected';
            state.tools = [];
            state.error = null;
        },
    },
});

export const {
    setConnectionStatus,
    setTools,
    setError,
    setBridgeUrl,
    updatePingTime,
    disconnect,
} = mcpSlice.actions;

export default mcpSlice.reducer;
