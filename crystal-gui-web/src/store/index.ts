/**
 * Redux Store Configuration
 */

import { configureStore } from '@reduxjs/toolkit';
import chatReducer from './chatSlice';
import structureReducer from './structureSlice';
import mcpReducer from './mcpSlice';

export const store = configureStore({
    reducer: {
        chat: chatReducer,
        structure: structureReducer,
        mcp: mcpReducer,
    },
    middleware: (getDefaultMiddleware) =>
        getDefaultMiddleware({
            serializableCheck: {
                // Ignore these paths in the state for serializable check
                ignoredPaths: ['structure.structures'],
            },
        }),
});

// Expose Redux state to window for E2E testing
if (typeof window !== 'undefined') {
    Object.defineProperty(window, '__REDUX_STATE__', {
        get: () => store.getState(),
        configurable: true
    });
}

// Infer types from store
export type RootState = ReturnType<typeof store.getState>;
export type AppDispatch = typeof store.dispatch;
