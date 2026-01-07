/**
 * ToolCallCard - Visualizes MCP tool execution with status
 */

import type { ToolExecution } from '../../types';

interface Props {
    execution: ToolExecution;
}

export default function ToolCallCard({ execution }: Props) {
    const statusColors = {
        pending: 'border-slate-500',
        running: 'border-yellow-500',
        success: 'border-green-500',
        error: 'border-red-500',
    };

    const statusIcons = {
        pending: '⏳',
        running: '⚙️',
        success: '✅',
        error: '❌',
    };

    const duration = execution.endTime
        ? ((execution.endTime - execution.startTime) / 1000).toFixed(1)
        : ((Date.now() - execution.startTime) / 1000).toFixed(1);

    return (
        <div
            className={`tool-call border-l-4 ${statusColors[execution.status]} 
                  bg-slate-800/50 rounded-r-lg p-3 my-2`}
        >
            {/* Header */}
            <div className="flex items-center justify-between mb-2">
                <div className="flex items-center gap-2">
                    <span className="text-lg">{statusIcons[execution.status]}</span>
                    <span className="font-mono text-sm text-blue-400">{execution.toolName}</span>
                </div>
                <span className="text-xs text-slate-500">{duration}s</span>
            </div>

            {/* Arguments */}
            <div className="text-xs text-slate-400 mb-2">
                <span className="text-slate-500">Args: </span>
                <code className="bg-slate-900/50 px-1 rounded">
                    {JSON.stringify(execution.arguments).slice(0, 80)}...
                </code>
            </div>

            {/* Result/Error */}
            {execution.result && (
                <div className="text-sm text-green-400 mt-2">
                    {execution.result}
                </div>
            )}
            {execution.error && (
                <div className="text-sm text-red-400 mt-2">
                    Error: {execution.error}
                </div>
            )}

            {/* Loading animation */}
            {execution.status === 'running' && (
                <div className="flex gap-1 mt-2">
                    <div className="w-2 h-2 bg-yellow-500 rounded-full animate-bounce" style={{ animationDelay: '0ms' }} />
                    <div className="w-2 h-2 bg-yellow-500 rounded-full animate-bounce" style={{ animationDelay: '150ms' }} />
                    <div className="w-2 h-2 bg-yellow-500 rounded-full animate-bounce" style={{ animationDelay: '300ms' }} />
                </div>
            )}
        </div>
    );
}
