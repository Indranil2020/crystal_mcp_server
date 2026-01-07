/**
 * StatusBar - Bottom status bar showing connection and progress
 */

import { useAppSelector } from '../../store/hooks';

export default function StatusBar() {
    const { connectionStatus, tools, error } = useAppSelector(state => state.mcp);
    const { structures } = useAppSelector(state => state.structure);
    const { isProcessing, toolExecutions } = useAppSelector(state => state.chat);

    const runningTools = toolExecutions.filter(e => e.status === 'running');

    const statusColor = {
        disconnected: 'text-red-400',
        connecting: 'text-yellow-400',
        connected: 'text-green-400',
        error: 'text-red-400',
    }[connectionStatus];

    return (
        <div className="status-bar flex items-center justify-between text-xs">
            {/* Left: Connection status */}
            <div className="flex items-center gap-4">
                <div className={`flex items-center gap-1 ${statusColor}`}>
                    <span className="w-2 h-2 rounded-full bg-current" />
                    <span className="capitalize">{connectionStatus}</span>
                </div>

                {connectionStatus === 'connected' && (
                    <span className="text-slate-500">
                        {tools.length} tools available
                    </span>
                )}

                {error && (
                    <span className="text-red-400 truncate max-w-md">
                        ⚠ {error}
                    </span>
                )}
            </div>

            {/* Center: Progress */}
            <div className="flex items-center gap-4">
                {isProcessing && (
                    <span className="text-blue-400 animate-pulse">
                        Processing...
                    </span>
                )}

                {runningTools.length > 0 && (
                    <span className="text-yellow-400">
                        Running: {runningTools.map(t => t.toolName).join(', ')}
                    </span>
                )}
            </div>

            {/* Right: Stats */}
            <div className="flex items-center gap-4 text-slate-500">
                <span>{structures.length} structure{structures.length !== 1 ? 's' : ''} loaded</span>
                <span>Crystal GUI Web v1.0.0</span>
            </div>
        </div>
    );
}
