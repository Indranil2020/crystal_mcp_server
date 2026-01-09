import { Component, type ErrorInfo, type ReactNode } from 'react';

// Debug logger
function debugLog(action: string, details?: string): void {
    const timestamp = new Date().toISOString().split('T')[1].slice(0, 12);
    console.log(`%c[ERROR_BOUNDARY] ${timestamp} ${action}${details ? ': ' + details : ''}`, 'color: #F44336; font-weight: bold');
}

interface Props {
    children: ReactNode;
    fallback?: ReactNode;
    onRetry?: () => void;
}

interface State {
    hasError: boolean;
    error?: Error;
    errorInfo?: ErrorInfo;
}

/**
 * Error Boundary for Viewer Components
 * 
 * Enhanced version with:
 * - Retry button that resets state and re-mounts children
 * - Specific recovery suggestions based on error type
 * - Error logging and clear error log button
 */
export class ViewerErrorBoundary extends Component<Props, State> {
    constructor(props: Props) {
        super(props);
        this.state = { hasError: false };
    }

    static getDerivedStateFromError(error: Error): State {
        debugLog('CAUGHT', error.message);
        return { hasError: true, error };
    }

    componentDidCatch(error: Error, errorInfo: ErrorInfo) {
        debugLog('ERROR_DETAILS', `${error.name}: ${error.message}`);
        debugLog('STACK', error.stack || 'No stack trace');
        debugLog('COMPONENT_STACK', errorInfo.componentStack || 'No component stack');

        this.setState({ errorInfo });
    }

    handleRetry = () => {
        debugLog('RETRY', 'User requested retry - resetting error state');
        this.setState({ hasError: false, error: undefined, errorInfo: undefined });
        if (this.props.onRetry) {
            this.props.onRetry();
        }
    };

    handleClearLog = () => {
        debugLog('CLEAR_LOG', 'User cleared error log');
        console.clear();
    };

    // Get specific recovery suggestions based on error type
    getRecoverySuggestions(): { title: string; suggestions: string[] } {
        const { error } = this.state;
        const errorMessage = error?.message?.toLowerCase() || '';
        const errorName = error?.name?.toLowerCase() || '';

        // WebGL-related errors
        if (errorMessage.includes('webgl') || errorMessage.includes('gl') || errorMessage.includes('context')) {
            return {
                title: 'WebGL Error',
                suggestions: [
                    'Check that GPU drivers are up to date',
                    'Try disabling hardware acceleration in browser settings',
                    'Refresh the page and try again',
                    'Use a different browser (Chrome/Firefox recommended)',
                ]
            };
        }

        // Memory-related errors
        if (errorMessage.includes('memory') || errorMessage.includes('out of') || errorName.includes('rangeerror')) {
            return {
                title: 'Memory Error',
                suggestions: [
                    'The structure may be too large - try a smaller molecule',
                    'Close other browser tabs to free memory',
                    'Refresh the page and try with fewer atoms',
                    'Switch to a lower quality rendering mode',
                ]
            };
        }

        // Network/MCP errors
        if (errorMessage.includes('network') || errorMessage.includes('fetch') || errorMessage.includes('mcp') || errorMessage.includes('connection')) {
            return {
                title: 'Connection Error',
                suggestions: [
                    'Check that the MCP Bridge is running (python bridge/server.py)',
                    'Verify your network connection',
                    'Restart the MCP server and refresh the page',
                    'Check browser console for detailed error messages',
                ]
            };
        }

        // Parsing/format errors
        if (errorMessage.includes('parse') || errorMessage.includes('format') || errorMessage.includes('invalid')) {
            return {
                title: 'Data Format Error',
                suggestions: [
                    'The structure data may be corrupted or invalid',
                    'Try generating the structure again',
                    'Check that the input format is correct (SMILES, CIF, etc.)',
                    'Report this issue with the exact input used',
                ]
            };
        }

        // Generic fallback
        return {
            title: 'Viewer Error',
            suggestions: [
                'Try refreshing the page',
                'Clear browser cache and reload',
                'Check the browser console for more details',
                'If the problem persists, try a different browser',
            ]
        };
    }

    render() {
        if (this.state.hasError) {
            if (this.props.fallback) {
                return this.props.fallback;
            }

            const recovery = this.getRecoverySuggestions();

            return (
                <div className="flex-1 flex items-center justify-center bg-slate-900 p-4">
                    <div className="max-w-md w-full bg-slate-800 rounded-xl p-6 shadow-xl border border-slate-700">
                        {/* Error Header */}
                        <div className="flex items-center gap-3 mb-4">
                            <div className="w-10 h-10 rounded-full bg-red-900/50 flex items-center justify-center">
                                <span className="text-red-400 text-xl">⚠️</span>
                            </div>
                            <div>
                                <h3 className="text-lg font-semibold text-amber-400">{recovery.title}</h3>
                                <p className="text-sm text-slate-400">Something went wrong in the viewer</p>
                            </div>
                        </div>

                        {/* Error Message */}
                        <div className="bg-slate-900 rounded-lg p-3 mb-4 font-mono text-xs text-red-300 overflow-x-auto">
                            {this.state.error?.message || 'Unknown error'}
                        </div>

                        {/* Recovery Suggestions */}
                        <div className="mb-4">
                            <p className="text-sm text-slate-300 font-medium mb-2">Try these steps:</p>
                            <ul className="text-sm text-slate-400 space-y-1">
                                {recovery.suggestions.map((suggestion, index) => (
                                    <li key={index} className="flex items-start gap-2">
                                        <span className="text-blue-400 mt-0.5">•</span>
                                        <span>{suggestion}</span>
                                    </li>
                                ))}
                            </ul>
                        </div>

                        {/* Action Buttons */}
                        <div className="flex gap-3">
                            <button
                                onClick={this.handleRetry}
                                className="flex-1 px-4 py-2 bg-blue-600 hover:bg-blue-500 text-white rounded-lg font-medium transition-colors"
                            >
                                🔄 Retry
                            </button>
                            <button
                                onClick={this.handleClearLog}
                                className="px-4 py-2 bg-slate-700 hover:bg-slate-600 text-slate-300 rounded-lg transition-colors"
                                title="Clear browser console"
                            >
                                Clear Log
                            </button>
                        </div>

                        {/* Console hint */}
                        <p className="text-xs text-slate-500 mt-4 text-center">
                            Check browser console (F12) for detailed error information
                        </p>
                    </div>
                </div>
            );
        }

        return this.props.children;
    }
}
