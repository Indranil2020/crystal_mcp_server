/**
 * LoadingOverlay - Full-screen loading indicator for long operations
 * 
 * Features:
 * - Animated chemistry-themed spinner
 * - Progress text with optional percentage
 * - Cancel button for interruptible operations
 */

import { useAppSelector } from '../store/hooks';

interface LoadingOverlayProps {
    onCancel?: () => void;
}

function debugLog(action: string, details?: string): void {
    const timestamp = new Date().toISOString().split('T')[1].slice(0, 12);
    console.log(`%c[LOADING] ${timestamp} ${action}${details ? ': ' + details : ''}`, 'color: #00BCD4; font-weight: bold');
}

export function LoadingOverlay({ onCancel }: LoadingOverlayProps) {
    const { isLoading, loadingMessage, loadingProgress } = useAppSelector(state => state.structure);

    if (!isLoading) {
        return null;
    }

    debugLog('RENDER', `message="${loadingMessage}" progress=${loadingProgress}`);

    return (
        <div className="fixed inset-0 z-50 flex items-center justify-center bg-slate-900/80 backdrop-blur-sm">
            <div className="bg-slate-800 rounded-xl p-8 shadow-2xl border border-slate-700 min-w-[300px] text-center">
                {/* Chemistry-themed spinner - rotating molecule icon */}
                <div className="relative w-16 h-16 mx-auto mb-6">
                    {/* Outer ring */}
                    <div className="absolute inset-0 border-4 border-blue-500/30 rounded-full" />
                    {/* Spinning ring */}
                    <div className="absolute inset-0 border-4 border-transparent border-t-blue-500 rounded-full animate-spin" />
                    {/* Inner atom representation */}
                    <div className="absolute inset-3 flex items-center justify-center">
                        <div className="w-4 h-4 bg-blue-500 rounded-full animate-pulse" />
                    </div>
                    {/* Orbiting electrons */}
                    <div className="absolute inset-0 animate-spin" style={{ animationDuration: '2s' }}>
                        <div className="absolute top-0 left-1/2 -translate-x-1/2 w-2 h-2 bg-cyan-400 rounded-full" />
                    </div>
                    <div className="absolute inset-0 animate-spin" style={{ animationDuration: '3s', animationDirection: 'reverse' }}>
                        <div className="absolute bottom-0 left-1/2 -translate-x-1/2 w-2 h-2 bg-purple-400 rounded-full" />
                    </div>
                </div>

                {/* Loading text */}
                <p className="text-lg text-slate-200 font-medium mb-2">
                    {loadingMessage || 'Loading...'}
                </p>

                {/* Progress bar (if progress is available) */}
                {loadingProgress !== null && loadingProgress >= 0 && (
                    <div className="w-full bg-slate-700 rounded-full h-2 mb-4 overflow-hidden">
                        <div
                            className="bg-gradient-to-r from-blue-500 to-cyan-400 h-full rounded-full transition-all duration-300"
                            style={{ width: `${Math.min(100, Math.max(0, loadingProgress))}%` }}
                        />
                    </div>
                )}

                {/* Progress percentage */}
                {loadingProgress !== null && loadingProgress >= 0 && (
                    <p className="text-sm text-slate-400 mb-4">
                        {Math.round(loadingProgress)}% complete
                    </p>
                )}

                {/* Cancel button (if cancellation is supported) */}
                {onCancel && (
                    <button
                        onClick={() => {
                            debugLog('CANCEL', 'User requested cancellation');
                            onCancel();
                        }}
                        className="px-4 py-2 text-sm bg-slate-700 hover:bg-slate-600 text-slate-300 rounded-lg transition-colors"
                    >
                        Cancel
                    </button>
                )}
            </div>
        </div>
    );
}

export default LoadingOverlay;
