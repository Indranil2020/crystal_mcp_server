import { Component, type ErrorInfo, type ReactNode } from 'react';

interface Props {
    children: ReactNode;
    fallback?: ReactNode;
}

interface State {
    hasError: boolean;
    error?: Error;
}

/**
 * Error Boundary for Viewer Components
 * Prevents crashes in 3D/2D viewers from taking down entire app
 */
export class ViewerErrorBoundary extends Component<Props, State> {
    constructor(props: Props) {
        super(props);
        this.state = { hasError: false };
    }

    static getDerivedStateFromError(error: Error): State {
        return { hasError: true, error };
    }

    componentDidCatch(error: Error, errorInfo: ErrorInfo) {
        console.error('[ViewerErrorBoundary] Caught error:', error);
        console.error('[ViewerErrorBoundary] Error info:', errorInfo);
    }

    render() {
        if (this.state.hasError) {
            if (this.props.fallback) {
                return this.props.fallback;
            }

            return (
                <div className="flex-1 flex items-center justify-center bg-slate-900">
                    <div className="text-center text-amber-400 p-4">
                        <p className="text-lg">⚠️ Viewer Error</p>
                        <p className="text-sm mt-2 text-slate-400">
                            {this.state.error?.message || 'Component failed to load'}
                        </p>
                        <p className="text-xs mt-1 text-slate-500">
                            (Check console for details)
                        </p>
                    </div>
                </div>
            );
        }

        return this.props.children;
    }
}
