/**
 * Crystal GUI Web - Main App Component
 * 
 * Layout: 
 * - Top: Toolbar
 * - Left: Chat Panel with LLM + Tool visualization
 * - Center: 3D Viewer (MolStar) 
 * - Right: 2D Editor (Kekule.js)
 * - Bottom: Status Bar
 */

import { useEffect, useState, Suspense, lazy, useMemo } from 'react';
import { useAppDispatch } from './store/hooks';
import { setConnectionStatus, setTools, setError } from './store/mcpSlice';
import { mcpClient, llmClient, toolOrchestrator } from './services';
import ChatPanel from './components/ChatPanel/ChatPanel';
import Toolbar from './components/Toolbar/Toolbar';
import StatusBar from './components/StatusBar/StatusBar';

// Lazy load heavy viewer components (only if needed)
const MolStarViewer = lazy(() => import('./components/Viewers/MolStarViewer'));
const KekuleEditor = lazy(() => import('./components/Viewers/KekuleEditor'));

// Check WebGL support (runs once at module load)
function checkWebGLSupport(): boolean {
  try {
    const canvas = document.createElement('canvas');
    const gl = canvas.getContext('webgl') || canvas.getContext('experimental-webgl');
    return !!gl;
  } catch {
    return false;
  }
}

// Loading placeholder
function ViewerLoading({ label }: { label: string }) {
  return (
    <div className="flex-1 flex items-center justify-center bg-slate-900">
      <div className="text-center text-slate-400">
        <div className="w-8 h-8 border-2 border-blue-500 border-t-transparent rounded-full animate-spin mx-auto mb-2" />
        <p className="text-sm">Loading {label}...</p>
      </div>
    </div>
  );
}

// Fallback when WebGL not available
function NoWebGLFallback() {
  return (
    <div className="flex-1 flex items-center justify-center bg-slate-900">
      <div className="text-center text-amber-400 p-4">
        <p className="text-lg">⚠️ 3D Viewer Unavailable</p>
        <p className="text-sm mt-2 text-slate-400">WebGL is not supported in this browser.</p>
        <p className="text-xs mt-1 text-slate-500">(Headless browser or GPU disabled)</p>
      </div>
    </div>
  );
}

function App() {
  const dispatch = useAppDispatch();
  const [show2DEditor] = useState(true);

  // Check WebGL support once at mount
  const hasWebGL = useMemo(() => checkWebGLSupport(), []);

  // Initialize connections on mount
  useEffect(() => {
    async function initializeConnections() {
      dispatch(setConnectionStatus('connecting'));

      try {
        // Check LLM connection
        const llmConnected = await llmClient.checkConnection();
        console.log('[App] LLM connection:', llmConnected ? 'OK' : 'Failed');

        // Check MCP bridge connection
        const mcpConnected = await mcpClient.checkConnection();
        if (!mcpConnected) {
          dispatch(setError('MCP Bridge not available. Start: cd bridge && python server.py'));
          dispatch(setConnectionStatus('error'));
          return;
        }

        // Initialize MCP and fetch tools
        await mcpClient.initialize();
        const tools = await mcpClient.listTools();
        dispatch(setTools(tools));

        // Also initialize the tool orchestrator with the same tools
        toolOrchestrator.setTools(tools);

        dispatch(setConnectionStatus('connected'));

        console.log('[App] Initialized with', tools.length, 'tools');
        console.log('[App] WebGL support:', hasWebGL);
      } catch (error) {
        console.error('[App] Initialization error:', error);
        dispatch(setError(error instanceof Error ? error.message : 'Connection failed'));
        dispatch(setConnectionStatus('error'));
      }
    }

    initializeConnections();
  }, [dispatch, hasWebGL]);

  return (
    <div className="h-screen flex flex-col bg-slate-900">
      {/* Toolbar */}
      <Toolbar />

      {/* Main Content Area */}
      <div className="flex-1 flex min-h-0">
        {/* Left Panel: Chat */}
        <div className="w-96 flex-shrink-0 border-r border-slate-700 flex flex-col">
          <ChatPanel />
        </div>

        {/* Center: 3D Viewer (MolStar) - only load if WebGL available */}
        <div className="flex-1 flex flex-col min-w-0">
          {hasWebGL ? (
            <Suspense fallback={<ViewerLoading label="3D Viewer" />}>
              <MolStarViewer />
            </Suspense>
          ) : (
            <NoWebGLFallback />
          )}
        </div>

        {/* Right Panel: 2D Editor (Kekule.js) */}
        {show2DEditor && (
          <div className="w-96 flex-shrink-0 border-l border-slate-700 flex flex-col">
            <Suspense fallback={<ViewerLoading label="2D Editor" />}>
              <KekuleEditor />
            </Suspense>
          </div>
        )}
      </div>

      {/* Status Bar */}
      <StatusBar />
    </div>
  );
}

export default App;
