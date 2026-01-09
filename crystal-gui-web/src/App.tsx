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

import { useEffect, useState, Suspense, lazy, useMemo, useCallback } from 'react';
import {
  Panel,
  Group as PanelGroup,
  Separator as PanelResizeHandle
} from 'react-resizable-panels';
import { useAppDispatch } from './store/hooks';
import { setConnectionStatus, setTools, setError } from './store/mcpSlice';
import { setLoading } from './store/structureSlice';
import { mcpClient, llmClient, toolOrchestrator } from './services';
import ChatPanel from './components/ChatPanel/ChatPanel';
import Toolbar from './components/Toolbar/Toolbar';
import StatusBar from './components/StatusBar/StatusBar';
import { ViewerErrorBoundary } from './components/ViewerErrorBoundary';
import { LoadingOverlay } from './components/LoadingOverlay';
import { useKeyboardShortcuts } from './hooks/useKeyboardShortcuts';


// Startup banner - runs immediately when module loads
console.log('%c==============================================================', 'color: #4CAF50; font-weight: bold');
console.log('%c  [CRYSTAL GUI WEB] Starting...', 'color: #4CAF50; font-size: 16px; font-weight: bold');
console.log('%c==============================================================', 'color: #4CAF50; font-weight: bold');
console.log('%c  Required services:', 'color: #888');
console.log('%c    1. MCP Bridge: python bridge/server.py (port 8080)', 'color: #FFB74D');
console.log('%c    2. Ollama:     ollama serve (port 11434)', 'color: #FFB74D');
console.log('%c==============================================================', 'color: #4CAF50; font-weight: bold');

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
        <p className="text-lg">[WARNING] 3D Viewer Unavailable</p>
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

  // Cancel loading handler
  const handleCancelLoading = useCallback(() => {
    console.log('%c[App] [LOADING] User cancelled loading operation', 'color: #FF9800');
    dispatch(setLoading(false));
  }, [dispatch]);

  // Global keyboard shortcuts
  useKeyboardShortcuts({
    enabled: true,
    // Note: onExportScreenshot, onResetCamera, onToggleMeasurementMode 
    // are handled by MolStarViewer which exposes these via refs or callbacks
  });

  // Initialize connections on mount
  useEffect(() => {
    async function initializeConnections() {
      console.log('%c[App] [INIT] Initializing connections...', 'color: #2196F3; font-weight: bold');
      dispatch(setConnectionStatus('connecting'));

      // Step 1: Check LLM connection
      console.log('%c[App] Step 1/3: Checking Ollama LLM...', 'color: #888');
      const llmConnected = await llmClient.checkConnection();
      if (llmConnected) {
        console.log('%c[App] [OK] Ollama LLM: Connected', 'color: #4CAF50');
      } else {
        console.log('%c[App] [WARN] Ollama LLM: Not available (chat will not work)', 'color: #FF9800');
      }

      // Step 2: Check MCP bridge connection
      console.log('%c[App] Step 2/3: Checking MCP Bridge at http://localhost:8080...', 'color: #888');
      const mcpConnected = await mcpClient.checkConnection();
      if (!mcpConnected) {
        console.log('%c[App] [ERROR] MCP Bridge: Connection REFUSED!', 'color: #F44336; font-weight: bold');
        console.log('%c[App] [INFO] To fix: Open a new terminal and run:', 'color: #FFB74D');
        console.log('%c        cd crystal-gui-web/bridge && python server.py', 'color: #FFB74D; font-weight: bold');
        dispatch(setError('MCP Bridge not available. Start: cd bridge && python server.py'));
        dispatch(setConnectionStatus('error'));
        return;
      }
      console.log('%c[App] [OK] MCP Bridge: Connected', 'color: #4CAF50');

      // Step 3: Initialize MCP and fetch tools
      console.log('%c[App] Step 3/3: Initializing MCP and fetching tools...', 'color: #888');
      await mcpClient.initialize();
      const tools = await mcpClient.listTools();
      dispatch(setTools(tools));

      // Also initialize the tool orchestrator with the same tools
      toolOrchestrator.setTools(tools);

      dispatch(setConnectionStatus('connected'));

      console.log('%c[App] [SUCCESS] Initialization complete!', 'color: #4CAF50; font-weight: bold');
      console.log(`%c[App] [INFO] Loaded ${tools.length} tools`, 'color: #4CAF50');
      console.log(`%c[App] [INFO] WebGL: ${hasWebGL ? 'Supported' : 'Not available'}`, hasWebGL ? 'color: #4CAF50' : 'color: #FF9800');
      console.log('%c==============================================================', 'color: #4CAF50; font-weight: bold');
    }

    initializeConnections();
  }, [dispatch, hasWebGL]);

  return (
    <div className="h-screen flex flex-col bg-slate-900">
      {/* Toolbar */}
      <Toolbar />

      {/* Main Content Area - Resizable Panels */}
      <PanelGroup orientation="horizontal" className="flex-1">
        {/* Left Panel: Chat */}
        <Panel
          defaultSize="25%"
          minSize="15%"
          maxSize="40%"
          className="flex flex-col overflow-hidden"
        >
          <ChatPanel />
        </Panel>

        {/* Resize Handle */}
        <PanelResizeHandle className="w-2 bg-slate-700 hover:bg-blue-500 transition-colors cursor-col-resize flex items-center justify-center">
          <div className="w-0.5 h-8 bg-slate-500 rounded-full" />
        </PanelResizeHandle>

        {/* Center: 3D Viewer (MolStar) - only load if WebGL available */}
        <Panel
          defaultSize="50%"
          minSize="30%"
          className="flex flex-col overflow-hidden"
        >
          {hasWebGL ? (
            <ViewerErrorBoundary>
              <Suspense fallback={<ViewerLoading label="3D Viewer" />}>
                <MolStarViewer />
              </Suspense>
            </ViewerErrorBoundary>
          ) : (
            <NoWebGLFallback />
          )}
        </Panel>

        {/* Resize Handle */}
        {show2DEditor && (
          <PanelResizeHandle className="w-2 bg-slate-700 hover:bg-blue-500 transition-colors cursor-col-resize flex items-center justify-center">
            <div className="w-0.5 h-8 bg-slate-500 rounded-full" />
          </PanelResizeHandle>
        )}

        {/* Right Panel: 2D Editor (Kekule.js) */}
        {show2DEditor && (
          <Panel
            defaultSize="25%"
            minSize="15%"
            maxSize="40%"
            className="flex flex-col overflow-hidden"
          >
            <ViewerErrorBoundary>
              <Suspense fallback={<ViewerLoading label="2D Editor" />}>
                <KekuleEditor />
              </Suspense>
            </ViewerErrorBoundary>
          </Panel>
        )}
      </PanelGroup>

      {/* Status Bar */}
      <StatusBar />

      {/* Global Loading Overlay */}
      <LoadingOverlay onCancel={handleCancelLoading} />
    </div>
  );
}

export default App;
