"""
MCP Bridge Server - FastAPI HTTP/WebSocket bridge to MCP Server

Bridges browser frontend to the MCP server which uses stdio communication.
Provides both HTTP endpoints (fallback) and WebSocket (streaming) support.
"""

import asyncio
import json
import subprocess
import sys
from pathlib import Path
from typing import Optional

from fastapi import FastAPI, HTTPException, WebSocket, WebSocketDisconnect
from fastapi.middleware.cors import CORSMiddleware
from pydantic import BaseModel

# Path to the MCP server
MCP_SERVER_PATH = Path(__file__).parent.parent.parent / "dist" / "index.js"

app = FastAPI(
    title="MCP Bridge Server",
    description="HTTP/WebSocket bridge for MCP Server communication",
    version="1.0.0",
)

# CORS for local development - allow all localhost ports (Vite may use different ports)
app.add_middleware(
    CORSMiddleware,
    allow_origin_regex=r"^http://(localhost|127\.0\.0\.1):\d+$",  # Any localhost port
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)


class McpProcess:
    """Manages a single MCP server subprocess."""
    
    def __init__(self):
        self.process: Optional[subprocess.Popen] = None
        self.request_id = 0
        self.initialized = False
        self._lock = asyncio.Lock()
    
    async def start(self):
        """Start the MCP server process."""
        if self.process is not None:
            return
        
        self.process = subprocess.Popen(
            ["node", str(MCP_SERVER_PATH)],
            stdin=subprocess.PIPE,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            bufsize=1,
        )
    
    async def stop(self):
        """Stop the MCP server process."""
        if self.process:
            self.process.terminate()
            self.process = None
            self.initialized = False
    
    def _next_id(self) -> int:
        self.request_id += 1
        return self.request_id
    
    async def send_request(self, method: str, params: dict = None) -> dict:
        """Send a JSON-RPC request and wait for response."""
        async with self._lock:
            if not self.process:
                await self.start()
            
            request = {
                "jsonrpc": "2.0",
                "id": self._next_id(),
                "method": method,
                "params": params or {},
            }
            
            request_str = json.dumps(request) + "\n"
            self.process.stdin.write(request_str)
            self.process.stdin.flush()
            
            # Read response
            response_line = self.process.stdout.readline()
            if not response_line:
                raise HTTPException(status_code=500, detail="Empty response from MCP server")
            
            response = json.loads(response_line)
            
            if "error" in response:
                raise HTTPException(
                    status_code=400,
                    detail=response["error"].get("message", "MCP error")
                )
            
            return response.get("result", {})
    
    async def initialize(self):
        """Initialize the MCP connection."""
        if self.initialized:
            return
        
        await self.send_request("initialize", {
            "protocolVersion": "2024-11-05",
            "capabilities": {},
            "clientInfo": {
                "name": "crystal-gui-bridge",
                "version": "1.0.0",
            },
        })
        
        # Send initialized notification
        if self.process:
            notification = json.dumps({
                "jsonrpc": "2.0",
                "method": "notifications/initialized",
            }) + "\n"
            self.process.stdin.write(notification)
            self.process.stdin.flush()
        
        self.initialized = True


# Global MCP process
mcp = McpProcess()


# Pydantic models
class InitializeRequest(BaseModel):
    protocolVersion: str = "2024-11-05"
    capabilities: dict = {}
    clientInfo: dict = {"name": "crystal-gui-web", "version": "1.0.0"}


class ToolCallRequest(BaseModel):
    name: str
    arguments: dict = {}


# Health check
@app.get("/health")
async def health_check():
    """Health check endpoint."""
    return {"status": "ok", "mcp_initialized": mcp.initialized}


# MCP endpoints
@app.post("/mcp/initialize")
async def mcp_initialize(request: InitializeRequest):
    """Initialize MCP connection."""
    await mcp.initialize()
    return {"success": True, "initialized": True}


@app.get("/mcp/tools")
async def mcp_list_tools():
    """List available MCP tools."""
    if not mcp.initialized:
        await mcp.initialize()
    
    result = await mcp.send_request("tools/list", {})
    return {"tools": result.get("tools", [])}


@app.post("/mcp/call")
async def mcp_call_tool(request: ToolCallRequest):
    """Call an MCP tool."""
    if not mcp.initialized:
        await mcp.initialize()
    
    result = await mcp.send_request("tools/call", {
        "name": request.name,
        "arguments": request.arguments,
    })
    
    return result


# WebSocket for streaming (future enhancement)
@app.websocket("/mcp/ws")
async def mcp_websocket(websocket: WebSocket):
    """WebSocket endpoint for streaming MCP communication."""
    await websocket.accept()
    
    try:
        while True:
            data = await websocket.receive_json()
            
            method = data.get("method")
            params = data.get("params", {})
            
            if method == "tools/list":
                if not mcp.initialized:
                    await mcp.initialize()
                result = await mcp.send_request("tools/list", {})
                await websocket.send_json({"type": "tools", "data": result})
            
            elif method == "tools/call":
                if not mcp.initialized:
                    await mcp.initialize()
                result = await mcp.send_request("tools/call", params)
                await websocket.send_json({"type": "result", "data": result})
            
            else:
                await websocket.send_json({"type": "error", "message": f"Unknown method: {method}"})
    
    except WebSocketDisconnect:
        pass


# Startup/shutdown events
@app.on_event("startup")
async def startup():
    """Start MCP process on server startup."""
    await mcp.start()


@app.on_event("shutdown")
async def shutdown():
    """Stop MCP process on server shutdown."""
    await mcp.stop()


if __name__ == "__main__":
    import uvicorn
    uvicorn.run(app, host="0.0.0.0", port=8080)
