#!/usr/bin/env python3
"""
Crystal MCP Server - Web Interface

A Flask-based web interface that:
1. Provides a chat-like UI for natural language input
2. Connects to the Crystal MCP Server via stdio
3. Displays generated structures with 3D visualization

This is a REAL interface - it runs the actual MCP server.
"""

import subprocess
import json
import sys
import os
import threading
import queue
import time
from pathlib import Path
from flask import Flask, render_template_string, request, jsonify

PROJECT_ROOT = Path(__file__).parent.parent

app = Flask(__name__)

# Global MCP client
mcp_client = None


class MCPWebClient:
    """MCP Client for web interface."""

    def __init__(self):
        self.process = None
        self.request_id = 0
        self.response_queue = queue.Queue()
        self.running = False
        self.initialized = False
        self.logs = []

    def start(self):
        """Start the MCP server."""
        server_path = PROJECT_ROOT / "dist" / "index.js"

        if not server_path.exists():
            raise FileNotFoundError(f"Server not built: {server_path}")

        env = os.environ.copy()
        env["PYTHON_PATH"] = sys.executable

        self.process = subprocess.Popen(
            ["node", str(server_path)],
            stdin=subprocess.PIPE,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            cwd=str(PROJECT_ROOT),
            env=env,
            text=True,
            bufsize=1
        )

        self.running = True

        # Start reader threads
        threading.Thread(target=self._read_stderr, daemon=True).start()
        threading.Thread(target=self._read_stdout, daemon=True).start()

        # Wait and initialize
        time.sleep(2)
        self._initialize()

        return True

    def _read_stderr(self):
        while self.running and self.process:
            line = self.process.stderr.readline()
            if line:
                self.logs.append(f"[SERVER] {line.rstrip()}")
            elif self.process.poll() is not None:
                break

    def _read_stdout(self):
        while self.running and self.process:
            line = self.process.stdout.readline()
            if line:
                try:
                    response = json.loads(line.strip())
                    self.response_queue.put(response)
                except json.JSONDecodeError:
                    pass
            elif self.process.poll() is not None:
                break

    def _send(self, method, params=None):
        self.request_id += 1
        request = {"jsonrpc": "2.0", "id": self.request_id, "method": method}
        if params:
            request["params"] = params

        self.process.stdin.write(json.dumps(request) + "\n")
        self.process.stdin.flush()

        try:
            return self.response_queue.get(timeout=120)
        except queue.Empty:
            return {"error": "Timeout"}

    def _initialize(self):
        response = self._send("initialize", {
            "protocolVersion": "2024-11-05",
            "capabilities": {"roots": {"listChanged": True}},
            "clientInfo": {"name": "crystal-web-client", "version": "1.0.0"}
        })

        if "result" in response:
            self.process.stdin.write(json.dumps({
                "jsonrpc": "2.0",
                "method": "notifications/initialized"
            }) + "\n")
            self.process.stdin.flush()
            self.initialized = True

        return response

    def call_tool(self, tool_name, arguments):
        """Call an MCP tool."""
        return self._send("tools/call", {"name": tool_name, "arguments": arguments})

    def stop(self):
        self.running = False
        if self.process:
            self.process.terminate()


def parse_natural_language(query: str) -> dict:
    """
    Parse natural language to MCP tool call.
    In production, this would use an LLM. Here we use pattern matching.
    """
    query_lower = query.lower()

    # Pattern matching for structure requests
    if "silicon" in query_lower and "diamond" in query_lower:
        return {
            "tool": "generate_crystal",
            "params": {"composition": ["Si"] * 8, "space_group": 227, "seed": 42},
            "description": "Silicon diamond cubic structure (Fd-3m)"
        }
    elif "nacl" in query_lower or "rocksalt" in query_lower or ("sodium" in query_lower and "chlor" in query_lower):
        return {
            "tool": "generate_crystal",
            "params": {"composition": ["Na"] * 4 + ["Cl"] * 4, "space_group": 225, "seed": 42},
            "description": "NaCl rocksalt structure (Fm-3m)"
        }
    elif "perovskite" in query_lower or "batio3" in query_lower:
        return {
            "tool": "generate_crystal",
            "params": {"composition": ["Ba", "Ti", "O", "O", "O"], "space_group": 221, "seed": 42},
            "description": "BaTiO3 perovskite structure (Pm-3m)"
        }
    elif "graphite" in query_lower:
        return {
            "tool": "generate_crystal",
            "params": {"composition": ["C"] * 4, "space_group": 194, "seed": 42},
            "description": "Graphite structure (P6_3/mmc)"
        }
    elif "hbn" in query_lower or "boron nitride" in query_lower:
        return {
            "tool": "generate_crystal",
            "params": {"composition": ["B", "B", "N", "N"], "space_group": 194, "seed": 42},
            "description": "Hexagonal boron nitride (P6_3/mmc)"
        }
    elif "iron" in query_lower and "bcc" in query_lower:
        return {
            "tool": "generate_crystal",
            "params": {"composition": ["Fe", "Fe"], "space_group": 229, "seed": 42},
            "description": "BCC Iron (Im-3m)"
        }
    elif "gold" in query_lower or "fcc" in query_lower:
        return {
            "tool": "generate_crystal",
            "params": {"composition": ["Au"] * 4, "space_group": 225, "seed": 42},
            "description": "FCC Gold (Fm-3m)"
        }
    elif "zno" in query_lower or "zinc oxide" in query_lower:
        return {
            "tool": "generate_crystal",
            "params": {"composition": ["Zn", "Zn", "O", "O"], "space_group": 186, "seed": 42},
            "description": "Wurtzite ZnO (P6_3mc)"
        }
    elif "tio2" in query_lower or "rutile" in query_lower:
        return {
            "tool": "generate_crystal",
            "params": {"composition": ["Ti", "Ti", "O", "O", "O", "O"], "space_group": 136, "seed": 42},
            "description": "TiO2 Rutile (P4_2/mnm)"
        }
    else:
        return None


def generate_xyz_string(structure: dict) -> str:
    """Generate XYZ string from structure."""
    atoms = structure.get("atoms", [])
    lines = [str(len(atoms)), "Generated by Crystal MCP Server"]
    for atom in atoms:
        x, y, z = atom.get("cartesian", [0, 0, 0])
        lines.append(f"{atom['element']}  {x:.6f}  {y:.6f}  {z:.6f}")
    return "\n".join(lines)


HTML_TEMPLATE = '''
<!DOCTYPE html>
<html>
<head>
    <meta charset="utf-8">
    <title>Crystal MCP Server - Web Interface</title>
    <script src="https://3Dmol.csb.pitt.edu/build/3Dmol-min.js"></script>
    <style>
        * { margin: 0; padding: 0; box-sizing: border-box; }
        body {
            font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', sans-serif;
            background: linear-gradient(135deg, #1a1a2e 0%, #16213e 100%);
            min-height: 100vh;
            color: #fff;
        }
        .header {
            background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
            padding: 20px;
            text-align: center;
            box-shadow: 0 4px 20px rgba(0,0,0,0.3);
        }
        .header h1 { font-size: 28px; margin-bottom: 5px; }
        .header p { opacity: 0.9; font-size: 14px; }
        .status {
            display: inline-block;
            padding: 4px 12px;
            border-radius: 20px;
            font-size: 12px;
            margin-top: 10px;
        }
        .status.connected { background: #4caf50; }
        .status.disconnected { background: #f44336; }

        .container {
            display: flex;
            height: calc(100vh - 100px);
            padding: 20px;
            gap: 20px;
        }

        .chat-panel {
            flex: 1;
            display: flex;
            flex-direction: column;
            background: rgba(255,255,255,0.05);
            border-radius: 16px;
            overflow: hidden;
        }

        .chat-messages {
            flex: 1;
            padding: 20px;
            overflow-y: auto;
        }

        .message {
            margin-bottom: 16px;
            animation: fadeIn 0.3s ease;
        }

        @keyframes fadeIn {
            from { opacity: 0; transform: translateY(10px); }
            to { opacity: 1; transform: translateY(0); }
        }

        .message.user { text-align: right; }

        .message-content {
            display: inline-block;
            padding: 12px 18px;
            border-radius: 18px;
            max-width: 80%;
        }

        .message.user .message-content {
            background: linear-gradient(135deg, #667eea, #764ba2);
        }

        .message.assistant .message-content {
            background: rgba(255,255,255,0.1);
            text-align: left;
        }

        .chat-input {
            display: flex;
            padding: 20px;
            background: rgba(0,0,0,0.2);
            gap: 10px;
        }

        .chat-input input {
            flex: 1;
            padding: 15px 20px;
            border: none;
            border-radius: 25px;
            font-size: 16px;
            background: rgba(255,255,255,0.1);
            color: #fff;
        }

        .chat-input input::placeholder { color: rgba(255,255,255,0.5); }
        .chat-input input:focus { outline: none; background: rgba(255,255,255,0.15); }

        .chat-input button {
            padding: 15px 30px;
            border: none;
            border-radius: 25px;
            background: linear-gradient(135deg, #667eea, #764ba2);
            color: #fff;
            font-size: 16px;
            cursor: pointer;
            transition: transform 0.2s, box-shadow 0.2s;
        }

        .chat-input button:hover {
            transform: translateY(-2px);
            box-shadow: 0 5px 20px rgba(102, 126, 234, 0.4);
        }

        .viewer-panel {
            width: 500px;
            background: #fff;
            border-radius: 16px;
            overflow: hidden;
            display: flex;
            flex-direction: column;
        }

        #viewer { flex: 1; min-height: 400px; }

        .structure-info {
            padding: 20px;
            background: #f5f5f5;
            color: #333;
        }

        .structure-info h3 { color: #667eea; margin-bottom: 10px; }

        .property {
            display: flex;
            justify-content: space-between;
            padding: 5px 0;
            border-bottom: 1px solid #eee;
        }

        .property-label { color: #666; }
        .property-value { font-weight: 600; color: #333; }

        .examples {
            padding: 10px 20px;
            background: rgba(0,0,0,0.1);
            font-size: 13px;
            color: rgba(255,255,255,0.7);
        }

        .example-tag {
            display: inline-block;
            background: rgba(255,255,255,0.1);
            padding: 4px 10px;
            border-radius: 12px;
            margin: 3px;
            cursor: pointer;
            transition: background 0.2s;
        }

        .example-tag:hover { background: rgba(255,255,255,0.2); }

        .loading {
            display: inline-block;
            width: 20px;
            height: 20px;
            border: 2px solid rgba(255,255,255,0.3);
            border-top-color: #fff;
            border-radius: 50%;
            animation: spin 1s linear infinite;
        }

        @keyframes spin { to { transform: rotate(360deg); } }
    </style>
</head>
<body>
    <div class="header">
        <h1>Crystal MCP Server</h1>
        <p>Generate crystal structures using natural language via MCP Protocol</p>
        <div id="status-indicator" class="status {{ 'connected' if connected else 'disconnected' }}">
            {{ 'Connected to MCP Server' if connected else 'Disconnected' }}
        </div>
    </div>

    <div class="container">
        <div class="chat-panel">
            <div class="examples">
                Try:
                <span class="example-tag" data-query="Generate silicon diamond structure">silicon diamond</span>
                <span class="example-tag" data-query="Create NaCl rocksalt">NaCl rocksalt</span>
                <span class="example-tag" data-query="BaTiO3 perovskite">BaTiO3 perovskite</span>
                <span class="example-tag" data-query="Graphite carbon">graphite</span>
                <span class="example-tag" data-query="Iron BCC magnetic">iron BCC</span>
            </div>
            <div class="chat-messages" id="messages"></div>
            <div class="chat-input">
                <input type="text" id="query" placeholder="Describe the crystal structure you want...">
                <button id="generate-btn">Generate</button>
            </div>
        </div>

        <div class="viewer-panel">
            <div id="viewer"></div>
            <div class="structure-info" id="info">
                <h3>Structure Properties</h3>
                <p style="color: #999;">Generate a structure to see its properties</p>
            </div>
        </div>
    </div>

    <script>
        let viewer = null;

        // Initialize viewer on load
        document.addEventListener('DOMContentLoaded', function() {
            const element = document.getElementById('viewer');
            viewer = $3Dmol.createViewer(element, {backgroundColor: 'white'});
            viewer.render();

            // Add welcome message
            addMessage('Welcome! I am the Crystal MCP Server. Describe the crystal structure you want to generate in natural language, and I will create it for you with 3D visualization.', false);

            // Setup event listeners
            document.querySelectorAll('.example-tag').forEach(function(tag) {
                tag.addEventListener('click', function() {
                    document.getElementById('query').value = this.getAttribute('data-query');
                });
            });

            document.getElementById('generate-btn').addEventListener('click', sendQuery);
            document.getElementById('query').addEventListener('keypress', function(e) {
                if (e.key === 'Enter') sendQuery();
            });
        });

        function addMessage(text, isUser) {
            const messages = document.getElementById('messages');
            const div = document.createElement('div');
            div.className = 'message ' + (isUser ? 'user' : 'assistant');

            const content = document.createElement('div');
            content.className = 'message-content';
            content.textContent = text;

            div.appendChild(content);
            messages.appendChild(div);
            messages.scrollTop = messages.scrollHeight;

            return div;
        }

        function addLoadingMessage() {
            const messages = document.getElementById('messages');
            const div = document.createElement('div');
            div.className = 'message assistant';
            div.id = 'loading-msg';

            const content = document.createElement('div');
            content.className = 'message-content';

            const spinner = document.createElement('span');
            spinner.className = 'loading';
            content.appendChild(spinner);

            const text = document.createTextNode(' Generating structure via MCP...');
            content.appendChild(text);

            div.appendChild(content);
            messages.appendChild(div);
            messages.scrollTop = messages.scrollHeight;
        }

        function removeLoadingMessage() {
            const loading = document.getElementById('loading-msg');
            if (loading) loading.remove();
        }

        function updateViewer(xyzData, structure) {
            viewer.removeAllModels();
            viewer.addModel(xyzData, "xyz");
            viewer.setStyle({}, {stick: {radius: 0.15}, sphere: {scale: 0.25}});
            viewer.addUnitCell(viewer.getModel());
            viewer.zoomTo();
            viewer.render();

            // Update info panel using safe DOM methods
            const info = document.getElementById('info');
            info.textContent = '';  // Clear safely

            const h3 = document.createElement('h3');
            h3.textContent = structure.metadata.formula;
            info.appendChild(h3);

            const properties = [
                ['Space Group', structure.space_group.symbol + ' (#' + structure.space_group.number + ')'],
                ['Crystal System', structure.space_group.crystal_system],
                ['Atoms', structure.metadata.natoms],
                ['Volume', structure.lattice.volume.toFixed(2) + ' A^3'],
                ['a', structure.lattice.a.toFixed(4) + ' A'],
                ['b', structure.lattice.b.toFixed(4) + ' A'],
                ['c', structure.lattice.c.toFixed(4) + ' A']
            ];

            properties.forEach(function(prop) {
                const div = document.createElement('div');
                div.className = 'property';

                const label = document.createElement('span');
                label.className = 'property-label';
                label.textContent = prop[0];

                const value = document.createElement('span');
                value.className = 'property-value';
                value.textContent = prop[1];

                div.appendChild(label);
                div.appendChild(value);
                info.appendChild(div);
            });
        }

        async function sendQuery() {
            const input = document.getElementById('query');
            const query = input.value.trim();
            if (!query) return;

            addMessage(query, true);
            input.value = '';
            addLoadingMessage();

            try {
                const response = await fetch('/generate', {
                    method: 'POST',
                    headers: {'Content-Type': 'application/json'},
                    body: JSON.stringify({query: query})
                });

                const data = await response.json();
                removeLoadingMessage();

                if (data.success) {
                    const msg = 'Generated ' + data.structure.metadata.formula +
                        ' | Space Group: ' + data.structure.space_group.symbol +
                        ' | ' + data.structure.metadata.natoms + ' atoms' +
                        ' [via MCP tool: ' + data.tool + ']';
                    addMessage(msg, false);
                    updateViewer(data.xyz, data.structure);
                } else {
                    addMessage('Error: ' + (data.error || 'Failed to generate structure'), false);
                }
            } catch (error) {
                removeLoadingMessage();
                addMessage('Error: ' + error.message, false);
            }
        }
    </script>
</body>
</html>
'''


@app.route('/')
def index():
    connected = mcp_client and mcp_client.initialized
    return render_template_string(HTML_TEMPLATE, connected=connected)


@app.route('/generate', methods=['POST'])
def generate():
    global mcp_client

    data = request.json
    query = data.get('query', '')

    # Parse natural language
    parsed = parse_natural_language(query)

    if not parsed:
        return jsonify({
            "success": False,
            "error": f"Could not understand: '{query}'. Try: 'silicon diamond', 'NaCl rocksalt', 'BaTiO3 perovskite'"
        })

    # Make sure client is connected
    if not mcp_client or not mcp_client.initialized:
        return jsonify({"success": False, "error": "MCP Server not connected"})

    # Call MCP tool
    try:
        response = mcp_client.call_tool(parsed["tool"], parsed["params"])

        if "result" in response and "content" in response["result"]:
            for item in response["result"]["content"]:
                if item.get("type") == "text":
                    result = json.loads(item["text"])
                    if result.get("success"):
                        structure = result["structure"]
                        xyz = generate_xyz_string(structure)
                        return jsonify({
                            "success": True,
                            "structure": structure,
                            "xyz": xyz,
                            "tool": parsed["tool"],
                            "description": parsed["description"]
                        })
                    else:
                        return jsonify({
                            "success": False,
                            "error": result.get("error", {}).get("message", "Generation failed")
                        })

        return jsonify({"success": False, "error": "Invalid response from MCP server"})

    except Exception as e:
        return jsonify({"success": False, "error": str(e)})


def main():
    global mcp_client

    print("="*60)
    print("   Crystal MCP Server - Web Interface")
    print("   REAL MCP Protocol Communication")
    print("="*60)

    # Start MCP client
    print("\nStarting MCP Server...")
    mcp_client = MCPWebClient()

    try:
        mcp_client.start()
        print("MCP Server connected and initialized")
    except FileNotFoundError as e:
        print(f"\nERROR: {e}")
        print("Please run 'npm run build' first!")
        return

    print("\n" + "="*60)
    print("   Web Interface Ready!")
    print("   Open http://localhost:5000 in your browser")
    print("="*60 + "\n")

    try:
        app.run(host='0.0.0.0', port=5000, debug=False)
    finally:
        if mcp_client:
            mcp_client.stop()


if __name__ == '__main__':
    main()
