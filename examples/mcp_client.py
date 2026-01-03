#!/usr/bin/env python3
"""
Real MCP Client for Crystal MCP Server

This client connects to the Crystal MCP Server via stdio using the actual
MCP (Model Context Protocol) JSON-RPC format.

This is NOT a shortcut - this sends real MCP protocol messages to the server.
"""

import subprocess
import json
import sys
import os
from pathlib import Path
from typing import Optional, Dict, Any
import threading
import queue

PROJECT_ROOT = Path(__file__).parent.parent


class MCPClient:
    """Real MCP Client that communicates with the server via stdio."""

    def __init__(self):
        self.process: Optional[subprocess.Popen] = None
        self.request_id = 0
        self.response_queue = queue.Queue()
        self.reader_thread: Optional[threading.Thread] = None
        self.running = False

    def start_server(self):
        """Start the MCP server as a subprocess."""
        server_path = PROJECT_ROOT / "dist" / "index.js"

        if not server_path.exists():
            raise FileNotFoundError(f"Server not built. Run 'npm run build' first. Missing: {server_path}")

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

        # Start reader thread for stderr (server logs)
        self.stderr_thread = threading.Thread(target=self._read_stderr, daemon=True)
        self.stderr_thread.start()

        # Start reader thread for stdout (MCP responses)
        self.reader_thread = threading.Thread(target=self._read_responses, daemon=True)
        self.reader_thread.start()

        print("✓ MCP Server started")

    def _read_stderr(self):
        """Read server logs from stderr."""
        while self.running and self.process:
            line = self.process.stderr.readline()
            if line:
                print(f"[SERVER] {line.rstrip()}")
            elif self.process.poll() is not None:
                break

    def _read_responses(self):
        """Read MCP responses from stdout."""
        while self.running and self.process:
            line = self.process.stdout.readline()
            if line:
                try:
                    response = json.loads(line.strip())
                    self.response_queue.put(response)
                except json.JSONDecodeError:
                    print(f"[RAW] {line.rstrip()}")
            elif self.process.poll() is not None:
                break

    def _send_request(self, method: str, params: Optional[Dict] = None) -> Dict[str, Any]:
        """Send an MCP JSON-RPC request and wait for response."""
        self.request_id += 1

        request = {
            "jsonrpc": "2.0",
            "id": self.request_id,
            "method": method,
        }
        if params:
            request["params"] = params

        request_json = json.dumps(request)
        print(f"\n>>> SENDING MCP REQUEST:")
        print(f"    {request_json[:200]}{'...' if len(request_json) > 200 else ''}")

        self.process.stdin.write(request_json + "\n")
        self.process.stdin.flush()

        # Wait for response
        try:
            response = self.response_queue.get(timeout=60)
            print(f"\n<<< RECEIVED MCP RESPONSE:")
            response_str = json.dumps(response)
            print(f"    {response_str[:500]}{'...' if len(response_str) > 500 else ''}")
            return response
        except queue.Empty:
            return {"error": "Timeout waiting for response"}

    def initialize(self):
        """Send MCP initialize request."""
        print("\n" + "="*60)
        print("STEP 1: MCP INITIALIZE")
        print("="*60)

        response = self._send_request("initialize", {
            "protocolVersion": "2024-11-05",
            "capabilities": {
                "roots": {"listChanged": True}
            },
            "clientInfo": {
                "name": "crystal-mcp-test-client",
                "version": "1.0.0"
            }
        })

        if "result" in response:
            print("\n✓ Server initialized successfully")
            print(f"  Server: {response['result'].get('serverInfo', {}).get('name', 'Unknown')}")
            print(f"  Version: {response['result'].get('serverInfo', {}).get('version', 'Unknown')}")

            # Send initialized notification
            self._send_notification("notifications/initialized", {})

        return response

    def _send_notification(self, method: str, params: Dict = None):
        """Send an MCP notification (no response expected)."""
        notification = {
            "jsonrpc": "2.0",
            "method": method,
        }
        if params:
            notification["params"] = params

        self.process.stdin.write(json.dumps(notification) + "\n")
        self.process.stdin.flush()

    def list_tools(self):
        """List available MCP tools."""
        print("\n" + "="*60)
        print("STEP 2: LIST AVAILABLE TOOLS")
        print("="*60)

        response = self._send_request("tools/list", {})

        if "result" in response and "tools" in response["result"]:
            tools = response["result"]["tools"]
            print(f"\n✓ Found {len(tools)} MCP tools:")
            for tool in tools[:10]:  # Show first 10
                print(f"  • {tool['name']}: {tool.get('description', '')[:60]}...")
            if len(tools) > 10:
                print(f"  ... and {len(tools) - 10} more")

        return response

    def call_tool(self, tool_name: str, arguments: Dict[str, Any]):
        """Call an MCP tool."""
        print("\n" + "="*60)
        print(f"STEP 3: CALL TOOL - {tool_name}")
        print("="*60)
        print(f"Arguments: {json.dumps(arguments, indent=2)}")

        response = self._send_request("tools/call", {
            "name": tool_name,
            "arguments": arguments
        })

        return response

    def stop(self):
        """Stop the MCP server."""
        self.running = False
        if self.process:
            self.process.terminate()
            self.process.wait()
            print("\n✓ MCP Server stopped")


def save_output_files(result: Dict, output_dir: Path, name: str):
    """Save the generated structure to files."""
    output_dir.mkdir(exist_ok=True)

    if "result" in result and "content" in result["result"]:
        content = result["result"]["content"]
        for item in content:
            if item.get("type") == "text":
                text = item.get("text", "")
                # Parse JSON from <json-data> tags (new format)
                if "<json-data>" in text:
                    start = text.find("<json-data>") + 11
                    end = text.find("</json-data>")
                    try:
                        data = json.loads(text[start:end].strip())
                        if data.get("success") and "structure" in data:
                            # Save JSON
                            json_file = output_dir / f"{name}.json"
                            with open(json_file, 'w') as f:
                                json.dump(data, f, indent=2)
                            print(f"  ✓ Saved: {json_file}")

                            # Save CIF file if available
                            if data.get("files", {}).get("cif"):
                                cif_file = output_dir / f"{name}.cif"
                                with open(cif_file, 'w') as f:
                                    f.write(data["files"]["cif"])
                                print(f"  ✓ Saved: {cif_file}")

                            # Save XYZ file if available
                            if data.get("files", {}).get("xyz"):
                                xyz_file = output_dir / f"{name}.xyz"
                                with open(xyz_file, 'w') as f:
                                    f.write(data["files"]["xyz"])
                                print(f"  ✓ Saved: {xyz_file}")

                            return data
                    except json.JSONDecodeError as e:
                        print(f"  ✗ JSON parse error: {e}")
                # Fallback: try direct JSON parse
                else:
                    try:
                        data = json.loads(text)
                        if data.get("success") and "structure" in data:
                            json_file = output_dir / f"{name}.json"
                            with open(json_file, 'w') as f:
                                json.dump(data, f, indent=2)
                            print(f"  ✓ Saved: {json_file}")
                            return data
                    except json.JSONDecodeError:
                        pass
    return None


def main():
    """Run the real MCP client test."""
    print("="*70)
    print("   REAL MCP CLIENT TEST - Crystal Structure Generator")
    print("   Using actual MCP Protocol over stdio")
    print("="*70)

    client = MCPClient()
    output_dir = Path(__file__).parent / "mcp_real_outputs"

    try:
        # Start the MCP server
        client.start_server()

        # Wait for server to initialize
        import time
        time.sleep(2)

        # Step 1: Initialize
        init_response = client.initialize()
        if "error" in init_response:
            print(f"ERROR: {init_response['error']}")
            return

        time.sleep(0.5)

        # Step 2: List tools
        tools_response = client.list_tools()

        time.sleep(0.5)

        # Step 3: Generate a crystal structure using REAL MCP call
        print("\n" + "="*60)
        print("GENERATING CRYSTAL STRUCTURE VIA MCP")
        print("="*60)

        # Silicon Diamond structure
        print("\n>>> Natural Language Intent: 'Generate silicon in diamond structure'")
        print(">>> Translated to MCP tool call: generate_crystal")

        result = client.call_tool("generate_crystal", {
            "composition": ["Si", "Si", "Si", "Si", "Si", "Si", "Si", "Si"],
            "space_group": 227,
            "seed": 42
        })

        if result:
            data = save_output_files(result, output_dir, "silicon_diamond_mcp")
            if data and data.get("success"):
                structure = data["structure"]
                print(f"\n✓ Structure generated via MCP:")
                print(f"  Formula: {structure['metadata']['formula']}")
                print(f"  Space Group: {structure['space_group']['symbol']} (#{structure['space_group']['number']})")
                print(f"  Atoms: {structure['metadata']['natoms']}")
                print(f"  Volume: {structure['lattice']['volume']:.2f} Å³")

        # Step 4: Generate another structure
        print("\n>>> Natural Language Intent: 'Create NaCl rocksalt for ionic study'")
        print(">>> Translated to MCP tool call: generate_crystal")

        result2 = client.call_tool("generate_crystal", {
            "composition": ["Na", "Na", "Na", "Na", "Cl", "Cl", "Cl", "Cl"],
            "space_group": 225,
            "seed": 42
        })

        if result2:
            data2 = save_output_files(result2, output_dir, "nacl_rocksalt_mcp")
            if data2 and data2.get("success"):
                structure = data2["structure"]
                print(f"\n✓ Structure generated via MCP:")
                print(f"  Formula: {structure['metadata']['formula']}")
                print(f"  Space Group: {structure['space_group']['symbol']} (#{structure['space_group']['number']})")

        print("\n" + "="*60)
        print("MCP CLIENT TEST COMPLETE")
        print("="*60)
        print(f"\nOutput files saved to: {output_dir}")

    except FileNotFoundError as e:
        print(f"\nERROR: {e}")
        print("Please run 'npm run build' first to compile the TypeScript server.")
    except Exception as e:
        print(f"\nERROR: {e}")
        import traceback
        traceback.print_exc()
    finally:
        client.stop()


if __name__ == "__main__":
    main()
