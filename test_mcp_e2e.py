
import subprocess
import json
import sys
import time
import os

def run_e2e_test():
    print("=== MCP End-to-End Integration Test ===")
    
    # Path to server
    server_path = os.path.abspath("dist/index.js")
    if not os.path.exists(server_path):
        print("Error: Server dist/index.js not found.")
        return False

    # Start server process
    print(f"Starting server: node {server_path}")
    process = subprocess.Popen(
        ["node", server_path],
        stdin=subprocess.PIPE,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
        bufsize=0
    )

    def rpc_request(method, params=None, req_id=1):
        req = {
            "jsonrpc": "2.0",
            "method": method,
            "id": req_id
        }
        if params is not None:
            req["params"] = params
        
        json_req = json.dumps(req)
        print(f"-> Sending: {method} (id={req_id})")
        process.stdin.write(json_req + "\n")
        process.stdin.flush()
        
        # Read response (blocking read line)
        response_line = process.stdout.readline()
        print(f"<- Received: {len(response_line)} bytes: {response_line.strip()}")
        if not response_line:
            return None
        return json.loads(response_line)

    try:
        # 1. Initialize
        res = rpc_request("initialize", {
            "protocolVersion": "2024-11-05", # Dummy version
            "capabilities": {},
            "clientInfo": {"name": "test-client", "version": "1.0"}
        }, 1)
        
        if not res or "result" not in res:
            print("FAIL: Init failed")
            print(res)
            return False
        print("   Init OK. Server: " + res["result"]["serverInfo"]["name"])
        
        # 2. Tools List
        res = rpc_request("tools/list", {}, 2)
        tools = res["result"]["tools"]
        print(f"   Found {len(tools)} tools.")
        
        comp_tool = next((t for t in tools if t["name"] == "comprehensive_generate"), None)
        if not comp_tool:
            print("FAIL: comprehensive_generate tool not found")
            return False
        print("   Verified 'comprehensive_generate' exists.")
        
        # 3. Call Tool: Generate Structure (Si)
        print("   Calling comprehensive_generate (generate_from_spacegroup)...")
        res = rpc_request("tools/call", {
            "name": "comprehensive_generate",
            "arguments": {
                "operation": "generate_from_spacegroup",
                "spacegroup": 227,
                "elements": ["Si"],
                "composition": [8],
                "a": 5.43
            }
        }, 3)
        
        if "error" in res:
            print("FAIL: Generation error")
            print(res["error"])
            return False
            
        tool_res = res["result"]
        # MCP tool result is list of content items. Python returns JSON text?
        # My typescript handler returns content: [{ type: "text", text: JSON_STRING }]
        
        if len(tool_res["content"]) > 1:
            content_text = tool_res["content"][1]["text"]
        else:
            # Fallback if only one item (e.g. error or list)
            content_text = tool_res["content"][0]["text"]
            
        structure_data = json.loads(content_text)
        
        if not structure_data.get("success"):
            print("FAIL: Python generation failed")
            print(structure_data)
            return False
            
        print("   Generation OK. Got structure.")
        structure_dict = structure_data["structure"]
        
        # 4. Call Tool: Export Structure (VASP)
        print("   Calling comprehensive_generate (export_vasp)...")
        res = rpc_request("tools/call", {
            "name": "comprehensive_generate",
            "arguments": {
                "operation": "export_vasp",
                "structure": structure_dict
            }
        }, 4)
        
        tool_res = res["result"]
        if len(tool_res["content"]) > 1:
            content_text = tool_res["content"][1]["text"]
        else:
            content_text = tool_res["content"][0]["text"]
            
        try:
            export_data = json.loads(content_text)
        except json.JSONDecodeError:
            print(f"FAIL: Could not parse JSON from: {content_text[:100]}...")
            return False
        
        if not export_data.get("success"):
            print("FAIL: Export failed")
            print(export_data)
            return False
            
        print("   Export OK.")
        poscar = export_data["content"]
        if "Si" in poscar and "Direct" in poscar:
            print("   Verified POSCAR content.")
        else:
            print("WARNING: POV content suspicious.")
            print(poscar)

        # Cleanup
        process.terminate()
        print("\nSUCCESS: End-to-End Test Passed.")
        return True

    except Exception as e:
        print(f"\nCRASH: {e}")
        process.terminate()
        return False

if __name__ == "__main__":
    if run_e2e_test():
        sys.exit(0)
    else:
        sys.exit(1)
