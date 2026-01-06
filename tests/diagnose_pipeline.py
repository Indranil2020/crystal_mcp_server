#!/usr/bin/env python3
"""
diagnose_pipeline.py - Systematic Pipeline Diagnosis

This script probes each component of the NL → Molecule Visualization pipeline
to identify exactly where failures occur.

Components tested:
1. Backend (Python): molecule_generator.py direct call
2. MCP Server: JSON-RPC tool call
3. LLM (Ollama): Tool calling capability
4. Full Integration: NL → LLM → MCP → Structure

Usage:
    python tests/diagnose_pipeline.py
"""

import subprocess
import json
import requests
import sys
import os
import tempfile
from pathlib import Path

# Configuration
MCP_SERVER_DIR = Path(__file__).parent.parent.absolute()
OLLAMA_URL = "http://localhost:11434"
OLLAMA_MODEL = "qwen2.5:7b"

# Colors for output
class Colors:
    GREEN = '\033[92m'
    RED = '\033[91m'
    YELLOW = '\033[93m'
    BLUE = '\033[94m'
    BOLD = '\033[1m'
    END = '\033[0m'

def print_header(text):
    print(f"\n{Colors.BOLD}{Colors.BLUE}{'='*70}")
    print(f"  {text}")
    print(f"{'='*70}{Colors.END}\n")

def print_success(text):
    print(f"{Colors.GREEN}✓ {text}{Colors.END}")

def print_failure(text):
    print(f"{Colors.RED}✗ {text}{Colors.END}")

def print_info(text):
    print(f"{Colors.YELLOW}→ {text}{Colors.END}")


# ============================================================================
# PROBE 1: Direct Python Backend Test
# ============================================================================
def probe_python_backend():
    """Test if the Python molecule_generator.py works directly."""
    print_header("PROBE 1: Python Backend (molecule_generator.py)")
    
    # Create temp input file
    input_data = {"name": "benzene", "input_type": "auto", "optimize": True, "vacuum": 10.0}
    
    with tempfile.NamedTemporaryFile(mode='w', suffix='.json', delete=False) as f:
        json.dump(input_data, f)
        input_file = f.name
    
    try:
        python_script = MCP_SERVER_DIR / "src" / "python" / "molecule_generator.py"
        result = subprocess.run(
            ["python3", str(python_script), input_file],
            capture_output=True,
            text=True,
            timeout=30
        )
        
        print_info(f"Command: python3 molecule_generator.py {input_file}")
        print_info(f"Input: {json.dumps(input_data)}")
        
        if result.returncode != 0:
            print_failure(f"Process failed with return code {result.returncode}")
            print(f"STDERR: {result.stderr[:500]}")
            return False, None
        
        try:
            output = json.loads(result.stdout)
        except json.JSONDecodeError as e:
            print_failure(f"Invalid JSON output: {e}")
            print(f"STDOUT: {result.stdout[:500]}")
            return False, None
        
        if output.get("success"):
            structure = output.get("structure", {})
            atoms = structure.get("atoms", structure.get("sites", []))
            print_success(f"Backend works! Generated {len(atoms)} atoms")
            print_info(f"Formula: {structure.get('metadata', {}).get('formula', 'Unknown')}")
            
            # Validate structure has required fields for GUI
            lattice = structure.get("lattice", {})
            space_group = structure.get("space_group", {})
            
            issues = []
            if not lattice.get("a"):
                issues.append("Missing lattice.a")
            if not space_group.get("number"):
                issues.append("Missing space_group.number")
            if not atoms:
                issues.append("No atoms/sites array")
            else:
                first_atom = atoms[0]
                if "cartesian" not in first_atom or first_atom["cartesian"] == [0, 0, 0]:
                    issues.append("Atoms missing cartesian coordinates")
            
            if issues:
                print_failure(f"Structure issues: {', '.join(issues)}")
                return False, output
            
            return True, output
        else:
            print_failure(f"Backend returned error: {output.get('error')}")
            return False, output
            
    finally:
        os.unlink(input_file)


# ============================================================================
# PROBE 2: MCP Server Tool Call
# ============================================================================
def probe_mcp_server():
    """Test if the MCP server correctly handles build_molecule tool."""
    print_header("PROBE 2: MCP Server (tools/call)")
    
    init_msg = json.dumps({
        "jsonrpc": "2.0",
        "id": 0,
        "method": "initialize",
        "params": {
            "protocolVersion": "2024-11-05",
            "capabilities": {},
            "clientInfo": {"name": "diagnose", "version": "1.0"}
        }
    })
    
    call_msg = json.dumps({
        "jsonrpc": "2.0",
        "id": 1,
        "method": "tools/call",
        "params": {
            "name": "build_molecule",
            "arguments": {"name": "benzene"}
        }
    })
    
    script = f"""cd {MCP_SERVER_DIR}
cat <<'MCPINPUT' | timeout 60 node dist/index.js 2>&1
{init_msg}
{call_msg}
MCPINPUT"""
    
    print_info("Sending MCP request for 'build_molecule' with name='benzene'")
    
    try:
        result = subprocess.run(
            ["bash", "-c", script],
            capture_output=True,
            text=True,
            timeout=120
        )
        
        output = result.stdout
        
        # Find the response with id=1
        for line in output.split('\n'):
            line = line.strip()
            if not line:
                continue
            
            if '"id":1' in line or '"id": 1' in line:
                try:
                    rpc = json.loads(line)
                    
                    if "error" in rpc:
                        print_failure(f"MCP Error: {rpc['error']}")
                        return False, None
                    
                    rpc_result = rpc.get('result', {})
                    is_error = rpc_result.get('isError', False)
                    content = rpc_result.get('content', [])
                    
                    if is_error:
                        error_text = content[0].get('text', '') if content else 'Unknown error'
                        print_failure(f"Tool returned error: {error_text}")
                        return False, None
                    
                    # Look for <json-data> in content
                    for item in content:
                        text = item.get('text', '')
                        if '<json-data>' in text:
                            start = text.find('<json-data>') + len('<json-data>')
                            end = text.find('</json-data>')
                            if end > start:
                                json_str = text[start:end].strip()
                                wrapper = json.loads(json_str)
                                structure = wrapper.get('structure', wrapper)
                                atoms = structure.get('atoms', structure.get('sites', []))
                                
                                print_success(f"MCP Server works! Got {len(atoms)} atoms")
                                
                                # Check structure format
                                print_info(f"Structure keys: {list(structure.keys())}")
                                print_info(f"Lattice keys: {list(structure.get('lattice', {}).keys())}")
                                print_info(f"Space group: {structure.get('space_group', {})}")
                                
                                return True, structure
                    
                    # No <json-data> found
                    print_failure("No <json-data> tag in response")
                    for item in content:
                        print_info(f"Content: {item.get('text', '')[:200]}")
                    return False, None
                    
                except json.JSONDecodeError:
                    continue
        
        print_failure("No valid response found")
        print_info(f"Raw output (first 500 chars): {output[:500]}")
        return False, None
        
    except subprocess.TimeoutExpired:
        print_failure("MCP Server timeout")
        return False, None
    except Exception as e:
        print_failure(f"Error: {e}")
        return False, None


# ============================================================================
# PROBE 3: Ollama Connection and Model
# ============================================================================
def probe_ollama_connection():
    """Test if Ollama is running and the model is available."""
    print_header("PROBE 3: Ollama Connection")
    
    try:
        response = requests.get(f"{OLLAMA_URL}/api/tags", timeout=10)
        if response.status_code != 200:
            print_failure(f"Ollama not responding: HTTP {response.status_code}")
            return False, None
        
        models = response.json().get("models", [])
        model_names = [m.get("name", "") for m in models]
        
        print_success(f"Ollama is running with {len(models)} models")
        
        # Check if our model is available
        if any(OLLAMA_MODEL in name for name in model_names):
            print_success(f"Model '{OLLAMA_MODEL}' is available")
        else:
            print_failure(f"Model '{OLLAMA_MODEL}' not found")
            print_info(f"Available models: {model_names[:5]}")
            return False, model_names
        
        return True, model_names
        
    except requests.exceptions.ConnectionError:
        print_failure("Cannot connect to Ollama. Is it running?")
        return False, None
    except Exception as e:
        print_failure(f"Error: {e}")
        return False, None


# ============================================================================
# PROBE 4: LLM Native Tool Calling
# ============================================================================
def probe_llm_tool_calling():
    """Test if Ollama can correctly generate tool calls."""
    print_header("PROBE 4: LLM Native Tool Calling")
    
    # Simplified tool definition matching what the GUI sends
    tools = [
        {
            "type": "function",
            "function": {
                "name": "build_molecule",
                "description": "Generate a single molecule from any identifier",
                "parameters": {
                    "type": "object",
                    "properties": {
                        "name": {
                            "type": "string",
                            "description": "Molecule identifier (name, SMILES, IUPAC, CID)"
                        },
                        "input_type": {
                            "type": "string",
                            "enum": ["auto", "name", "smiles", "iupac", "cid"],
                            "default": "auto"
                        }
                    },
                    "required": ["name"]
                }
            }
        }
    ]
    
    # System prompt similar to GUI's llm_client.rs
    system_prompt = """You are a tool-calling assistant for crystal and molecule generation.

CRITICAL TOOL SELECTION RULES:
1. For MOLECULES (water, aspirin, benzene, caffeine, PTCDA, H2O, CO2, ANY chemical compound):
   → Use: build_molecule with {"name": "<molecule_name>"}
   
EXAMPLES:
- "create water molecule" → build_molecule {"name": "water"}
- "generate H2O" → build_molecule {"name": "H2O"}
- "show benzene" → build_molecule {"name": "benzene"}

ALWAYS use build_molecule for any molecule request."""
    
    messages = [
        {"role": "system", "content": system_prompt},
        {"role": "user", "content": "Generate a benzene molecule"}
    ]
    
    payload = {
        "model": OLLAMA_MODEL,
        "messages": messages,
        "stream": False,
        "tools": tools
    }
    
    print_info("Sending to Ollama with native tool calling...")
    print_info(f"User message: 'Generate a benzene molecule'")
    
    try:
        response = requests.post(f"{OLLAMA_URL}/api/chat", json=payload, timeout=60)
        
        if response.status_code != 200:
            print_failure(f"Ollama error: HTTP {response.status_code}")
            return False, None
        
        result = response.json()
        message = result.get("message", {})
        
        # Check for tool calls
        tool_calls = message.get("tool_calls", [])
        content = message.get("content", "")
        
        print_info(f"Response content: {content[:200] if content else '(empty)'}")
        print_info(f"Tool calls: {len(tool_calls)}")
        
        if tool_calls:
            for tc in tool_calls:
                func = tc.get("function", {})
                name = func.get("name", "")
                args = func.get("arguments", {})
                print_success(f"LLM generated tool call: {name}")
                print_info(f"Arguments: {json.dumps(args)}")
                
                if name == "build_molecule" and args.get("name"):
                    print_success("Tool call is correct for molecule generation!")
                    return True, tc
                else:
                    print_failure(f"Unexpected tool call: {name}")
                    return False, tc
        else:
            print_failure("LLM did NOT generate a tool call!")
            print_info("This is the root cause: the LLM is returning text instead of tool calls")
            
            # Try to extract JSON from text as fallback
            if content:
                import re
                json_match = re.search(r'\{.*\}', content, re.DOTALL)
                if json_match:
                    try:
                        parsed = json.loads(json_match.group(0))
                        print_info(f"Found JSON in text response: {parsed}")
                        if parsed.get("tool") or parsed.get("name"):
                            print_info("LLM is outputting JSON text, but NOT using native tool calling")
                    except:
                        pass
            
            return False, content
            
    except Exception as e:
        print_failure(f"Error: {e}")
        return False, None


# ============================================================================
# PROBE 5: Full Integration (GUI Simulation)
# ============================================================================
def probe_full_integration():
    """Simulate the exact flow the GUI uses."""
    print_header("PROBE 5: Full Integration (GUI Flow Simulation)")
    
    # Step 1: Get tool list from MCP
    print_info("Step 1: Getting tools from MCP server...")
    
    init_msg = json.dumps({
        "jsonrpc": "2.0",
        "id": 0,
        "method": "initialize",
        "params": {
            "protocolVersion": "2024-11-05",
            "capabilities": {},
            "clientInfo": {"name": "gui-sim", "version": "1.0"}
        }
    })
    
    list_msg = json.dumps({
        "jsonrpc": "2.0",
        "id": 1,
        "method": "tools/list",
        "params": {}
    })
    
    script = f"""cd {MCP_SERVER_DIR}
cat <<'MCPIN' | timeout 30 node dist/index.js 2>&1
{init_msg}
{list_msg}
MCPIN"""
    
    try:
        result = subprocess.run(["bash", "-c", script], capture_output=True, text=True, timeout=60)
        
        # Parse tools list  
        tools = []
        for line in result.stdout.split('\n'):
            if '"id":1' in line or '"id": 1' in line:
                try:
                    rpc = json.loads(line)
                    tools_list = rpc.get('result', {}).get('tools', [])
                    tools = tools_list
                    break
                except:
                    continue
        
        if not tools:
            print_failure("Failed to get tools list from MCP")
            return False
        
        print_success(f"Got {len(tools)} tools from MCP")
        
        # Find build_molecule
        build_molecule_tool = None
        for t in tools:
            if t.get("name") == "build_molecule":
                build_molecule_tool = t
                break
        
        if not build_molecule_tool:
            print_failure("build_molecule tool not found!")
            return False
        
        print_success("Found build_molecule tool")
        
        # Step 2: Convert to Ollama format (like llm_client.rs does)
        print_info("Step 2: Converting tools to Ollama format...")
        
        ollama_tools = []
        for tool in tools[:5]:  # Limit to first 5 to reduce prompt size
            ollama_tools.append({
                "type": "function",
                "function": {
                    "name": tool.get("name"),
                    "description": tool.get("description", "")[:200],
                    "parameters": tool.get("inputSchema", {"type": "object", "properties": {}})
                }
            })
        
        # Step 3: Send to Ollama
        print_info("Step 3: Sending to Ollama with tools...")
        
        payload = {
            "model": OLLAMA_MODEL,
            "messages": [
                {"role": "system", "content": "You are a tool-calling assistant for crystal and molecule generation. For molecules, use build_molecule."},
                {"role": "user", "content": "Generate a water molecule"}
            ],
            "stream": False,
            "tools": ollama_tools
        }
        
        response = requests.post(f"{OLLAMA_URL}/api/chat", json=payload, timeout=60)
        result = response.json()
        message = result.get("message", {})
        tool_calls = message.get("tool_calls", [])
        
        if not tool_calls:
            print_failure("Ollama did not generate tool calls with MCP tools!")
            print_info(f"Response: {message.get('content', '')[:300]}")
            return False
        
        print_success(f"Ollama generated {len(tool_calls)} tool call(s)")
        
        # Check the tool call
        first_call = tool_calls[0]
        tool_name = first_call.get("function", {}).get("name", "")
        tool_args = first_call.get("function", {}).get("arguments", {})
        
        print_info(f"Tool: {tool_name}")
        print_info(f"Args: {json.dumps(tool_args)}")
        
        if tool_name != "build_molecule":
            print_failure(f"Wrong tool selected: {tool_name}")
            return False
        
        print_success("Correct tool selected!")
        return True
        
    except Exception as e:
        print_failure(f"Error: {e}")
        return False


# ============================================================================
# Main Diagnostic Flow
# ============================================================================
def main():
    print(f"\n{Colors.BOLD}{'#'*70}")
    print("#  CRYSTAL-MCP PIPELINE DIAGNOSTIC TOOL")
    print("#  Probing: NL Input → LLM → MCP Server → Molecule Visualization")
    print(f"{'#'*70}{Colors.END}\n")
    
    results = {}
    
    # Probe 1: Python Backend
    results["backend"] = probe_python_backend()[0]
    
    # Probe 2: MCP Server
    results["mcp_server"] = probe_mcp_server()[0]
    
    # Probe 3: Ollama Connection
    results["ollama"] = probe_ollama_connection()[0]
    
    # Probe 4: LLM Tool Calling
    if results["ollama"]:
        results["llm_tools"] = probe_llm_tool_calling()[0]
    else:
        results["llm_tools"] = False
    
    # Probe 5: Full Integration
    if all([results["backend"], results["mcp_server"], results["ollama"]]):
        results["integration"] = probe_full_integration()
    else:
        results["integration"] = False
    
    # Summary
    print_header("DIAGNOSTIC SUMMARY")
    
    all_pass = True
    for name, passed in results.items():
        if passed:
            print_success(f"{name}: PASSED")
        else:
            print_failure(f"{name}: FAILED")
            all_pass = False
    
    print("\n" + "-"*70)
    
    if all_pass:
        print_success("All probes passed! The pipeline should be working.")
    else:
        print("\n" + Colors.BOLD + "DIAGNOSIS:" + Colors.END)
        
        if not results["backend"]:
            print("  → Python molecule_generator.py is broken. Fix the backend first.")
        elif not results["mcp_server"]:
            print("  → MCP Server is not working correctly. Check dist/index.js compilation.")
        elif not results["ollama"]:
            print("  → Ollama is not running. Start it with 'ollama serve'.")
        elif not results["llm_tools"]:
            print("  → LLM is not generating tool calls. This is likely the root cause!")
            print("    The model may not support native tool calling, or the prompt needs adjustment.")
            print("    Consider:")
            print("      1. Using a model that supports function calling (llama3.1, qwen2.5)")
            print("      2. Adjusting the system prompt to be more explicit")
            print("      3. Falling back to JSON text parsing (already implemented in app.rs)")
        elif not results["integration"]:
            print("  → Full integration fails. Check tool schema compatibility.")
    
    return 0 if all_pass else 1


if __name__ == "__main__":
    sys.exit(main())
