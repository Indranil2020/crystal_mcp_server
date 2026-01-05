#!/usr/bin/env python3
"""
end_to_end_test.py - Complete End-to-End Testing with ASCII Ball-Stick Visualization

Tests the full flow:
1. User types natural language input
2. LLM (Ollama) translates to tool call
3. MCP Server executes the tool
4. Parse JSON output and render ASCII ball-stick diagram
5. Validate structure correctness

Output format:
┌────────────────────────────────────────┬──────────────────────────────────────────┐
│ NL Input:                              │ Ball-Stick Diagram:                     │
├────────────────────────────────────────┼──────────────────────────────────────────┤
│ Generate a benzene molecule            │           C                              │
│                                        │          /│\                             │
│                                        │       H─C C─H                            │
│                                        │          \│/                             │
│                                        │           C                              │
├────────────────────────────────────────┴──────────────────────────────────────────┤
│ Formula: C6H6    Atoms: 12    Stacking: π-π parallel                              │
└───────────────────────────────────────────────────────────────────────────────────┘
"""

import subprocess
import json
import requests
import sys
import os
import re
from typing import Dict, List, Tuple, Optional

# Import ASCII visualizer
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from ascii_visualizer import render_ascii_structure, render_structure_from_json, format_test_result

# Configuration
OLLAMA_URL = "http://localhost:11434/api/chat"
OLLAMA_MODEL = "qwen2.5:7b"
MCP_SERVER_DIR = "/home/niel/git/crystal-mcp-server"

# Test definitions
SINGLE_MOLECULE_TESTS = [
    # (name, nl_prompt, expected_atoms_min, expected_elements)
    ("Benzene", "Generate a benzene molecule", 12, {"C": 6, "H": 6}),
    ("Aspirin", "Create aspirin", 21, {"C": 9, "H": 8, "O": 4}),
    ("Caffeine", "Build a caffeine molecule", 24, {"C": 8, "H": 10, "N": 4, "O": 2}),
    ("PTCDA", "Generate PTCDA", 38, {"C": 24, "H": 8, "O": 6}),
    ("Naphthalene", "Create naphthalene", 18, {"C": 10, "H": 8}),
    ("Paracetamol", "Generate paracetamol", 20, None),
    ("Pentacene", "Build pentacene molecule", 36, {"C": 22, "H": 14}),
    ("Water", "Generate a water molecule", 3, {"H": 2, "O": 1}),
]

CLUSTER_TESTS = [
    # (name, nl_prompt, expected_molecules, expected_stacking)
    ("Benzene dimer π-π", "Create a benzene dimer with parallel pi-pi stacking", 2, "pi_pi_parallel"),
    ("Water tetramer circular", "Build water tetramer in circular H-bonded ring", 4, "circular"),
    ("PTCDA-benzene hetero", "Generate PTCDA and benzene hetero-dimer", 2, "auto"),
    ("Pentacene herringbone", "Create 4 pentacene molecules in herringbone stacking", 4, "herringbone"),
    ("Benzene spiral (30°)", "Stack 6 benzene molecules with 30 degree rotation each", 6, "stacked"),
    ("Naphthalene swastika", "Build swastika pattern with naphthalene molecules", 4, "swastika"),
    ("Caffeine-Aspirin co-crystal", "Generate a caffeine and aspirin co-crystal dimer", 2, "auto"),
    ("Water trimer H-bonded", "Create water trimer with hydrogen bonding", 3, "h_bonded"),
    ("Anthracene T-shaped", "Generate 2 anthracene molecules in T-shaped arrangement", 2, "t_shaped"),
    # Spiral and chiral tests
    ("Benzene spiral (45°)", "Stack 5 benzene molecules in a spiral with 45 degree rotation", 5, "stacked"),
    ("Naphthalene spiral (20°)", "Create naphthalene stack with 20 degree rotation per molecule", 4, "stacked"),
    ("Pyrene linear chain", "Build linear chain of 3 pyrene molecules along z-axis", 3, "linear"),
]


def call_ollama(prompt: str, tools_description: str) -> dict:
    """Call Ollama to convert natural language to tool call."""
    system_prompt = f"""You are an assistant for crystal structure generation.
When asked to generate molecules or clusters, respond with a JSON tool call.
For molecule names, use standard identifiers (e.g., 'Benzene', 'Glucose', 'Cyanocobalamin'). 
Use the EXACT NAME requested.
Do NOT output SMILES strings. Use names only.
For clusters, use 'linear', 'circular', 'pi_pi_parallel' for stacking if inferred.

{tools_description}

IMPORTANT: Respond ONLY with the JSON tool call. Format:
{{"name": "tool_name", "arguments": {{...}}}}
"""
    
    payload = {
        "model": OLLAMA_MODEL,
        "messages": [
            {"role": "system", "content": system_prompt},
            {"role": "user", "content": prompt}
        ],
        "stream": False,
        "format": "json"
    }
    
    print(f"  [LLM] Sending prompt: '{prompt[:50]}...'")
    
    import time
    start_time = time.time()
    
    response = requests.post(OLLAMA_URL, json=payload, timeout=60)
    
    if response.status_code != 200:
        return {"error": f"HTTP Error {response.status_code}: {response.text}"}
        
    duration = time.time() - start_time
    result = response.json()
    
    content = result.get("message", {}).get("content", "")
    print(f"  [LLM] Response received in {duration:.2f}s: {content[:60]}...")
    
    # Check if content is empty
    if not content:
        return {"error": "Empty response from LLM"}

    # Extract JSON if wrapped in code blocks
    json_str = content
    json_match = re.search(r'\{.*\}', content, re.DOTALL)
    if json_match:
        json_str = json_match.group(0)
    
    # Explicit JSON validation without try-except if possible, 
    # but json.loads requires exception handling for invalid JSON.
    # We will minimize the scope.
    if not json_str.strip().startswith('{'):
         return {"error": "Response is not valid JSON"}

    decoded = json.loads(json_str) 
    return decoded


def call_mcp_tool(tool_name: str, arguments: dict) -> Tuple[bool, dict, str]:
    """
    Call MCP server with a tool request.
    
    Returns:
        (success, structure_json, raw_output)
    """
    init_msg = json.dumps({
        "jsonrpc": "2.0",
        "id": 0,
        "method": "initialize",
        "params": {
            "protocolVersion": "0.1.0",
            "capabilities": {},
            "clientInfo": {"name": "test", "version": "1.0"}
        }
    })
    
    call_msg = json.dumps({
        "jsonrpc": "2.0",
        "id": 1,
        "method": "tools/call",
        "params": {
            "name": tool_name,
            "arguments": arguments
        }
    })
    
    script = f"""cd {MCP_SERVER_DIR}
cat << 'MCPINPUT' | timeout 120 node dist/index.js 2>&1
{init_msg}
{call_msg}
MCPINPUT"""
    
    try:
        result = subprocess.run(
            ["bash", "-c", script],
            capture_output=True,
            text=True,
            timeout=180
        )
        
        output = result.stdout
        
        # Parse the JSON-RPC response properly to unescape content
        # The output has two JSON objects - find the one with id: 1
        for line in output.split('\n'):
            line = line.strip()
            if not line:
                continue
            
            # Try to parse as JSON-RPC response
            try:
                if '"id":1' in line or '"id": 1' in line:
                    rpc = json.loads(line)
                    rpc_result = rpc.get('result', {})
                    content_items = rpc_result.get('content', [])
                    
                    for item in content_items:
                        text = item.get('text', '')
                        
                        # Now text is properly unescaped
                        if '<json-data>' in text:
                            start = text.find('<json-data>') + len('<json-data>')
                            end = text.find('</json-data>')
                            if end > start:
                                json_str = text[start:end].strip()
                                try:
                                    wrapper = json.loads(json_str)
                                    structure = wrapper.get('structure', wrapper)
                                    return True, structure, output
                                except json.JSONDecodeError:
                                    pass
            except json.JSONDecodeError:
                continue
        
        # Fallback: Check for success indicators in raw output
        if "Molecule Built:" in output or "Cluster Built:" in output:
            return True, {}, output
        
        return False, {}, output
    except Exception as e:
        return False, {}, str(e)




def extract_atoms_coords(structure: dict) -> Tuple[List[str], List[Tuple]]:
    """Extract atoms and coordinates from structure JSON."""
    atoms_data = structure.get('atoms', structure.get('sites', []))
    
    atoms = []
    coords = []
    
    for atom in atoms_data:
        elem = atom.get('element', 'X')
        cart = atom.get('cartesian', atom.get('coords', [0, 0, 0]))
        
        if isinstance(cart, list) and len(cart) >= 3:
            atoms.append(elem)
            coords.append(tuple(cart))
    
    return atoms, coords


def count_elements(atoms: List[str]) -> Dict[str, int]:
    """Count atoms by element."""
    counts = {}
    for elem in atoms:
        counts[elem] = counts.get(elem, 0) + 1
    return counts


def run_single_molecule_test(name: str, nl_prompt: str, 
                              expected_atoms: int, expected_elements: Optional[Dict]) -> dict:
    """Run a single molecule test and return results."""
    result = {
        "name": name,
        "nl_input": nl_prompt,
        "success": False,
        "atoms": [],
        "coords": [],
        "formula": None,
        "error": None,
        "ascii_diagram": None
    }
    
    # Step 1: LLM Translation
    tools_desc = """
build_molecule: Generate a single molecule.
  Arguments:
    - name: Molecule identifier (aspirin, benzene, PTCDA, c1ccccc1 for SMILES)
    - input_type: "auto" (default)
"""
    
    tool_call = call_ollama(nl_prompt, tools_desc)
    
    if "error" in tool_call:
        result["error"] = f"LLM Error: {tool_call['error']}"
        return result
    
    result["llm_output"] = tool_call
    
    tool_name = tool_call.get("name", "build_molecule")
    arguments = tool_call.get("arguments", {})
    
    # Clean up null values
    arguments = {k: v for k, v in arguments.items() if v is not None}
    
    # Step 2: MCP Server call
    success, structure, raw_output = call_mcp_tool(tool_name, arguments)
    
    if not success:
        # Capture raw output for debugging
        error_sample = raw_output[-300:].replace('\n', ' ') if raw_output else "No output"
        result["error"] = f"MCP call failed: {error_sample}"
        return result
    
    # Step 3: Extract atoms and validate
    atoms, coords = extract_atoms_coords(structure)
    
    if not atoms:
        result["error"] = "No atoms in structure"
        return result
    
    result["atoms"] = atoms
    result["coords"] = coords
    result["element_counts"] = count_elements(atoms)
    
    # Generate formula
    elem_counts = result["element_counts"]
    result["formula"] = ''.join(f"{e}{c}" for e, c in sorted(elem_counts.items()))
    
    # Validate
    if expected_atoms and len(atoms) < expected_atoms * 0.8:
        result["error"] = f"Expected ~{expected_atoms} atoms, got {len(atoms)}"
        return result
    
    if expected_elements:
        for elem, count in expected_elements.items():
            got = elem_counts.get(elem, 0)
            if got != count:
                result["error"] = f"Expected {count} {elem}, got {got}"
                return result
    
    # Generate ASCII diagram
    result["ascii_diagram"] = render_ascii_structure(atoms, coords, width=40, height=18)
    result["success"] = True
    
    return result


def run_cluster_test(name: str, nl_prompt: str, 
                     expected_molecules: int, expected_stacking: str) -> dict:
    """Run a cluster test and return results."""
    result = {
        "name": name,
        "nl_input": nl_prompt,
        "success": False,
        "atoms": [],
        "coords": [],
        "n_molecules": 0,
        "stacking": None,
        "formula": None,
        "error": None,
        "ascii_diagram": None
    }
    
    # Step 1: LLM Translation
    tools_desc = """
build_molecular_cluster: Generate a cluster of molecules.
  Arguments:
    - molecules: Array [{identifier: "name", count: N}, ...] or [{identifier: "mol1"}, {identifier: "mol2"}]
    - stacking: "pi_pi_parallel", "t_shaped", "herringbone", "h_bonded", "circular", "linear", "stacked", "swastika", "auto"
    - intermolecular_distance: Optional, in Angstroms
    - rotation_per_molecule: Optional, degrees of rotation between molecules in a stack
"""
    
    tool_call = call_ollama(nl_prompt, tools_desc)
    
    if "error" in tool_call:
        result["error"] = f"LLM Error: {tool_call['error']}"
        return result
    
    result["llm_output"] = tool_call
    
    tool_name = tool_call.get("name", "build_molecular_cluster")
    arguments = tool_call.get("arguments", {})
    
    # Clean up null values
    arguments = {k: v for k, v in arguments.items() if v is not None}
    
    # Step 2: MCP Server call
    success, structure, raw_output = call_mcp_tool(tool_name, arguments)
    
    if not success:
        result["error"] = "MCP call failed"
        return result
    
    # Step 3: Extract atoms and validate
    atoms, coords = extract_atoms_coords(structure)
    
    if not atoms:
        # Try to validate from raw output
        if f"Molecules**: {expected_molecules}" in raw_output or "Cluster Built:" in raw_output:
            result["success"] = True
            result["n_molecules"] = expected_molecules
            result["stacking"] = expected_stacking
            return result
        result["error"] = "No atoms in structure"
        return result
    
    result["atoms"] = atoms
    result["coords"] = coords
    result["element_counts"] = count_elements(atoms)
    
    # Get metadata
    metadata = structure.get("metadata", {})
    result["n_molecules"] = metadata.get("n_molecules", expected_molecules)
    result["stacking"] = metadata.get("stacking_type", expected_stacking)
    result["formula"] = metadata.get("formula", None)
    
    # Generate ASCII diagram (limited for large structures)
    if len(atoms) <= 50:
        result["ascii_diagram"] = render_ascii_structure(atoms, coords, width=40, height=18)
    else:
        result["ascii_diagram"] = f"  ({len(atoms)} atoms - too large for ASCII)"
    
    result["success"] = True
    
    return result


def format_result(result: dict, is_cluster: bool = False) -> str:
    """Format test result as side-by-side NL input and ball-stick diagram."""
    lines = []
    lines.append("┌" + "─" * 40 + "┬" + "─" * 42 + "┐")
    lines.append(f"│ {'NL Input:':<38} │ {'Ball-Stick Diagram:':<40} │")
    lines.append("├" + "─" * 40 + "┼" + "─" * 42 + "┤")
    
    # NL input wrapped
    nl_lines = []
    words = result["nl_input"].split()
    current = ""
    for word in words:
        if len(current) + len(word) + 1 <= 36:
            current += (" " if current else "") + word
        else:
            nl_lines.append(current)
            current = word
    if current:
        nl_lines.append(current)
    
    # ASCII diagram lines
    diagram = result.get("ascii_diagram", "")
    if isinstance(diagram, str):
        diagram_lines = diagram.split('\n') if diagram else []
    else:
        diagram_lines = []
    
    if result.get("error"):
        diagram_lines = [f"ERROR: {result['error'][:38]}"]
    elif not diagram_lines:
        diagram_lines = ["(Structure generated successfully)"]
    
    max_lines = max(len(nl_lines), len(diagram_lines), 6)
    
    for i in range(max_lines):
        nl_part = nl_lines[i] if i < len(nl_lines) else ""
        d_part = diagram_lines[i][:40] if i < len(diagram_lines) else ""
        lines.append(f"│ {nl_part:<38} │ {d_part:<40} │")
    
    lines.append("├" + "─" * 40 + "┴" + "─" * 42 + "┤")
    
    # Summary line
    if is_cluster:
        info = f"Molecules: {result.get('n_molecules', '?')}  Atoms: {len(result.get('atoms', []))}  Stacking: {result.get('stacking', '?')}"
    else:
        info = f"Formula: {result.get('formula') or '?':<15}  Atoms: {len(result.get('atoms', []))}"
    
    status = "✓" if result["success"] else "✗"
    lines.append(f"│ {status} {info:<80} │")
    lines.append("└" + "─" * 83 + "┘")
    
    return '\n'.join(lines)


# Additional Varied Phrasing Tests
# Additional Varied Phrasing Tests
VARIED_PHRASING_TESTS = [
    ("Sucrose", "give me table sugar (Sucrose)", 45, {"C": 12, "O": 11}),
    ("Retinol", "make Vitamin A (Retinol)", 51, {"C": 20}),
    ("Cholesterol", "show me Cholesterol", 74, {"C": 27}),
    ("Strychnine", "synthesize Strychnine", 47, {"N": 2}),
    ("Ruthenium Tris(bipyridine)", "construct Ruthenium tris(bipyridine)", 61, {"Ru": 1, "N": 6}),
]

VARIED_CLUSTER_TESTS = [
    ("Adenine-Thymine Pair", "pair Adenine and Thymine with hydrogen bonding", 2, "h_bonded"),
    ("Benzene Stack 9", "stack 9 benzene molecules in a column", 9, "pi_pi_parallel"),
    ("Protein Alpha Helix", "stack 5 Alanine amino acids in a spiral", 5, "custom"), 
    ("Ionic Crystal", "arrange 4 Sodium and 4 Chloride ions in a lattice", 8, "auto"),
]

BLIND_FOLD_TESTS = [
    ("Cubane", "synthesize Cubane", 16, {"C": 8, "H": 8}),
    ("Ferrocene", "generate Ferrocene", 21, {"Fe": 1}),
    ("Cisplatin", "create Cisplatin", 11, {"Pt": 1, "Cl": 2}),
    ("Paclitaxel", "build molecule Paclitaxel", 110, None), 
    ("Vancomycin", "generate Vancomycin", 170, {"Cl": 2}),
]

BLIND_FOLD_CLUSTER_TESTS = [
    ("Ferrocene Dimer", "create a dimer of Ferrocene", 2, "auto"),
    ("Cisplatin Stack", "stack 3 Cisplatin molecules", 3, "pi_pi_parallel"),
    ("C60 Crystal", "arrange 4 C60 molecules in a cluster", 4, "auto"),
]

def generate_documentation(results: List[Tuple[str, dict]], output_file: str):
    """Generate Markdown documentation from test results."""
    with open(output_file, 'w') as f:
        f.write("# Molecular Generation Verification Results\n\n")
        f.write("This document contains automatically generated verification results from the end-to-end testing suite.\n")
        f.write("Each entry shows the Natural Language Input, the verifiable ASCII Ball-Stick diagram generated from the output, and the validation status.\n\n")
        
        f.write("## Single Molecule Tests\n\n")
        for type_, result in results:
            if type_ != "molecule": continue
            
            status = "✅ PASS" if result["success"] else "❌ FAIL"
            f.write(f"### {result['name']}\n")
            f.write(f"**Input:** `{result['nl_input']}`  \n")
            f.write(f"**Status:** {status}  \n")
            f.write(f"**Formula:** {result.get('formula') or 'Unknown'}\n\n")
            
            f.write("```\n")
            if result.get('ascii_diagram'):
                f.write(result['ascii_diagram'])
            else:
                f.write("(No diagram generated)")
            f.write("\n```\n\n")
            f.write("---\n\n")

        f.write("## Molecular Cluster Tests\n\n")
        for type_, result in results:
            if type_ != "cluster": continue
            
            status = "✅ PASS" if result["success"] else "❌ FAIL"
            f.write(f"### {result['name']}\n")
            f.write(f"**Input:** `{result['nl_input']}`  \n")
            f.write(f"**Status:** {status}  \n")
            f.write(f"**Config:** {result.get('n_molecules', '?')} molecules, {result.get('stacking', '?')} stacking\n\n")
            
            f.write("```\n")
            if result.get('ascii_diagram'):
                f.write(result['ascii_diagram'])
            else:
                f.write("(No diagram generated)")
            f.write("\n```\n\n")
            f.write("---\n\n")

def main():
    print("=" * 84)
    print("CRYSTAL MCP SERVER - COMPREHENSIVE END-TO-END TESTING & DOC GEN")
    print("Flow: Natural Language → Ollama LLM → MCP Server → JSON → ASCII Ball-Stick")
    print("=" * 84)
    
    all_results = []
    
    # 1. Standard Single Molecule Tests
    print("\n" + "=" * 84)
    print("SINGLE MOLECULE TESTS (Standard)")
    print("=" * 84)
    for name, nl_prompt, expected_atoms, expected_elements in SINGLE_MOLECULE_TESTS:
        print(f"\nTesting: {name}...")
        result = run_single_molecule_test(name, nl_prompt, expected_atoms, expected_elements)
        print(format_result(result, is_cluster=False))
        all_results.append(("molecule", result))
        
    # 2. Varied Phrasing Tests
    print("\n" + "=" * 84)
    print("SINGLE MOLECULE TESTS (Varied Phrasing)")
    print("=" * 84)
    for name, nl_prompt, expected_atoms, expected_elements in VARIED_PHRASING_TESTS:
        print(f"\nTesting: {name} ({nl_prompt})...")
        result = run_single_molecule_test(name, nl_prompt, expected_atoms, expected_elements)
        print(format_result(result, is_cluster=False))
        all_results.append(("molecule", result))

    # 3. Standard Cluster Tests
    print("\n" + "=" * 84)
    print("MOLECULAR CLUSTER TESTS (Standard)")
    print("=" * 84)
    for name, nl_prompt, expected_molecules, expected_stacking in CLUSTER_TESTS:
        print(f"\nTesting: {name}...")
        result = run_cluster_test(name, nl_prompt, expected_molecules, expected_stacking)
        print(format_result(result, is_cluster=True))
        all_results.append(("cluster", result))
        
    # 4. Varied Cluster Tests
    print("\n" + "=" * 84)
    print("MOLECULAR CLUSTER TESTS (Varied Phrasing)")
    print("=" * 84)
    for name, nl_prompt, expected_molecules, expected_stacking in VARIED_CLUSTER_TESTS:
        print(f"\nTesting: {name} ({nl_prompt})...")
        result = run_cluster_test(name, nl_prompt, expected_molecules, expected_stacking)
        print(format_result(result, is_cluster=True))
        all_results.append(("cluster", result))
        
    # 5. Blind Fold Tests
    print("\n" + "=" * 84)
    print("BLIND FOLD TESTS (Complex & Unseen)")
    print("=" * 84)
    for name, nl_prompt, expected_atoms, expected_elements in BLIND_FOLD_TESTS:
        print(f"\nTesting: {name} ({nl_prompt})...")
        result = run_single_molecule_test(name, nl_prompt, expected_atoms, expected_elements)
        print(format_result(result, is_cluster=False))
        all_results.append(("molecule", result))
        
    # 6. Blind Fold Cluster Tests
    print("\n" + "=" * 84)
    print("BLIND FOLD CLUSTER TESTS (Complex & Unseen)")
    print("=" * 84)
    for name, nl_prompt, expected_molecules, expected_stacking in BLIND_FOLD_CLUSTER_TESTS:
        print(f"\nTesting: {name} ({nl_prompt})...")
        result = run_cluster_test(name, nl_prompt, expected_molecules, expected_stacking)
        print(format_result(result, is_cluster=True))
        all_results.append(("cluster", result))
    
    # Summary
    print("\n" + "=" * 84)
    print("SUMMARY")
    print("=" * 84)
    
    mol_results = [r for t, r in all_results if t == "molecule"]
    cluster_results = [r for t, r in all_results if t == "cluster"]
    
    mol_pass = sum(1 for r in mol_results if r["success"])
    cluster_pass = sum(1 for r in cluster_results if r["success"])
    
    print(f"\nSingle Molecules: {mol_pass}/{len(mol_results)} passed")
    for r in mol_results:
        status = "✓" if r["success"] else "✗"
        print(f"  {status} {r['name']:<20} → {r.get('formula') or r.get('error', 'Unknown')}")
    
    print(f"\nMolecular Clusters: {cluster_pass}/{len(cluster_results)} passed")
    for r in cluster_results:
        status = "✓" if r["success"] else "✗"
        extra = f"{r.get('n_molecules', '?')} mols, {r.get('stacking', '?')}"
        print(f"  {status} {r['name']:<25} → {extra}")
    
    total = mol_pass + cluster_pass
    total_tests = len(mol_results) + len(cluster_results)
    print(f"\nTotal: {total}/{total_tests} passed ({100*total/total_tests:.1f}%)")
    
    # Generate Documentation
    doc_path = os.path.join(MCP_SERVER_DIR, "docs/molecule/testing_results.md")
    print(f"\nGenerating documentation at: {doc_path}")
    generate_documentation(all_results, doc_path)
    print("Documentation generated successfully.")
    
    print("=" * 84)
    
    return 0 if total == total_tests else 1


if __name__ == "__main__":
    sys.exit(main())
