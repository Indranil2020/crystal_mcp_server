
"""
generate_docs.py

Auto-generates API Reference documentation from the live registry.
"""
import os
import sys
import inspect

# Add path to find generators
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

from generators import GENERATOR_REGISTRY

def generate_markdown_docs():
    lines = []
    lines.append("# Crystal MCP Generator API Reference")
    lines.append("")
    lines.append("Complete catalog of all **228** available generator operations, organized by category.")
    lines.append("")
    
    total_ops = 0
    
    # Sort categories
    categories = sorted(GENERATOR_REGISTRY.keys())
    
    # Table of Contents
    lines.append("## Categories")
    for cat in categories:
        count = len(GENERATOR_REGISTRY[cat]["operations"])
        total_ops += count
        desc = GENERATOR_REGISTRY[cat]["description"]
        lines.append(f"- [**{cat.upper()}**](#{cat}) ({count} ops): {desc}")
    
    lines.append("")
    lines.append("---")
    lines.append("")
    
    # Details
    for cat in categories:
        data = GENERATOR_REGISTRY[cat]
        lines.append(f"## {cat}")
        lines.append(f"_{data['description']}_")
        lines.append("")
        
        # Table header
        lines.append("| Operation | Description | Parameters |")
        lines.append("|-----------|-------------|------------|")
        
        for op in sorted(data["operations"].keys()):
            op_data = data["operations"][op]
            desc = op_data["description"]
            params = ", ".join([f"`{p}`" for p in op_data["params"]])
            lines.append(f"| `{op}` | {desc} | {params} |")
            
        lines.append("")
        lines.append(f"[Back to Top](#categories)")
        lines.append("")
        lines.append("---")
        lines.append("")

    return "\n".join(lines)

if __name__ == "__main__":
    content = generate_markdown_docs()
    
    # Write to docs/api_reference.md
    # assuming we run from project root or similar
    # path relative to this script: ../../../docs/api_reference.md
    output_path = os.path.abspath(os.path.join(os.path.dirname(__file__), "../../../docs/api_reference.md"))
    
    with open(output_path, "w") as f:
        f.write(content)
        
    print(f"Generated API docs at: {output_path}")
