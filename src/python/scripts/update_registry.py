"""
update_registry.py - Automating Registry Generation

This script walks the src/python/generators directory, imports modules,
inspects functions, and generates the definitive generators/__init__.py file.
"""

import os
import pkgutil
import inspect
import importlib
import sys
from typing import Dict, Any, List

# Add src/python to path
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

SKIP_FILES = ["__init__.py"]
SKIP_DIRS = ["__pycache__"]

def get_function_info(func) -> Dict[str, Any]:
    """Extract params and docstring from function."""
    sig = inspect.signature(func)
    params = list(sig.parameters.keys())
    
    doc = inspect.getdoc(func)
    description = doc.split('\n')[0] if doc else f"Generate {func.__name__}"
    
    return {
        "params": params,
        "description": description
    }

def scan_generators():
    registry = {}
    base_path = os.path.join(os.path.dirname(__file__), "../generators")
    
    # Walk directory
    for root, dirs, files in os.walk(base_path):
        # Exclude pycache
        dirs[:] = [d for d in dirs if d not in SKIP_DIRS]
        
        # Get category from folder name
        rel_path = os.path.relpath(root, base_path)
        if rel_path == ".":
            continue
            
        category = rel_path.split(os.sep)[0]
        
        if category not in registry:
            registry[category] = {
                "description": f"{category.replace('_', ' ').title()} structures",
                "operations": {}
            }
            
        for file in files:
            if file in SKIP_FILES or not file.endswith(".py"):
                continue
                
            module_name = f"generators.{rel_path.replace(os.sep, '.')}.{file[:-3]}"
            
            try:
                module = importlib.import_module(module_name)
                
                # specific manual fixes for known problematic modules if import fails
                # (handled by try/except)
                
                for name, obj in inspect.getmembers(module):
                    if (inspect.isfunction(obj) and 
                        (name.startswith("generate_") or 
                         name.startswith("make_") or 
                         name.startswith("create_") or
                         name.startswith("build_") or
                         name.startswith("apply_") or
                         name.startswith("export_"))):
                        
                        # Verify it's defined in this module (not imported)
                        if obj.__module__ != module_name:
                            continue
                            
                        info = get_function_info(obj)
                        
                        registry[category]["operations"][name] = {
                            "module": module_name,
                            "function": name,
                            "params": info["params"],
                            "description": info["description"]
                        }
                        
            except ImportError as e:
                print(f"Warning: Could not import {module_name}: {e}")
            except Exception as e:
                print(f"Warning: Error inspecting {module_name}: {e}")

    return registry

def generate_init_file(registry):
    lines = [
        '"""',
        'generators/__init__.py - Central Registry for All Structure Generators',
        '',
        'This module provides a unified registry of all generator functions,',
        'enabling dynamic discovery and routing from the MCP layer.',
        '"""',
        '',
        'from typing import Dict, Any, List, Callable',
        '',
        '# =============================================================================',
        '# GENERATOR REGISTRY',
        '# =============================================================================',
        '# Each category maps to its available operations with metadata',
        '',
        'GENERATOR_REGISTRY: Dict[str, Dict[str, Any]] = {'
    ]
    
    # Sort categories
    for category in sorted(registry.keys()):
        data = registry[category]
        lines.append(f'    # -------------------------------------------------------------------------')
        lines.append(f'    # {category.upper()} (Category)')
        lines.append(f'    # -------------------------------------------------------------------------')
        lines.append(f'    "{category}": {{')
        lines.append(f'        "description": "{data["description"]}",')
        lines.append(f'        "operations": {{')
        
        # Sort operations
        for op_name in sorted(data["operations"].keys()):
            op_data = data["operations"][op_name]
            lines.append(f'            "{op_name}": {{')
            lines.append(f'                "module": "{op_data["module"]}",')
            lines.append(f'                "function": "{op_data["function"]}",')
            lines.append(f'                "params": {op_data["params"]},')
            lines.append(f'                "description": "{op_data["description"].replace(chr(34), chr(39))}"') # escape quotes
            lines.append(f'            }},')
            
        lines.append(f'        }}')
        lines.append(f'    }},')
        lines.append(f'    ')

    lines.append('}')
    lines.append('')
    lines.append('')
    lines.append('def get_all_operations() -> Dict[str, List[str]]:')
    lines.append('    """Get all available operations organized by category."""')
    lines.append('    result = {}')
    lines.append('    for category, data in GENERATOR_REGISTRY.items():')
    lines.append('        result[category] = list(data["operations"].keys())')
    lines.append('    return result')
    lines.append('')
    lines.append('')
    lines.append('def get_operation_info(operation: str) -> Dict[str, Any]:')
    lines.append('    """Get information about a specific operation."""')
    lines.append('    for category, data in GENERATOR_REGISTRY.items():')
    lines.append('        if operation in data["operations"]:')
    lines.append('            op_info = data["operations"][operation].copy()')
    lines.append('            op_info["category"] = category')
    lines.append('            return {"success": True, "operation": operation, "info": op_info}')
    lines.append('    ')
    lines.append('    return {')
    lines.append('        "success": False,')
    lines.append('        "error": {')
    lines.append('            "code": "UNKNOWN_OPERATION",')
    lines.append('            "message": f"Operation \'{operation}\' not found",')
    lines.append('            "available_operations": get_all_operations()')
    lines.append('        }')
    lines.append('    }')
    lines.append('')
    lines.append('')
    lines.append('def list_category_operations(category: str) -> Dict[str, Any]:')
    lines.append('    """List all operations in a category."""')
    lines.append('    if category not in GENERATOR_REGISTRY:')
    lines.append('        return {')
    lines.append('            "success": False,')
    lines.append('            "error": {')
    lines.append('                "code": "UNKNOWN_CATEGORY",')
    lines.append('                "message": f"Category \'{category}\' not found",')
    lines.append('                "available_categories": list(GENERATOR_REGISTRY.keys())')
    lines.append('            }')
    lines.append('        }')
    lines.append('    ')
    lines.append('    data = GENERATOR_REGISTRY[category]')
    lines.append('    return {')
    lines.append('        "success": True,')
    lines.append('        "category": category,')
    lines.append('        "description": data["description"],')
    lines.append('        "operations": {')
    lines.append('            name: info["description"] ')
    lines.append('            for name, info in data["operations"].items()')
    lines.append('        }')
    lines.append('    }')

    return "\n".join(lines)

if __name__ == "__main__":
    print("Scanning generators...")
    registry = scan_generators()
    
    total_ops = sum(len(d["operations"]) for d in registry.values())
    print(f"Found {len(registry)} categories and {total_ops} operations.")
    
    content = generate_init_file(registry)
    
    output_path = os.path.join(os.path.dirname(__file__), "../generators/__init__.py")
    with open(output_path, "w") as f:
        f.write(content)
        
    print(f"Written to {output_path}")
