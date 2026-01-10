#!/usr/bin/env python3
"""
molecule_suggester.py - Molecule Name Suggestion Tool

Provides intelligent molecule name suggestions when exact match fails.
Designed to work with both:
- Standalone MCP server (Claude, GPT, etc.)
- Web GUI with local LLM (Qwen Instruct)

Uses server-side intelligence so even simple LLMs can benefit.
"""

import sys
import json
import logging

# Setup logging to stderr
logging.basicConfig(
    level=logging.DEBUG,
    format='[%(asctime)s] %(levelname)s: %(message)s',
    stream=sys.stderr
)
logger = logging.getLogger(__name__)


def main():
    """Main entry point - reads JSON from stdin, outputs JSON to stdout."""
    # Read input from stdin
    input_data = json.load(sys.stdin)

    query = input_data.get("query", "")
    max_results = input_data.get("max_results", 10)
    include_smiles = input_data.get("include_smiles", True)
    min_similarity = input_data.get("min_similarity", 0.3)

    logger.debug(f"Suggesting molecules for: {query}")

    # Import the suggestion function
    from generators.molecule.universal_molecule import suggest_molecules

    # Get suggestions
    result = suggest_molecules(
        query=query,
        max_results=max_results,
        include_smiles=include_smiles,
    )

    # Filter by minimum similarity
    if result.get("suggestions"):
        result["suggestions"] = [
            s for s in result["suggestions"]
            if s.get("score", 0) >= min_similarity
        ]

    # Output result as JSON
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
