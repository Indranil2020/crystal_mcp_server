#!/usr/bin/env python3
"""
Comprehensive Test Report Generator

Purpose: Generate extensive analysis report with ZERO truncation
Shows: Full Prompt → LLM Interpretation → MCP Output
Goal: Identify where to improve (LLM vs Backend)
"""

import json
import sys
import re
from typing import Dict, List, Any
from collections import defaultdict

def load_test_results(filename: str) -> List[Dict]:
    """Load test results JSON"""
    try:
        with open(filename, 'r') as f:
            return json.load(f)
    except Exception as e:
        print(f"ERROR loading {filename}: {e}")
        sys.exit(1)

def format_json(obj: Any, indent=2) -> str:
    """Format JSON with proper indentation"""
    return json.dumps(obj, indent=indent)

def analyze_llm_interpretation(prompt: str, llm_args: Dict) -> Dict[str, Any]:
    """Analyze how LLM interpreted the prompt"""
    analysis = {
        'understood_correctly': [],
        'misunderstood': [],
        'missing': [],
        'unexpected': []
    }
    
    prompt_lower = prompt.lower()
    
    # Check distance interpretation
    if 'angstrom' in prompt_lower or 'gap' in prompt_lower:
        for num in ['0.5', '1', '2', '2.8', '3.4', '4', '5', '10', '15']:
            if num in prompt:
                if llm_args and llm_args.get('intermolecular_distance') ==float(num):
                    analysis['understood_correctly'].append(f'distance={num}Å')
                elif llm_args:
                    analysis['misunderstood'].append(f'Expected distance={num}Å, got {llm_args.get("intermolecular_distance")}')
    
    # Check axis/direction
    if 'along x' in prompt_lower:
        if llm_args and llm_args.get('axis') == 'x':
            analysis['understood_correctly'].append('axis=x')
        elif llm_args:
            analysis['misunderstood'].append(f'Expected axis=x, got {llm_args.get("axis")}')
    
    if 'along y' in prompt_lower:
        if llm_args and llm_args.get('axis') == 'y':
            analysis['understood_correctly'].append('axis=y')
        elif llm_args:
            analysis['misunderstood'].append(f'Expected axis=y, got {llm_args.get("axis")}')
    
    if 'along z' in prompt_lower:
        if llm_args and llm_args.get('axis') == 'z':
            analysis['understood_correctly'].append('axis=z')
        elif llm_args:
            analysis['misunderstood'].append(f'Expected axis=z, got {llm_args.get("axis")}')
    
    # Check angle interpretation - KEY FAILURE POINT
    if 'degree direction' in prompt_lower:
        # Extract angle from prompt
        angle_match = re.search(r'(\d+)\s+degree', prompt_lower)
        if angle_match:
            angle = int(angle_match.group(1))
            if llm_args:
                if llm_args.get('direction_angle') == angle:
                    analysis['understood_correctly'].append(f'direction_angle={angle}°')
                elif 'rotation' in str(llm_args):
                    analysis['misunderstood'].append(f'CRITICAL: LLM interpreted {angle}° DIRECTION as ROTATION instead of placement axis')
                else:
                    analysis['missing'].append(f'direction_angle={angle}° not captured')
    
    # Check vector notation
    vector_match = re.search(r'\[(\d+),\s*(\d+),\s*(\d+)\]', prompt)
    if vector_match:
        vec = [int(vector_match.group(1)), int(vector_match.group(2)), int(vector_match.group(3))]
        if llm_args and llm_args.get('direction_vector') == vec:
            analysis['understood_correctly'].append(f'direction_vector={vec}')
        elif llm_args:
            analysis['missing'].append(f'direction_vector={vec} not captured')
    
    # Check molecule count
    mcount_match = re.search(r'(\d+)\s+(benzene|ptcda|water|c60)', prompt_lower)
    if mcount_match:
        count = int(mcount_match.group(1))
        if llm_args and llm_args.get('molecules'):
            molecules = llm_args['molecules']
            # Handle both list of dicts and list of strings
            if isinstance(molecules, list):
                if len(molecules) > 0 and isinstance(molecules[0], dict):
                    actual_count = sum(m.get('count', 1) for m in molecules)
                else:
                    actual_count = len(molecules)  # List of strings
                
                if actual_count == count:
                    analysis['understood_correctly'].append(f'molecule_count={count}')
                else:
                    analysis['misunderstood'].append(f'Expected count={count}, got {actual_count}')
    
    return analysis

def generate_markdown_report(results: List[Dict]) -> str:
    """Generate comprehensive Markdown report with ZERO truncation"""
    
    report = [
        "# COMPREHENSIVE TEST REPORT - MOLECULAR ARRANGEMENT SYSTEM",
        "",
        "**Generated:** " + results[0]['timestamp'] if results else "N/A",
        f"**Total Tests:** {len(results)}",
        "",
        "---",
        ""
    ]
    
    # Summary Statistics
    successes = [r for r in results if r['status'] == 'SUCCESS']
    mcp_failed = [r for r in results if r['status'] == 'MCP_FAILED']
    no_request = [r for r in results if r['status'] == 'NO_REQUEST']
    test_errors = [r for r in results if r['status'] == 'TEST_ERROR']
    
    report.extend([
        "## Executive Summary",
        "",
        f"| Status | Count | Rate |",
        f"|:---|:---:|:---:|",
        f"| ✅ SUCCESS | {len(successes)} | {len(successes)/len(results)*100:.1f}% |",
        f"| ❌ MCP_FAILED | {len(mcp_failed)} | {len(mcp_failed)/len(results)*100:.1f}% |",
        f"| ⚠️ NO_REQUEST | {len(no_request)} | {len(no_request)/len(results)*100:.1f}% |",
        f"| 🔥 TEST_ERROR | {len(test_errors)} | {len(test_errors)/len(results)*100:.1f}% |",
        "",
        "---",
        ""
    ])
    
    # Category-wise breakdown
    categories = defaultdict(lambda: {'success': 0, 'failed': 0, 'no_request': 0})
    for r in results:
        cat = r.get('category', 'Unknown')
        if r['status'] == 'SUCCESS':
            categories[cat]['success'] += 1
        elif r['status'] == 'MCP_FAILED':
            categories[cat]['failed'] += 1
        elif r['status'] == 'NO_REQUEST':
            categories[cat]['no_request'] += 1
    
    report.extend([
        "## Category Breakdown",
        "",
        "| Category | Success | MCP Failed | No Request |",
        "|:---|:---:|:---:|:---:|"
    ])
    
    for cat in sorted(categories.keys()):
        stats = categories[cat]
        report.append(f"| {cat} | {stats['success']} | {stats['failed']} | {stats['no_request']} |")
    
    report.extend(["", "---", ""])
    
    # Detailed Analysis - NO TRUNCATION
    report.extend([
        "## DETAILED TEST ANALYSIS",
        "",
        "### Format:",
        "- **Prompt:** Full user input (no truncation)",
        "- **LLM Interpretation:** Complete tool call with all arguments",
        "- **LLM Analysis:** What LLM understood vs misunderstood",
        "- **MCP Output:** Full response including metadata",
        "- **Verdict:** SUCCESS/FAILURE with specific failure reason",
        "",
        "---",
        ""
    ])
    
    for test in results:
        report.extend([
            f"## Test {test['id']}: {test['category']}",
            "",
            f"**Notes:** {test.get('notes', 'N/A')}",
            "",
            "### 📝 PROMPT (User Input)",
            "```",
            test['prompt'],  # FULL PROMPT - NO TRUNCATION
            "```",
            ""
        ])
        
        # LLM Interpretation
        report.extend([
            "### 🤖 LLM INTERPRETATION",
            "",
            f"**Tool Selected:** `{test.get('llm_tool', 'NONE')}`",
            ""
        ])
        
        if test.get('llm_args'):
            report.extend([
                "**Full Arguments:**",
                "```json",
                format_json(test['llm_args']),  # FULL ARGS - NO TRUNCATION
                "```",
                ""
            ])
            
            # Analyze LLM interpretation
            analysis = analyze_llm_interpretation(test['prompt'], test['llm_args'])
            report.extend([
                "**LLM Understanding Analysis:**",
                ""
            ])
            
            if analysis['understood_correctly']:
                report.append("✅ **Correctly Understood:**")
                for item in analysis['understood_correctly']:
                    report.append(f"  - {item}")
                report.append("")
            
            if analysis['misunderstood']:
                report.append("❌ **Misunderstood:**")
                for item in analysis['misunderstood']:
                    report.append(f"  - {item}")
                report.append("")
            
            if analysis['missing']:
                report.append("⚠️ **Missing:**")
                for item in analysis['missing']:
                    report.append(f"  - {item}")
                report.append("")
        else:
            report.extend([
                "**⚠️ NO LLM REQUEST GENERATED**",
                "",
                "This indicates the LLM failed to produce ANY tool call.",
                ""
            ])
        
        # MCP Response
        report.extend([
            "### ⚙️ MCP SERVER RESPONSE",
            "",
            f"**Success:** {test.get('mcp_success', False)}",
            ""
        ])
        
        if test.get('mcp_metadata'):
            report.extend([
                "**Metadata:**",
                "```json",
                format_json(test['mcp_metadata']),  # FULL METADATA
                "```",
                ""
            ])
        
        if test.get('mcp_error'):
            report.extend([
                "**Error:**",
                "```",
                str(test['mcp_error']),  # FULL ERROR
                "```",
                ""
            ])
        
        # Verdict
        report.extend([
            "### 🎯 VERDICT",
            ""
        ])
        
        if test['status'] == 'SUCCESS':
            report.append("✅ **SUCCESS** - Structure generated correctly")
        elif test['status'] == 'MCP_FAILED':
            report.extend([
                "❌ **MCP_FAILED**",
                "",
                "**Root Cause Analysis:**",
                test['llm_args'] is not None and "- LLM produced a tool call, but MCP backend rejected/failed it" or "- Unknown",
                "- Check if backend validation is too strict",
                "- Check if required parameters were missing",
                "- Check backend logs for specific error"
            ])
        elif test['status'] == 'NO_REQUEST':
            report.extend([
                "⚠️ **NO_REQUEST**",
                "",
                "**Root Cause Analysis:**",
                "- LLM failed to generate ANY tool call from this prompt",
                "- Possible causes:",
                "  1. Prompt too ambiguous/complex for current LLM",
                "  2. LLM doesn't understand this phrasing",
                "  3. Temperature too high causing inconsistent behavior",
                "",
                "**Recommendation:** Improve LLM model or refine prompt understanding"
            ])
        
        report.extend(["", "---", ""])
    
    # Root Cause Summary
    report.extend([
        "## 🔬 ROOT CAUSE ANALYSIS SUMMARY",
        "",
        "### LLM Failures",
        ""
    ])
    
    llm_issues = defaultdict(int)
    for test in results:
        if test['status'] == 'NO_REQUEST':
            llm_issues['No tool call generated'] += 1
        elif test.get('llm_args'):
            analysis = analyze_llm_interpretation(test['prompt'], test['llm_args'])
            for issue in analysis['misunderstood']:
                if 'ROTATION instead of placement' in issue:
                    llm_issues['Angle→Rotation confusion'] += 1
                else:
                    llm_issues[issue] += 1
    
    if llm_issues:
        for issue, count in sorted(llm_issues.items(), key=lambda x: -x[1]):
            report.append(f"- **{issue}:** {count} occurrences")
    else:
        report.append("- No LLM issues detected")
    
    report.extend(["", "### Backend Failures", ""])
    
    backend_issues = defaultdict(int)
    for test in mcp_failed:
        if test.get('llm_args'):
            # Try to identify backend issue
            args = test['llm_args']
            if 'molecules' not in args:
                backend_issues['Missing molecules parameter'] += 1
            elif not args.get('molecules'):
                backend_issues['Empty molecules array'] += 1
            else:
                backend_issues['Unknown backend rejection'] += 1
    
    if backend_issues:
        for issue, count in sorted(backend_issues.items(), key=lambda x: -x[1]):
            report.append(f"- **{issue}:** {count} occurrences")
    else:
        report.append("- No backend issues detected")
    
    report.extend([
        "",
        "---",
        "",
        "## 💡 RECOMMENDATIONS",
        "",
        "### LLM Improvements Needed:",
        ""
    ])
    
    if llm_issues:
        if llm_issues.get('Angle→Rotation confusion', 0) > 0:
            report.extend([
                "1. **CRITICAL: Fix Angle Interpretation**",
                "   - LLM confuses 'N degree direction' with rotation",
                "   - Need to add examples to system prompt:",
                "     - '30 degree direction' → `direction_angle: 30`",
                "     - NOT `rotation_per_molecule: 30`",
                ""
            ])
        
        if llm_issues.get('No tool call generated', 0) > 0:
            report.extend([
                "2. **Improve Prompt Understanding**",
                "   - Some prompts produce NO tool calls",
                "   - Consider:",
                "     - Lower temperature (set to 0 for determinism)",
                "     - Better examples in system prompt",
                "     - More capable LLM model",
                ""
            ])
    
    report.extend([
        "### Backend Improvements Needed:",
        ""
    ])
    
    if backend_issues:
        report.extend([
            "1. **Relax Parameter Validation**",
            "   - Backend may be rejecting valid LLM outputs",
            "   - Add better error messages",
            "   - Make optional parameters truly optional",
            ""
        ])
    
    report.extend([
        "---",
        "",
        "## 📊 TEST COVERAGE MATRIX",
        "",
        "| Feature | Tested | Success | Rate |",
        "|:---|:---:|:---:|:---:|"
    ])
    
    # Feature coverage
    features = {
        'Gap-based positioning': [r for r in results if 'Gap' in r.get('category', '')],
        'Direction angles': [r for r in results if 'Angle' in r.get('category', '')],
        'Direction vectors': [r for r in results if 'Vec' in r.get('category', '')],
        'Auto-relaxation': [r for r in results if 'Relax' in r.get('category', '')],
        'Patterns': [r for r in results if 'Pattern' in r.get('category', '')],
        'Constraints': [r for r in results if 'Constr' in r.get('category', '')],
        'Formulas': [r for r in results if 'Form' in r.get('category', '')],
        'Heterogeneous': [r for r in results if 'Hetero' in r.get('category', '')],
        'Rotation': [r for r in results if 'Rot' in r.get('category', '')],
        'Edge cases': [r for r in results if 'Edge' in r.get('category', '')],
    }
    
    for feature, tests in features.items():
        if tests:
            success = len([t for t in tests if t['status'] == 'SUCCESS'])
            rate = success / len(tests) * 100
            report.append(f"| {feature} | {len(tests)} | {success} | {rate:.1f}% |")
    
    report.extend(["", "---", ""])
    
    return "\n".join(report)

def main():
    """Main entry point"""
    results_file = 'tests/molecule_test/comprehensive_test_report.json'
    output_file = 'tests/molecule_test/COMPREHENSIVE_ANALYSIS_REPORT.md'
    
    print(f"Loading test results from: {results_file}")
    results = load_test_results(results_file)
    print(f"Loaded {len(results)} test results")
    
    print("Generating comprehensive report...")
    report = generate_markdown_report(results)
    
    print(f"Writing report to: {output_file}")
    with open(output_file, 'w') as f:
        f.write(report)
    
    print(f"✅ Report generated: {output_file}")
    print(f"   Total lines: {len(report.split(chr(10)))}")
    print(f"   Total bytes: {len(report)}")

if __name__ == '__main__':
    main()
