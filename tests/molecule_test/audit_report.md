# Prompt-Based Audit Report: "Ultrafine" Level Scenarios
**Date:** 10 Jan 2026
**Auditor:** Crystal Assistant (Agentic Mode)

## Overview
We executed a "One-by-One" Deep Audit of the stress-test prompts defined in `test_prompt.md`.
**Objective**: Determine "what is really producing" and identify specific failure points.
**Scope**: Level 1 (Topology), Level 2 (Supramolecular), Level 5 (Scale).

## Test Case 1: Polycatenane (Level 1)
**Prompt**: *"Generate a polycatenane consisting of three interlocked macrocycles... 180° apart..."*
- **Status**: ❌ **FAILED**
- **Capabilities Tested**: Complex Topology, Unknown Molecule Name.
- **Log Verification**:
  ```text
  [DEBUG arrangement_adapter] generate_molecular_cluster called
  [DEBUG arrangement_adapter]   stacking: custom
  [DEBUG arrangement_adapter]   use_new_engine: False  <-- FELL BACK TO LEGACY
  [DEBUG arrangement_adapter] Legacy result: success=False
  ```
- **Root Cause**:
  1.  **Frontend Gap**: The Assistant did NOT generate `constraints` or `formulas` for the topology. It defaulted to legacy `custom` stacking.
  2.  **Data Gap**: Backend failed to resolve molecule name `"1,4-phenylene-2,3-dicarboximide"`.

## Test Case 2: HOF Diamondoid Network (Level 2)
**Prompt**: *"Create a 3D hydrogen-bonded organic framework (HOF)... 18 guanidinium, 18 sulfonate..."*
- **Status**: ❌ **FAILED** (Incorrect Structure)
- **Visual Result**: Single dense cluster of 18 molecules (Guanidinium only).
- **Log Verification**:
  - `generate_molecular_cluster` received only ONE molecule type list.
  - `sulfonate` was seemingly dropped by Frontend parsing.
  - `ConstraintSolver` never ran.
- **Root Cause**:
  - **Frontend Gap**: Failed to parse multi-component list correctly.
  - **Frontend Gap**: Did not translate diamondoid network into `constraints`.

## Test Case 3: Lipid Raft (Level 5)
**Prompt**: *"Create a full lipid raft model: 50 POPC, 30 sphingomyelin, 20 cholesterol..."*
- **Status**: ❌ **TIMEOUT** (>60s)
- **Capabilities Tested**: Scale (100 biomolecules), New Engine.
- **Log Verification**:
  ```text
  [DEBUG arrangement_adapter] use_new_engine: True  <-- SUCCESS: NEW ENGINE USED
  [DEBUG arrangement_adapter] Resolving 'POPC' -> Success (RDKit)
  [DEBUG arrangement_adapter] Resolving 'sphingomyelin' -> FAILED
  [DEBUG arrangement_adapter] Calling engine.arrange_molecules... (TIMEOUT)
  ```
- **Root Cause**:
  1.  **Computation limit**: Python Solver cannot optimize 70+ large molecules in <60s.
  2.  **Data Gap**: Missing `sphingomyelin`.

## Conclusion
The System is **Hardened** for geometric tasks (verified in `test_arrangement_geometry.py`) but **NOT Ready** for the "Ultrafine" natural language prompts.
- **The Engine works** when inputs are correct (e.g. `use_new_engine: True` in Prompt 22).
- **The End-to-End Pipeline fails** because the Frontend LLM is not expertly translating complex English into the Engine's `constraints` API.
