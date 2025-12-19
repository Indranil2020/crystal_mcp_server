"""
quantum/ - Quantum Structures (Dots, Wires, Wells)

Covers Category 8 items from structure_catalogue.md.
"""

from .dots import (
    generate_quantum_dot,
    generate_coreshell_qd,
    calculate_bandgap_qd,
    QUANTUM_DOT_DATABASE,
    CORESHELL_DATABASE,
    LIGAND_DATABASE,
)


__all__ = [
    "generate_quantum_dot", "generate_coreshell_qd", "calculate_bandgap_qd",
    "QUANTUM_DOT_DATABASE", "CORESHELL_DATABASE", "LIGAND_DATABASE",
]
