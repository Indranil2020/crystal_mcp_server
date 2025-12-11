"""
quantum/ - Quantum Materials Generation

Category 15-16: Quantum dots, qubits, color centers.
"""

from .qubits import generate_qubit_host, QUBIT_DATABASE
from .dots import generate_quantum_dot_3d, QD_DATABASE


__all__ = ["generate_qubit_host", "QUBIT_DATABASE", "generate_quantum_dot_3d", "QD_DATABASE"]
