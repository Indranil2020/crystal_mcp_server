"""
molecule/biomolecules.py - Biological Macromolecules

Generates biomolecular structures:
- DNA/RNA sequences
- Peptides (alpha-helix, beta-sheet)
- G-quadruplexes
- Peptide libraries for surface binding
"""

from typing import Dict, Any, List, Optional, Tuple
import numpy as np


# Amino acid database
AMINO_ACIDS = {
    "A": {"name": "Alanine", "3letter": "Ala", "atoms": 5, "charge": 0, "polar": False},
    "R": {"name": "Arginine", "3letter": "Arg", "atoms": 11, "charge": 1, "polar": True},
    "N": {"name": "Asparagine", "3letter": "Asn", "atoms": 8, "charge": 0, "polar": True},
    "D": {"name": "Aspartate", "3letter": "Asp", "atoms": 8, "charge": -1, "polar": True},
    "C": {"name": "Cysteine", "3letter": "Cys", "atoms": 6, "charge": 0, "polar": True},
    "E": {"name": "Glutamate", "3letter": "Glu", "atoms": 9, "charge": -1, "polar": True},
    "Q": {"name": "Glutamine", "3letter": "Gln", "atoms": 9, "charge": 0, "polar": True},
    "G": {"name": "Glycine", "3letter": "Gly", "atoms": 4, "charge": 0, "polar": False},
    "H": {"name": "Histidine", "3letter": "His", "atoms": 10, "charge": 0, "polar": True},
    "I": {"name": "Isoleucine", "3letter": "Ile", "atoms": 8, "charge": 0, "polar": False},
    "L": {"name": "Leucine", "3letter": "Leu", "atoms": 8, "charge": 0, "polar": False},
    "K": {"name": "Lysine", "3letter": "Lys", "atoms": 9, "charge": 1, "polar": True},
    "M": {"name": "Methionine", "3letter": "Met", "atoms": 8, "charge": 0, "polar": False},
    "F": {"name": "Phenylalanine", "3letter": "Phe", "atoms": 11, "charge": 0, "polar": False},
    "P": {"name": "Proline", "3letter": "Pro", "atoms": 7, "charge": 0, "polar": False},
    "S": {"name": "Serine", "3letter": "Ser", "atoms": 6, "charge": 0, "polar": True},
    "T": {"name": "Threonine", "3letter": "Thr", "atoms": 7, "charge": 0, "polar": True},
    "W": {"name": "Tryptophan", "3letter": "Trp", "atoms": 14, "charge": 0, "polar": False},
    "Y": {"name": "Tyrosine", "3letter": "Tyr", "atoms": 12, "charge": 0, "polar": True},
    "V": {"name": "Valine", "3letter": "Val", "atoms": 7, "charge": 0, "polar": False},
}


# Nucleotide database
NUCLEOTIDES = {
    "A": {"name": "Adenine", "type": "purine", "pairs_with": "T", "atoms": 15},
    "T": {"name": "Thymine", "type": "pyrimidine", "pairs_with": "A", "atoms": 15},
    "G": {"name": "Guanine", "type": "purine", "pairs_with": "C", "atoms": 16},
    "C": {"name": "Cytosine", "type": "pyrimidine", "pairs_with": "G", "atoms": 14},
    "U": {"name": "Uracil", "type": "pyrimidine", "pairs_with": "A", "atoms": 12},
}


# Secondary structure parameters
HELIX_PARAMS = {
    "alpha": {"rise": 1.5, "rotation": 100, "radius": 2.3, "residues_per_turn": 3.6},
    "310": {"rise": 2.0, "rotation": 120, "radius": 1.9, "residues_per_turn": 3.0},
    "pi": {"rise": 1.15, "rotation": 87, "radius": 2.8, "residues_per_turn": 4.4},
    "PPII": {"rise": 3.12, "rotation": 120, "radius": 4.5, "residues_per_turn": 3.0},
}


def generate_peptide(
    sequence: str = "ALANINE",
    secondary_structure: str = "alpha",
    termini: str = "NH2-COOH"
) -> Dict[str, Any]:
    """
    Generate peptide structure.
    
    Args:
        sequence: Amino acid sequence (1-letter codes)
        secondary_structure: 'alpha', '310', 'pi', 'PPII', 'beta', 'coil'
        termini: Terminal groups 'NH2-COOH', 'acetyl-amide', 'neutral'
    
    Returns:
        Peptide structure
    """
    # Validate sequence
    for aa in sequence:
        if aa not in AMINO_ACIDS:
            return {
                "success": False,
                "error": {"code": "INVALID_AA", "message": f"Unknown amino acid '{aa}'",
                          "available": list(AMINO_ACIDS.keys())}
            }
    
    atoms = []
    n_residues = len(sequence)
    
    if secondary_structure in HELIX_PARAMS:
        params = HELIX_PARAMS[secondary_structure]
        rise = params["rise"]
        rotation = np.radians(params["rotation"])
        radius = params["radius"]
    else:
        rise = 3.5  # Extended
        rotation = np.radians(180)
        radius = 0
    
    # Generate backbone
    for i, aa in enumerate(sequence):
        z = i * rise
        angle = i * rotation
        
        # Backbone N
        x_n = radius * np.cos(angle)
        y_n = radius * np.sin(angle)
        atoms.append({"element": "N", "cartesian": [x_n, y_n, z], "residue": i+1, "aa": aa})
        
        # Backbone CA
        x_ca = radius * np.cos(angle + 0.3)
        y_ca = radius * np.sin(angle + 0.3)
        atoms.append({"element": "C", "cartesian": [x_ca, y_ca, z + 0.5], "residue": i+1, "aa": aa})
        
        # Backbone C
        x_c = radius * np.cos(angle + 0.6)
        y_c = radius * np.sin(angle + 0.6)
        atoms.append({"element": "C", "cartesian": [x_c, y_c, z + 1.0], "residue": i+1, "aa": aa})
        
        # Backbone O
        atoms.append({"element": "O", "cartesian": [x_c + 1.2, y_c, z + 1.0], "residue": i+1, "aa": aa})
        
        # Side chain beta-carbon (simplified)
        if aa != "G":  # Glycine has no side chain
            x_cb = x_ca + 1.5 * np.cos(angle + np.pi/2)
            y_cb = y_ca + 1.5 * np.sin(angle + np.pi/2)
            atoms.append({"element": "C", "cartesian": [x_cb, y_cb, z + 0.5], "residue": i+1, "aa": aa})
        
        # Add backbone H
        atoms.append({"element": "H", "cartesian": [x_n - 0.5, y_n - 0.5, z], "residue": i+1, "aa": aa})
    
    # Calculate charge
    total_charge = sum(AMINO_ACIDS[aa]["charge"] for aa in sequence)
    if termini == "NH2-COOH":
        total_charge += 1  # N-terminus
        total_charge -= 1  # C-terminus (net 0 from termini)
    
    return {
        "success": True,
        "sequence": sequence,
        "n_residues": n_residues,
        "secondary_structure": secondary_structure,
        "termini": termini,
        "total_charge": total_charge,
        "molecular_weight_approx": sum(110 for _ in sequence),  # Approximate MW
        "n_atoms": len(atoms),
        "residue_details": [{"position": i+1, "aa": aa, **AMINO_ACIDS[aa]} for i, aa in enumerate(sequence)],
        "structure": {"atoms": atoms}
    }


def generate_dna_strand(
    sequence: str = "ATCGATCG",
    form: str = "B",
    single_stranded: bool = False
) -> Dict[str, Any]:
    """
    Generate DNA structure.
    
    Args:
        sequence: Nucleotide sequence (5' to 3')
        form: DNA form ('A', 'B', 'Z')
        single_stranded: If True, generate single strand
    
    Returns:
        DNA structure
    """
    # DNA form parameters
    form_params = {
        "A": {"rise": 2.56, "rotation": 32.7, "radius": 10.0, "bp_per_turn": 11},
        "B": {"rise": 3.38, "rotation": 36.0, "radius": 10.0, "bp_per_turn": 10},
        "Z": {"rise": 3.72, "rotation": -60.0, "radius": 9.2, "bp_per_turn": 12},
    }
    
    if form not in form_params:
        return {"success": False, "error": {"code": "INVALID_FORM", "message": f"Unknown form '{form}'"}}
    
    params = form_params[form]
    rise = params["rise"]
    rotation = np.radians(params["rotation"])
    radius = params["radius"]
    
    atoms = []
    n_bases = len(sequence)
    
    for i, base in enumerate(sequence):
        if base not in NUCLEOTIDES:
            return {"success": False, "error": {"code": "INVALID_BASE", "message": f"Unknown base '{base}'"}}
        
        z = i * rise
        angle = i * rotation
        
        # Sugar-phosphate backbone (simplified)
        x_p = radius * np.cos(angle)
        y_p = radius * np.sin(angle)
        atoms.append({"element": "P", "cartesian": [x_p, y_p, z], "strand": "5to3", "base_idx": i})
        
        # Sugar
        x_s = (radius - 1.5) * np.cos(angle)
        y_s = (radius - 1.5) * np.sin(angle)
        atoms.append({"element": "O", "cartesian": [x_s, y_s, z], "strand": "5to3"})
        atoms.append({"element": "C", "cartesian": [x_s, y_s, z + 0.5], "strand": "5to3"})
        
        # Base (simplified as single N at center of base pair)
        x_b = (radius - 5) * np.cos(angle)
        y_b = (radius - 5) * np.sin(angle)
        atoms.append({"element": "N", "cartesian": [x_b, y_b, z], "strand": "5to3", "base": base})
        
        # Complementary strand
        if not single_stranded:
            comp_base = NUCLEOTIDES[base]["pairs_with"]
            angle_comp = angle + np.pi
            
            x_p2 = radius * np.cos(angle_comp)
            y_p2 = radius * np.sin(angle_comp)
            atoms.append({"element": "P", "cartesian": [x_p2, y_p2, z], "strand": "3to5"})
            
            x_b2 = (radius - 5) * np.cos(angle_comp)
            y_b2 = (radius - 5) * np.sin(angle_comp)
            atoms.append({"element": "N", "cartesian": [x_b2, y_b2, z], "strand": "3to5", "base": comp_base})
    
    # Complementary sequence
    comp_seq = "".join(NUCLEOTIDES[b]["pairs_with"] for b in sequence)[::-1]
    
    return {
        "success": True,
        "sequence_5to3": sequence,
        "complement_5to3": comp_seq if not single_stranded else None,
        "n_base_pairs": n_bases,
        "form": form,
        "single_stranded": single_stranded,
        "helix_rise_angstrom": rise,
        "helix_rotation_deg": params["rotation"],
        "bp_per_turn": params["bp_per_turn"],
        "n_atoms": len(atoms),
        "structure": {"atoms": atoms}
    }


def generate_g_quadruplex(
    sequence: str = "GGGTTAGGGTTAGGGTTAGGG",
    topology: str = "parallel",
    cation: str = "K"
) -> Dict[str, Any]:
    """
    Generate G-quadruplex structure.
    
    Args:
        sequence: G-rich sequence
        topology: 'parallel', 'antiparallel', 'hybrid'
        cation: Stabilizing cation (K+, Na+)
    
    Returns:
        G-quadruplex structure
    """
    atoms = []
    
    # G-quadruplex parameters
    g_quartet_radius = 3.5  # G-quartet radius
    quartet_spacing = 3.4   # Stacking distance
    
    # Count potential G-quartets
    g_runs = sequence.count("GGG")
    n_quartets = min(g_runs, 3)
    
    # Generate G-quartets
    for q in range(n_quartets):
        z = q * quartet_spacing
        
        # Four guanines per quartet
        for g in range(4):
            angle = g * np.pi / 2
            
            # Guanine base (simplified)
            x = g_quartet_radius * np.cos(angle)
            y = g_quartet_radius * np.sin(angle)
            
            atoms.append({"element": "N", "cartesian": [x, y, z], "quartet": q+1, "position": g+1})
            atoms.append({"element": "C", "cartesian": [x*0.7, y*0.7, z], "quartet": q+1})
            atoms.append({"element": "O", "cartesian": [x*0.4, y*0.4, z], "quartet": q+1})
        
        # Central cation
        if q < n_quartets - 1:
            z_cat = z + quartet_spacing / 2
            atoms.append({"element": cation, "cartesian": [0, 0, z_cat], "role": "stabilizing_cation"})
    
    return {
        "success": True,
        "sequence": sequence,
        "topology": topology,
        "n_quartets": n_quartets,
        "cation": cation,
        "n_cations": max(0, n_quartets - 1),
        "n_atoms": len(atoms),
        "structure": {"atoms": atoms}
    }


def generate_dipeptide_library(
    output_format: str = "list"
) -> Dict[str, Any]:
    """
    Generate all 400 dipeptide combinations.
    
    Args:
        output_format: 'list' or 'structures'
    
    Returns:
        Dipeptide library information
    """
    amino_acids = list(AMINO_ACIDS.keys())
    dipeptides = []
    
    for aa1 in amino_acids:
        for aa2 in amino_acids:
            seq = aa1 + aa2
            charge = AMINO_ACIDS[aa1]["charge"] + AMINO_ACIDS[aa2]["charge"]
            polar = AMINO_ACIDS[aa1]["polar"] or AMINO_ACIDS[aa2]["polar"]
            
            dipeptides.append({
                "sequence": seq,
                "name_3letter": f"{AMINO_ACIDS[aa1]['3letter']}-{AMINO_ACIDS[aa2]['3letter']}",
                "charge": charge,
                "polar": polar
            })
    
    # Categorize by properties
    cationic = [d for d in dipeptides if d["charge"] > 0]
    anionic = [d for d in dipeptides if d["charge"] < 0]
    neutral = [d for d in dipeptides if d["charge"] == 0]
    
    return {
        "success": True,
        "n_dipeptides": len(dipeptides),
        "n_cationic": len(cationic),
        "n_anionic": len(anionic),
        "n_neutral_polar": len([d for d in neutral if d["polar"]]),
        "n_neutral_hydrophobic": len([d for d in neutral if not d["polar"]]),
        "dipeptides": dipeptides if output_format == "list" else None,
        "sample_cationic": [d["sequence"] for d in cationic[:5]],
        "sample_anionic": [d["sequence"] for d in anionic[:5]],
    }


def generate_tripeptide_library(
    output_format: str = "summary"
) -> Dict[str, Any]:
    """
    Generate tripeptide library statistics.
    
    Args:
        output_format: 'summary' or 'partial'
    
    Returns:
        Tripeptide library information
    """
    amino_acids = list(AMINO_ACIDS.keys())
    n_tripeptides = len(amino_acids) ** 3  # 8000 combinations
    
    # Sample some categories
    samples = {
        "all_hydrophobic": ["ALV", "VLA", "ILV", "LIV", "FLV"],
        "all_charged": ["KRK", "KRD", "DED", "ERE", "KKK"],
        "mixed_polar": ["STS", "TNQ", "QNS", "YST", "HQN"],
        "cell_penetrating_like": ["RRR", "KRK", "RKR", "KKR", "RRK"],
    }
    
    return {
        "success": True,
        "n_tripeptides": n_tripeptides,
        "n_amino_acids": len(amino_acids),
        "sample_sequences": samples,
        "note": "Full library generation available on demand"
    }
