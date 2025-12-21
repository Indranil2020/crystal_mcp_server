"""
generators/__init__.py - Central Registry for All Structure Generators

This module provides a unified registry of all generator functions,
enabling dynamic discovery and routing from the MCP layer.
"""

from typing import Dict, Any, List, Callable

# =============================================================================
# GENERATOR REGISTRY
# =============================================================================
# Each category maps to its available operations with metadata

GENERATOR_REGISTRY: Dict[str, Dict[str, Any]] = {
    # -------------------------------------------------------------------------
    # ADSORPTION (Category)
    # -------------------------------------------------------------------------
    "adsorption": {
        "description": "Adsorption structures",
        "operations": {
            "generate_adsorbate_on_surface": {
                "module": "generators.adsorption.adsorbates",
                "function": "generate_adsorbate_on_surface",
                "params": ['surface_structure', 'adsorbate', 'site', 'coverage', 'height_offset', 'rotation_deg'],
                "description": "Place adsorbate molecule on surface."
            },
            "add_adsorbate": {
                "module": "generators.adsorption.adsorbates",
                "function": "generate_adsorbate_on_surface",
                "params": ['surface_structure', 'adsorbate', 'site', 'coverage', 'height_offset', 'rotation_deg'],
                "description": "Alias for generate_adsorbate_on_surface."
            },
            "generate_coadsorption": {
                "module": "generators.adsorption.adsorbates",
                "function": "generate_coadsorption",
                "params": ['surface_structure', 'adsorbates', 'sites', 'positions'],
                "description": "Generate co-adsorption configuration."
            },
            "generate_confined_il": {
                "module": "generators.adsorption.ionic_liquids",
                "function": "generate_confined_il",
                "params": ['ionic_liquid', 'confinement', 'pore_diameter_A'],
                "description": "Generate confined ionic liquid."
            },
            "generate_coverage_pattern": {
                "module": "generators.adsorption.fragments",
                "function": "generate_coverage_pattern",
                "params": ['surface', 'adsorbate', 'pattern', 'height'],
                "description": "Generate ordered adsorbate coverage pattern."
            },
            "generate_electrode_il_interface": {
                "module": "generators.adsorption.ionic_liquids",
                "function": "generate_electrode_il_interface",
                "params": ['electrode', 'ionic_liquid', 'electrode_layers', 'il_thickness_A', 'n_ion_pairs'],
                "description": "Generate electrode-ionic liquid interface."
            },
            "generate_fragment_adsorbate": {
                "module": "generators.adsorption.fragments",
                "function": "generate_fragment_adsorbate",
                "params": ['surface', 'fragment', 'site', 'coverage', 'site_index'],
                "description": "Add dissociated fragment to surface."
            },
            "generate_ice_layer": {
                "module": "generators.adsorption.water_layers",
                "function": "generate_ice_layer",
                "params": ['surface', 'ice_type', 'n_bilayers', 'proton_ordering'],
                "description": "Generate ice layer on surface."
            },
            "generate_ionic_liquid_bulk": {
                "module": "generators.adsorption.ionic_liquids",
                "function": "generate_ionic_liquid_bulk",
                "params": ['ionic_liquid', 'box_size_A', 'n_ion_pairs'],
                "description": "Generate bulk ionic liquid structure."
            },
            "generate_mixed_sam": {
                "module": "generators.adsorption.sams",
                "function": "generate_mixed_sam",
                "params": ['sam1', 'sam2', 'ratio', 'mixing_pattern', 'size'],
                "description": "Generate mixed SAM with two components."
            },
            "generate_sac_graphene": {
                "module": "generators.adsorption.sac",
                "function": "generate_sac_graphene",
                "params": ['metal', 'coordination', 'supercell', 'n_sac_sites'],
                "description": "Generate SAC on graphene-based support."
            },
            "generate_sac_oxide": {
                "module": "generators.adsorption.sac",
                "function": "generate_sac_oxide",
                "params": ['metal', 'support', 'supercell', 'anchoring_site'],
                "description": "Generate SAC on oxide support."
            },
            "generate_sac_tmd": {
                "module": "generators.adsorption.sac",
                "function": "generate_sac_tmd",
                "params": ['metal', 'support', 'defect_type', 'supercell'],
                "description": "Generate SAC on 2D TMD support."
            },
            "generate_sam": {
                "module": "generators.adsorption.sams",
                "function": "generate_sam",
                "params": ['sam_type', 'size', 'defect_density'],
                "description": "Generate self-assembled monolayer."
            },
            "generate_water_bilayer": {
                "module": "generators.adsorption.water_layers",
                "function": "generate_water_bilayer",
                "params": ['surface', 'orientation', 'coverage', 'dissociation_fraction'],
                "description": "Generate water bilayer on surface."
            },
        }
    },
    
    # -------------------------------------------------------------------------
    # BATTERY (Category)
    # -------------------------------------------------------------------------
    "battery": {
        "description": "Battery structures",
        "operations": {
            "generate_cathode": {
                "module": "generators.battery.cathodes",
                "function": "generate_cathode",
                "params": ['cathode_type', 'lithiation', 'supercell'],
                "description": "Generate cathode structure."
            },
            "generate_cathode_electrolyte_interface": {
                "module": "generators.battery.electrolytes",
                "function": "generate_cathode_electrolyte_interface",
                "params": ['cathode', 'electrolyte', 'interface_orientation'],
                "description": "Generate cathode-electrolyte interface."
            },
            "generate_electrolyte": {
                "module": "generators.battery.electrolytes",
                "function": "generate_electrolyte",
                "params": ['electrolyte_type', 'supercell', 'li_site_disorder'],
                "description": "Generate solid electrolyte structure."
            },
            "generate_lithiated_series": {
                "module": "generators.battery.cathodes",
                "function": "generate_lithiated_series",
                "params": ['cathode_type', 'n_points'],
                "description": "Generate series of lithiated structures."
            },
            "generate_surface_coated_cathode": {
                "module": "generators.battery.cathodes",
                "function": "generate_surface_coated_cathode",
                "params": ['cathode_type', 'coating', 'coating_thickness_nm', 'supercell'],
                "description": "Generate surface-coated cathode."
            },
            # New: Anode materials
            "generate_silicon_anode": {
                "module": "generators.battery.anodes",
                "function": "generate_silicon_anode",
                "params": ['phase', 'lithiation', 'supercell', 'amorphous', 'seed'],
                "description": "Generate silicon anode structure with Li insertion."
            },
            "generate_germanium_anode": {
                "module": "generators.battery.anodes",
                "function": "generate_germanium_anode",
                "params": ['phase', 'lithiation', 'supercell'],
                "description": "Generate germanium anode structure."
            },
            "generate_tin_anode": {
                "module": "generators.battery.anodes",
                "function": "generate_tin_anode",
                "params": ['phase', 'lithiation', 'supercell'],
                "description": "Generate tin anode structure."
            },
            "generate_phosphorus_anode": {
                "module": "generators.battery.anodes",
                "function": "generate_phosphorus_anode",
                "params": ['phase', 'lithiation', 'supercell'],
                "description": "Generate phosphorus anode (black P, red P)."
            },
            "generate_sulfur_cathode": {
                "module": "generators.battery.anodes",
                "function": "generate_sulfur_cathode",
                "params": ['phase', 'lithiation', 'host_material', 'supercell'],
                "description": "Generate sulfur cathode for Li-S batteries."
            },
            "generate_metal_anode": {
                "module": "generators.battery.anodes",
                "function": "generate_metal_anode",
                "params": ['metal', 'surface', 'thickness_layers', 'include_dendrite', 'include_sei'],
                "description": "Generate metal anode (Li, Na, K, Zn) with optional dendrite/SEI."
            },
        }
    },

    # -------------------------------------------------------------------------
    # BIOMOLECULAR (Category)
    # -------------------------------------------------------------------------
    "biomolecular": {
        "description": "Biomolecular structures",
        "operations": {
            "generate_dna_surface": {
                "module": "generators.biomolecular.interfaces",
                "function": "generate_dna_surface",
                "params": ['surface', 'dna_type', 'surface_size'],
                "description": "Generate DNA origami on 2D surface."
            },
            "generate_lipid_bilayer": {
                "module": "generators.biomolecular.lipids",
                "function": "generate_lipid_bilayer",
                "params": ['lipid', 'size', 'cholesterol_fraction', 'water_buffer', 'seed'],
                "description": "Generate lipid bilayer structure."
            },
            "generate_membrane_with_protein": {
                "module": "generators.biomolecular.lipids",
                "function": "generate_membrane_with_protein",
                "params": ['lipid', 'membrane_size', 'protein_type', 'seed'],
                "description": "Generate membrane with embedded protein."
            },
            "generate_peptide_surface": {
                "module": "generators.biomolecular.interfaces",
                "function": "generate_peptide_surface",
                "params": ['peptide', 'surface', 'surface_size', 'coverage_nmol_nm2'],
                "description": "Generate peptide-semiconductor interface."
            },
            "generate_protein_surface": {
                "module": "generators.biomolecular.interfaces",
                "function": "generate_protein_surface",
                "params": ['protein', 'surface', 'surface_size', 'orientation', 'coverage'],
                "description": "Generate protein-surface interface."
            },
        }
    },
    
    # -------------------------------------------------------------------------
    # BULK (Category)
    # -------------------------------------------------------------------------
    "bulk": {
        "description": "Bulk structures",
        "operations": {
            "apply_strain": {
                "module": "generators.bulk.strain",
                "function": "apply_strain",
                "params": ['structure_dict', 'strain', 'strain_type'],
                "description": "Apply strain to a structure."
            },
            "apply_strain_tensor": {
                "module": "generators.bulk.strain",
                "function": "apply_strain_tensor",
                "params": ['structure_dict', 'strain_tensor'],
                "description": "Apply full 3x3 strain tensor to structure."
            },
            "generate_all_cubic_prototypes": {
                "module": "generators.bulk.spacegroups",
                "function": "generate_all_cubic_prototypes",
                "params": ['element'],
                "description": "Generate all common cubic structure types for an element."
            },
            "generate_antiphase_boundary": {
                "module": "generators.bulk.boundaries",
                "function": "generate_antiphase_boundary",
                "params": ['structure_dict', 'apb_plane', 'shift_vector'],
                "description": "Generate an antiphase boundary (APB)."
            },
            "generate_antisite": {
                "module": "generators.bulk.defects",
                "function": "generate_antisite",
                "params": ['structure_dict', 'site_a', 'site_b'],
                "description": "Generate antisite defect by swapping two atoms."
            },
            "generate_bain_path": {
                "module": "generators.bulk.strain",
                "function": "generate_bain_path",
                "params": ['element', 'n_steps', 'a_bcc'],
                "description": "Generate structures along Bain transformation path (BCC ↔ FCC)."
            },
            "generate_breathing_mode": {
                "module": "generators.bulk.quantum_phases",
                "function": "generate_breathing_mode",
                "params": ['structure_dict', 'amplitude', 'pattern'],
                "description": "Generate breathing mode distortion (e.g., breathing Kagome)."
            },
            "generate_cdw_phase": {
                "module": "generators.bulk.quantum_phases",
                "function": "generate_cdw_phase",
                "params": ['material', 'phase', 'temperature'],
                "description": "Generate charge density wave structure."
            },
            "generate_clathrate": {
                "module": "generators.bulk.clathrates",
                "function": "generate_clathrate",
                "params": ['clathrate', 'empty_cages', 'supercell'],
                "description": "Generate clathrate structure."
            },
            "generate_cuprate": {
                "module": "generators.bulk.cuprates",
                "function": "generate_cuprate",
                "params": ['cuprate_type', 'n_cuo2_planes', 'oxygen_content', 'supercell'],
                "description": "Generate cuprate superconductor structure."
            },
            "generate_custom_cuprate": {
                "module": "generators.bulk.cuprates",
                "function": "generate_custom_cuprate",
                "params": ['n_cuo2_planes', 'charge_reservoir', 'spacer', 'a', 'c_per_layer'],
                "description": "Generate custom cuprate with specified layer sequence."
            },
            "generate_custom_stacking": {
                "module": "generators.bulk.polytypes",
                "function": "generate_custom_stacking",
                "params": ['element', 'sequence', 'a', 'c_layer', 'bilayer_element'],
                "description": "Generate a structure with custom stacking sequence."
            },
            "generate_dislocation_dipole": {
                "module": "generators.bulk.dislocations",
                "function": "generate_dislocation_dipole",
                "params": ['structure_dict', 'dislocation_type', 'separation', 'supercell'],
                "description": "Generate dislocation dipole structure."
            },
            "generate_dopant_distribution": {
                "module": "generators.bulk.defects",
                "function": "generate_dopant_distribution",
                "params": ['structure_dict', 'host_element', 'dopant_element', 'concentration', 'distribution', 'seed'],
                "description": "Generate a doped structure with specified dopant distribution."
            },
            "generate_double_perovskite": {
                "module": "generators.bulk.perovskites",
                "function": "generate_double_perovskite",
                "params": ['a_site', 'b_site_1', 'b_site_2', 'x_site', 'ordering', 'a'],
                "description": "Generate double perovskite A2BB'X6 structure."
            },
            "generate_elastic_dipole": {
                "module": "generators.bulk.dislocations",
                "function": "generate_elastic_dipole",
                "params": ['structure_dict', 'dipole_type', 'dipole_tensor'],
                "description": "Generate elastic dipole structure."
            },
            "generate_empty_clathrate": {
                "module": "generators.bulk.clathrates",
                "function": "generate_empty_clathrate",
                "params": ['framework', 'clathrate_type'],
                "description": "Generate empty clathrate framework (superconducting Si-46)."
            },
            "generate_energetic_cocrystal": {
                "module": "generators.bulk.energetics",
                "function": "generate_energetic_cocrystal",
                "params": ['component_a', 'component_b', 'ratio'],
                "description": "Generate energetic cocrystal structure."
            },
            "generate_energetic_material": {
                "module": "generators.bulk.energetics",
                "function": "generate_energetic_material",
                "params": ['material', 'polymorph'],
                "description": "Generate energetic material crystal structure."
            },
            "generate_epitaxial_strain": {
                "module": "generators.bulk.strain",
                "function": "generate_epitaxial_strain",
                "params": ['structure_dict', 'substrate_a', 'substrate_b', 'relax_c', 'poisson_ratio'],
                "description": "Apply epitaxial strain from substrate."
            },
            "generate_frank_loop": {
                "module": "generators.bulk.dislocations",
                "function": "generate_frank_loop",
                "params": ['structure_dict', 'loop_type', 'radius', 'habit_plane'],
                "description": "Generate Frank partial dislocation loop."
            },
            "generate_from_spacegroup": {
                "module": "generators.bulk.spacegroups",
                "function": "generate_from_spacegroup",
                "params": ['spacegroup', 'elements', 'composition', 'a', 'b', 'c', 'alpha', 'beta', 'gamma', 'factor'],
                "description": "Generate structure in any of 230 space groups using PyXtal."
            },
            "generate_from_spacegroup_symbol": {
                "module": "generators.bulk.spacegroups",
                "function": "generate_from_spacegroup_symbol",
                "params": ['symbol', 'elements', 'composition', 'factor'],
                "description": "Generate structure from space group symbol (e.g., 'Fm-3m', 'P6_3/mmc')."
            },
            "generate_graded_supercell": {
                "module": "generators.bulk.supercells",
                "function": "generate_graded_supercell",
                "params": ['base_structure', 'element_from', 'element_to', 'gradient_axis', 'supercell'],
                "description": "Generate compositionally graded supercell."
            },
            "generate_high_pressure_phase": {
                "module": "generators.bulk.special_phases",
                "function": "generate_high_pressure_phase",
                "params": ['phase_name'],
                "description": "Generate a known high-pressure phase."
            },
            # Note: generate_interstitial is in defect category, not bulk
            "generate_jahn_teller": {
                "module": "generators.bulk.quantum_phases",
                "function": "generate_jahn_teller",
                "params": ['structure_dict', 'jt_metal', 'jt_type', 'elongation'],
                "description": "Generate Jahn-Teller distorted structure."
            },
            "generate_layered_perovskite": {
                "module": "generators.bulk.perovskites",
                "function": "generate_layered_perovskite",
                "params": ['a_site', 'a_prime_site', 'b_site', 'x_site', 'n_layers', 'family'],
                "description": "Generate layered perovskite structure."
            },
            "generate_magnetic_ordering": {
                "module": "generators.bulk.magnetic",
                "function": "generate_magnetic_ordering",
                "params": ['structure_dict', 'magnetic_elements', 'ordering', 'moments', 'axis', 'q_vector'],
                "description": "Generate magnetic ordering on a structure."
            },
            "generate_nte_material": {
                "module": "generators.bulk.special_phases",
                "function": "generate_nte_material",
                "params": ['material'],
                "description": "Generate negative thermal expansion material."
            },
            "generate_pbx_model": {
                "module": "generators.bulk.energetics",
                "function": "generate_pbx_model",
                "params": ['energetic', 'binder', 'energetic_fraction'],
                "description": "Generate plastic-bonded explosive (PBX) model."
            },
            "generate_peierls_distortion": {
                "module": "generators.bulk.quantum_phases",
                "function": "generate_peierls_distortion",
                "params": ['structure_dict', 'chain_direction', 'dimerization'],
                "description": "Generate Peierls distorted 1D chain structure."
            },
            "generate_perovskite": {
                "module": "generators.bulk.perovskites",
                "function": "generate_perovskite",
                "params": ['a_site', 'b_site', 'x_site', 'a', 'phase', 'tilting'],
                "description": "Generate simple perovskite ABX3 structure."
            },
            "generate_ruddlesden_popper": {
                "module": "generators.bulk.perovskites",
                "function": "generate_ruddlesden_popper",
                "params": ['compound', 'n_layers'],
                "description": "Generate Ruddlesden-Popper layered perovskite (e.g., Sr2TiO4)."
            },
            "generate_polytype": {
                "module": "generators.bulk.polytypes",
                "function": "generate_polytype",
                "params": ['element', 'polytype', 'a', 'c_layer', 'bilayer_element', 'vacuum'],
                "description": "Generate a polytype structure with specified stacking sequence."
            },
            "generate_porous_energetic": {
                "module": "generators.bulk.energetics",
                "function": "generate_porous_energetic",
                "params": ['material', 'void_fraction', 'void_shape', 'void_radius'],
                "description": "Generate porous energetic material (e.g., for PBX modeling)."
            },
            "generate_prismatic_loop": {
                "module": "generators.bulk.dislocations",
                "function": "generate_prismatic_loop",
                "params": ['structure_dict', 'loop_type', 'radius', 'orientation'],
                "description": "Generate prismatic dislocation loop."
            },
            "generate_prototype": {
                "module": "generators.bulk.spacegroups",
                "function": "generate_prototype",
                "params": ['prototype', 'elements', 'a', 'c'],
                "description": "Generate common prototype structure."
            },
            "generate_prussian_blue": {
                "module": "generators.bulk.special_phases",
                "function": "generate_prussian_blue",
                "params": ['metal_a', 'metal_b', 'a'],
                "description": "Generate Prussian blue analog structure."
            },
            "generate_random_alloy": {
                "module": "generators.bulk.disorder",
                "function": "generate_random_alloy",
                "params": ['structure_dict', 'substitutions', 'seed'],
                "description": "Generate random alloy by random substitution."
            },
            "generate_random_alloy_supercell": {
                "module": "generators.bulk.supercells",
                "function": "generate_random_alloy_supercell",
                "params": ['base_structure', 'substitutions', 'supercell', 'seed'],
                "description": "Generate random alloy supercell."
            },
            "generate_sft": {
                "module": "generators.bulk.dislocations",
                "function": "generate_sft",
                "params": ['structure_dict', 'n_vacancies'],
                "description": "Generate stacking fault tetrahedron (SFT)."
            },
            "generate_slab_supercell": {
                "module": "generators.bulk.supercells",
                "function": "generate_slab_supercell",
                "params": ['base_structure', 'miller_index', 'slab_thickness', 'vacuum'],
                "description": "Generate slab supercell for surface calculations."
            },
            "generate_spin_peierls": {
                "module": "generators.bulk.quantum_phases",
                "function": "generate_spin_peierls",
                "params": ['structure_dict', 'magnetic_element', 'dimerization'],
                "description": "Generate spin-Peierls distorted structure."
            },
            "generate_sqs": {
                "module": "generators.bulk.disorder",
                "function": "generate_sqs",
                "params": ['structure_dict', 'target_composition', 'n_iterations', 'correlation_cutoff', 'seed'],
                "description": "Generate Special Quasirandom Structure (SQS) using simulated annealing."
            },
            "generate_stacking_fault": {
                "module": "generators.bulk.boundaries",
                "function": "generate_stacking_fault",
                "params": ['structure_dict', 'fault_type', 'fault_plane', 'position'],
                "description": "Generate a stacking fault structure."
            },
            "generate_superhydride": {
                "module": "generators.bulk.special_phases",
                "function": "generate_superhydride",
                "params": ['metal', 'h_content', 'a'],
                "description": "Generate superhydride structure (MHn)."
            },
            "generate_twin_boundary": {
                "module": "generators.bulk.boundaries",
                "function": "generate_twin_boundary",
                "params": ['structure_dict', 'twin_plane', 'n_layers_below', 'n_layers_above', 'twin_type'],
                "description": "Generate a twin boundary structure."
            },
            # Note: generate_vacancy is in defect category, not bulk
            "generate_zeolite": {
                "module": "generators.bulk.zeolites",
                "function": "generate_zeolite",
                "params": ['framework', 'cation', 'Si_Al_ratio', 'supercell'],
                "description": "Generate zeolite structure."
            },
            "generate_zeolite_with_guest": {
                "module": "generators.bulk.zeolites",
                "function": "generate_zeolite_with_guest",
                "params": ['framework', 'guest_molecule', 'loading'],
                "description": "Generate zeolite with guest molecules."
            },
            "make_supercell": {
                "module": "generators.bulk.supercells",
                "function": "make_supercell",
                "params": ['structure', 'scaling'],
                "description": "Create supercell from structure."
            },
        }
    },
    
    # -------------------------------------------------------------------------
    # CATALYST (Category)
    # -------------------------------------------------------------------------
    "catalyst": {
        "description": "Catalyst structures",
        "operations": {
            "generate_cluster": {
                "module": "generators.catalyst.clusters",
                "function": "generate_cluster",
                "params": ['cluster_type', 'return_atoms_only'],
                "description": "Generate isolated metal cluster."
            },
            "generate_hea_surface": {
                "module": "generators.catalyst.saa",
                "function": "generate_hea_surface",
                "params": ['hea_type', 'miller_index', 'supercell', 'termination', 'surface_segregation'],
                "description": "Generate High-Entropy Alloy surface."
            },
            "generate_photocatalyst_heterojunction": {
                "module": "generators.catalyst.clusters",
                "function": "generate_photocatalyst_heterojunction",
                "params": ['junction_type', 'materials', 'interface_orientation'],
                "description": "Generate photocatalyst heterojunction."
            },
            "generate_saa_surface": {
                "module": "generators.catalyst.saa",
                "function": "generate_saa_surface",
                "params": ['saa_type', 'miller_index', 'supercell', 'n_dopant', 'dopant_layer', 'vacuum'],
                "description": "Generate Single Atom Alloy surface."
            },
            "generate_sac_2d_support": {
                "module": "generators.catalyst.saa",
                "function": "generate_sac_2d_support",
                "params": ['metal', 'coordination_motif', 'support', 'supercell', 'n_sac_sites'],
                "description": "Generate 2D SAC with specific coordination motif."
            },
            "generate_strained_catalyst": {
                "module": "generators.catalyst.saa",
                "function": "generate_strained_catalyst",
                "params": ['catalyst_type', 'strain_type', 'strain_magnitude', 'supercell'],
                "description": "Generate strain-engineered catalyst surface."
            },
            "generate_supported_cluster": {
                "module": "generators.catalyst.clusters",
                "function": "generate_supported_cluster",
                "params": ['cluster', 'support', 'surface', 'supercell', 'adsorption_site'],
                "description": "Generate oxide-supported metal cluster."
            },
            "generate_zeolite_cluster": {
                "module": "generators.catalyst.clusters",
                "function": "generate_zeolite_cluster",
                "params": ['metal', 'zeolite', 'cluster_size', 'location'],
                "description": "Generate zeolite-confined metal cluster."
            },
        }
    },
    
    # -------------------------------------------------------------------------
    # DEFECT (Category)
    # -------------------------------------------------------------------------
    "defect": {
        "description": "Defect structures",
        "operations": {
            "generate_amorphous_interface": {
                "module": "generators.defect.amorphous",
                "function": "generate_amorphous_interface",
                "params": ['crystalline', 'amorphous', 'crystal_layers', 'amorphous_thickness_A'],
                "description": "Generate crystalline-amorphous interface."
            },
            "generate_amorphous_structure": {
                "module": "generators.defect.amorphous",
                "function": "generate_amorphous_structure",
                "params": ['material', 'n_atoms', 'density_gcm3', 'seed'],
                "description": "Generate amorphous structure using random placement with constraints."
            },
            "generate_antisite": {
                "module": "generators.defect.point_defects",
                "function": "generate_antisite",
                "params": ['host_structure', 'site1', 'site2'],
                "description": "Generate antisite defect (swap two atoms)."
            },
            "generate_crack_tip": {
                "module": "generators.defect.volume_defects",
                "function": "generate_crack_tip",
                "params": ['element', 'crack_mode', 'crack_length', 'supercell_size'],
                "description": "Generate structure with crack tip."
            },
            "create_defect": {
                "module": "generators.defect.point_defects",
                "function": "create_defect",
                "params": ['structure', 'defect_type', 'site_index', 'dopant', 'interstitial_species', 'interstitial_type'],
                "description": "Create a defect in a structure (general wrapper for vacancy, interstitial, substitution, etc.)."
            },
            "generate_dislocation": {
                "module": "generators.defect.extended_defects",
                "function": "generate_dislocation",
                "params": ['element', 'dislocation_type', 'supercell'],
                "description": "Generate structure with dislocation."
            },
            "generate_f_center": {
                "module": "generators.defect.point_defects",
                "function": "generate_f_center",
                "params": ['host_structure', 'anion_site'],
                "description": "Generate F-center (color center) by removing anion."
            },
            "generate_frenkel_pair": {
                "module": "generators.defect.point_defects",
                "function": "generate_frenkel_pair",
                "params": ['host_structure', 'vacancy_site', 'interstitial_position', 'element'],
                "description": "Generate Frenkel defect pair (vacancy + interstitial)."
            },
            "generate_grain_boundary": {
                "module": "generators.defect.extended_defects",
                "function": "generate_grain_boundary",
                "params": ['element', 'boundary_type', 'supercell'],
                "description": "Generate grain boundary structure."
            },
            "generate_hea_sqs": {
                "module": "generators.defect.sqs",
                "function": "generate_hea_sqs",
                "params": ['elements', 'composition', 'supercell_size', 'seed'],
                "description": "Generate High-Entropy Alloy SQS structure."
            },
            "generate_interstitial": {
                "module": "generators.defect.point_defects",
                "function": "generate_interstitial",
                "params": ['host_structure', 'interstitial_species', 'interstitial_type', 'position'],
                "description": "Generate interstitial defect."
            },
            "generate_melt_quench": {
                "module": "generators.defect.amorphous",
                "function": "generate_melt_quench",
                "params": ['elements', 'composition', 'n_atoms', 'quench_rate'],
                "description": "Generate amorphous structure using simulated melt-quench."
            },
            "generate_precipitate": {
                "module": "generators.defect.volume_defects",
                "function": "generate_precipitate",
                "params": ['matrix_element', 'precipitate_element', 'precipitate_type', 'precipitate_radius', 'supercell_size'],
                "description": "Generate structure with precipitate."
            },
            "generate_radiation_cascade": {
                "module": "generators.defect.volume_defects",
                "function": "generate_radiation_cascade",
                "params": ['element', 'pka_energy_keV', 'supercell_size', 'seed'],
                "description": "Generate radiation damage cascade structure."
            },
            "generate_sqs": {
                "module": "generators.defect.sqs",
                "function": "generate_sqs",
                "params": ['alloy', 'composition', 'supercell_size', 'n_iterations', 'seed'],
                "description": "Generate Special Quasirandom Structure."
            },
            "generate_sro_structure": {
                "module": "generators.defect.sqs",
                "function": "generate_sro_structure",
                "params": ['alloy', 'sro_parameter', 'supercell_size', 'seed'],
                "description": "Generate structure with specified short-range order."
            },
            "generate_stacking_fault": {
                "module": "generators.defect.extended_defects",
                "function": "generate_stacking_fault",
                "params": ['element', 'fault_type', 'fault_position', 'supercell'],
                "description": "Generate structure with stacking fault."
            },
            "generate_substitution": {
                "module": "generators.defect.point_defects",
                "function": "generate_substitution",
                "params": ['host_structure', 'site', 'dopant', 'dopant_type'],
                "description": "Generate substitutional defect (dopant)."
            },
            "generate_vacancy": {
                "module": "generators.defect.point_defects",
                "function": "generate_vacancy",
                "params": ['host_structure', 'vacancy_site', 'vacancy_type'],
                "description": "Generate vacancy defect."
            },
            "generate_void": {
                "module": "generators.defect.volume_defects",
                "function": "generate_void",
                "params": ['element', 'void_shape', 'void_radius', 'supercell_size', 'void_position'],
                "description": "Generate structure with void."
            },
        }
    },
    
    # -------------------------------------------------------------------------
    # ELECTRONIC (Category)
    # -------------------------------------------------------------------------
    "electronic": {
        "description": "Electronic structures",
        "operations": {
            "generate_2deg_interface": {
                "module": "generators.electronic.semiconductors",
                "function": "generate_2deg_interface",
                "params": ['system', 'spacer_nm', 'doping_nm', 'doping_density_cm2'],
                "description": "Generate 2DEG heterostructure interface."
            },
            "generate_cuprate_superconductor": {
                "module": "generators.electronic.superconductors",
                "function": "generate_cuprate_superconductor",
                "params": ['base', 'doping', 'n_CuO_planes', 'supercell'],
                "description": "Generate cuprate superconductor with doping control."
            },
            "generate_majorana_nanowire": {
                "module": "generators.electronic.topological",
                "function": "generate_majorana_nanowire",
                "params": ['platform', 'wire_length_nm', 'wire_diameter_nm', 'al_shell_thickness_nm', 'tunnel_barrier'],
                "description": "Generate Majorana nanowire structure."
            },
            "generate_quantum_well": {
                "module": "generators.electronic.semiconductors",
                "function": "generate_quantum_well",
                "params": ['well_material', 'barrier_material', 'well_thickness_nm', 'barrier_thickness_nm', 'n_wells', 'alloy_fraction'],
                "description": "Generate quantum well heterostructure."
            },
            "generate_superconductor": {
                "module": "generators.electronic.superconductors",
                "function": "generate_superconductor",
                "params": ['material', 'supercell'],
                "description": "Generate superconductor structure."
            },
            "generate_superlattice": {
                "module": "generators.electronic.semiconductors",
                "function": "generate_superlattice",
                "params": ['material1', 'material2', 'period1_nm', 'period2_nm', 'n_periods'],
                "description": "Generate superlattice structure."
            },
            "generate_topological_insulator": {
                "module": "generators.electronic.topological",
                "function": "generate_topological_insulator",
                "params": ['material', 'thickness_QL', 'supercell'],
                "description": "Generate topological insulator structure."
            },
            "generate_weyl_semimetal": {
                "module": "generators.electronic.topological",
                "function": "generate_weyl_semimetal",
                "params": ['material', 'supercell'],
                "description": "Generate Weyl or Dirac semimetal structure."
            },
        }
    },
    
    # -------------------------------------------------------------------------
    # EXTERNAL_FIELDS (Category)
    # -------------------------------------------------------------------------
    "external_fields": {
        "description": "External Fields structures",
        "operations": {
            "apply_strain_field": {
                "module": "generators.external_fields.strain_fields",
                "function": "apply_strain_field",
                "params": ['structure_dict', 'strain_profile', 'max_strain', 'direction'],
                "description": "Apply strain field to existing structure."
            },
            "generate_antiskyrmion": {
                "module": "generators.external_fields.magnetic_textures",
                "function": "generate_antiskyrmion",
                "params": ['material', 'size_nm', 'supercell_size'],
                "description": "Generate antiskyrmion structure."
            },
            "generate_magnetic_bobber": {
                "module": "generators.external_fields.magnetic_textures",
                "function": "generate_magnetic_bobber",
                "params": ['material', 'depth_nm', 'supercell_size'],
                "description": "Generate magnetic bobber (truncated skyrmion tube)."
            },
            "generate_rippled_structure": {
                "module": "generators.external_fields.strain_fields",
                "function": "generate_rippled_structure",
                "params": ['material', 'amplitude_angstrom', 'wavelength_angstrom', 'direction', 'size'],
                "description": "Generate rippled 2D structure."
            },
            "generate_skyrmion": {
                "module": "generators.external_fields.magnetic_textures",
                "function": "generate_skyrmion",
                "params": ['material', 'skyrmion_radius_nm', 'supercell_size', 'helicity', 'polarity'],
                "description": "Generate structure with skyrmion spin texture."
            },
            "generate_strained_structure": {
                "module": "generators.external_fields.strain_fields",
                "function": "generate_strained_structure",
                "params": ['element', 'strain_xx', 'strain_yy', 'strain_zz', 'strain_type', 'supercell'],
                "description": "Generate uniformly strained structure."
            },
            # New: Electric field effects
            "generate_polarized_cell": {
                "module": "generators.external_fields.electric_fields",
                "function": "generate_polarized_cell",
                "params": ['material', 'polarization_direction', 'polarization_magnitude', 'supercell'],
                "description": "Generate cell with electric polarization (ferroelectric, piezoelectric)."
            },
            "generate_ferroelectric_domain": {
                "module": "generators.external_fields.electric_fields",
                "function": "generate_ferroelectric_domain",
                "params": ['material', 'domain_type', 'domain_width_nm', 'supercell'],
                "description": "Generate ferroelectric domain structure (180°, 90° walls)."
            },
            "generate_polar_surface": {
                "module": "generators.external_fields.electric_fields",
                "function": "generate_polar_surface",
                "params": ['material', 'surface', 'polarity', 'thickness', 'vacuum', 'compensated'],
                "description": "Generate polar surface with optional compensation."
            },
            "generate_field_induced_phase": {
                "module": "generators.external_fields.electric_fields",
                "function": "generate_field_induced_phase",
                "params": ['material', 'field_strength_V_nm', 'temperature_K'],
                "description": "Generate field-induced phase structure (e.g., polarized SrTiO3)."
            },
            # New: Advanced field effects
            "generate_floquet_cell": {
                "module": "generators.external_fields.advanced_fields",
                "function": "generate_floquet_cell",
                "params": ['material', 'driving_frequency_THz', 'driving_amplitude', 'polarization', 'supercell'],
                "description": "Generate cell for Floquet engineering (periodic driving)."
            },
            "generate_phonon_pumped_cell": {
                "module": "generators.external_fields.advanced_fields",
                "function": "generate_phonon_pumped_cell",
                "params": ['material', 'phonon_mode', 'amplitude_A', 'supercell'],
                "description": "Generate cell with optically pumped phonon mode."
            },
            "generate_warm_dense_matter": {
                "module": "generators.external_fields.advanced_fields",
                "function": "generate_warm_dense_matter",
                "params": ['element', 'temperature_eV', 'density_g_cm3', 'n_atoms'],
                "description": "Generate warm dense matter configuration (high-T plasma)."
            },
            "generate_laser_shocked_cell": {
                "module": "generators.external_fields.advanced_fields",
                "function": "generate_laser_shocked_cell",
                "params": ['material', 'pressure_GPa', 'temperature_K', 'strain_rate', 'supercell'],
                "description": "Generate laser-shocked cell with compression."
            },
            "generate_bimeron": {
                "module": "generators.external_fields.advanced_fields",
                "function": "generate_bimeron",
                "params": ['material', 'bimeron_size_nm', 'supercell_size'],
                "description": "Generate bimeron magnetic texture (in-plane analog of skyrmion)."
            },
        }
    },

    # -------------------------------------------------------------------------
    # HIGH_PRESSURE (Category)
    # -------------------------------------------------------------------------
    "high_pressure": {
        "description": "High Pressure structures",
        "operations": {
            "generate_high_pressure_phase": {
                "module": "generators.high_pressure.phases",
                "function": "generate_high_pressure_phase",
                "params": ['phase', 'supercell'],
                "description": "Generate high-pressure phase structure."
            },
            "generate_polymeric_nitrogen": {
                "module": "generators.high_pressure.phases",
                "function": "generate_polymeric_nitrogen",
                "params": ['phase', 'supercell'],
                "description": "Generate polymeric nitrogen structure."
            },
            "generate_super_ionic_phase": {
                "module": "generators.high_pressure.phases",
                "function": "generate_super_ionic_phase",
                "params": ['compound', 'phase', 'temperature_K'],
                "description": "Generate super-ionic phase structure."
            },
            "generate_superhydride": {
                "module": "generators.high_pressure.phases",
                "function": "generate_superhydride",
                "params": ['compound', 'pressure_GPa', 'supercell'],
                "description": "Generate superhydride structure."
            },
        }
    },
    
    # -------------------------------------------------------------------------
    # MAGNETIC (Category)
    # -------------------------------------------------------------------------
    "magnetic": {
        "description": "Magnetic structures",
        "operations": {
            "generate_antiferromagnet": {
                "module": "generators.magnetic.materials",
                "function": "generate_antiferromagnet",
                "params": ['material', 'ordering', 'supercell'],
                "description": "Generate antiferromagnetic structure."
            },
            "generate_full_heusler": {
                "module": "generators.magnetic.heusler",
                "function": "generate_full_heusler",
                "params": ['formula', 'supercell'],
                "description": "Generate full Heusler alloy (X2YZ, L21 structure)."
            },
            "generate_half_heusler": {
                "module": "generators.magnetic.heusler",
                "function": "generate_half_heusler",
                "params": ['formula', 'supercell'],
                "description": "Generate half Heusler alloy (XYZ, C1b structure)."
            },
            "generate_heusler": {
                "module": "generators.magnetic.materials",
                "function": "generate_heusler",
                "params": ['compound', 'supercell'],
                "description": "Generate Heusler alloy structure (full or half)."
            },
            "generate_inverse_heusler": {
                "module": "generators.magnetic.heusler",
                "function": "generate_inverse_heusler",
                "params": ['formula', 'supercell'],
                "description": "Generate inverse Heusler alloy (Xa structure)."
            },
            "generate_magnetic_structure": {
                "module": "generators.magnetic.materials",
                "function": "generate_magnetic_structure",
                "params": ['material', 'ordering', 'supercell'],
                "description": "Generate magnetic material with specified ordering."
            },
            "generate_quaternary_heusler": {
                "module": "generators.magnetic.heusler",
                "function": "generate_quaternary_heusler",
                "params": ['formula', 'supercell'],
                "description": "Generate quaternary Heusler alloy (XX'YZ, LiMgPdSn-type)."
            },
        }
    },
    
    # -------------------------------------------------------------------------
    # MOLECULE (Category)
    # -------------------------------------------------------------------------
    "molecule": {
        "description": "Molecule structures",
        "operations": {
            "generate_cage_compound": {
                "module": "generators.molecule.cages",
                "function": "generate_cage_compound",
                "params": ['cage_type'],
                "description": "Wrapper to generate any cage compound."
            },
            "generate_cof": {
                "module": "generators.molecule.frameworks",
                "function": "generate_cof",
                "params": ['cof_name', 'n_layers', 'stacking'],
                "description": "Generate COF structure."
            },
            "generate_corrole": {
                "module": "generators.molecule.porphyrins",
                "function": "generate_corrole",
                "params": ['name', 'meso_substituents'],
                "description": "Generate corrole (contracted porphyrin with 3 meso positions)."
            },
            "generate_cyclopentadienyl_ring": {
                "module": "generators.molecule.organometallics",
                "function": "generate_cyclopentadienyl_ring",
                "params": ['radius', 'z', 'rotation'],
                "description": "Generate Cp ring atom positions."
            },
            "generate_dipeptide_library": {
                "module": "generators.molecule.biomolecules",
                "function": "generate_dipeptide_library",
                "params": ['output_format'],
                "description": "Generate all 400 dipeptide combinations."
            },
            "generate_dna_strand": {
                "module": "generators.molecule.biomolecules",
                "function": "generate_dna_strand",
                "params": ['sequence', 'form', 'single_stranded'],
                "description": "Generate DNA structure."
            },
            "generate_endohedral_fullerene": {
                "module": "generators.molecule.cages",
                "function": "generate_endohedral_fullerene",
                "params": ['system'],
                "description": "Generate endohedral fullerene."
            },
            "generate_fullerene": {
                "module": "generators.molecule.cages",
                "function": "generate_fullerene",
                "params": ['fullerene', 'optimize'],
                "description": "Generate fullerene structure."
            },
            "generate_g_quadruplex": {
                "module": "generators.molecule.biomolecules",
                "function": "generate_g_quadruplex",
                "params": ['sequence', 'topology', 'cation'],
                "description": "Generate G-quadruplex structure."
            },
            "generate_graphene_quantum_dot": {
                "module": "generators.molecule.carbon_nanostructures",
                "function": "generate_graphene_quantum_dot",
                "params": ['shape', 'size_nm', 'edge_type', 'edge_passivation'],
                "description": "Generate graphene quantum dot."
            },
            "generate_metal_cluster": {
                "module": "generators.molecule.cages",
                "function": "generate_metal_cluster",
                "params": ['cluster', 'ligand'],
                "description": "Generate metal cluster."
            },
            "generate_metallocene": {
                "module": "generators.molecule.organometallics",
                "function": "generate_metallocene",
                "params": ['name', 'conformation', 'add_hydrogens'],
                "description": "Generate metallocene structure."
            },
            "generate_mof": {
                "module": "generators.molecule.frameworks",
                "function": "generate_mof",
                "params": ['mof_name', 'supercell', 'defective', 'defect_fraction'],
                "description": "Generate MOF structure."
            },
            "generate_molecule": {
                "module": "generators.molecule.small_molecules",
                "function": "generate_molecule",
                "params": ['formula', 'optimize'],
                "description": "Generate small molecule structure."
            },
            "generate_nanobud": {
                "module": "generators.molecule.carbon_nanostructures",
                "function": "generate_nanobud",
                "params": ['tube_chirality', 'n_buds', 'bud_type', 'attachment_pattern'],
                "description": "Generate nanobud (CNT with fullerene attachments)."
            },
            "generate_nanohorn": {
                "module": "generators.molecule.carbon_nanostructures",
                "function": "generate_nanohorn",
                "params": ['horn_type', 'n_atoms', 'closed_tip'],
                "description": "Generate carbon nanohorn structure."
            },
            "generate_nhc_carbene": {
                "module": "generators.molecule.organometallics",
                "function": "generate_nhc_carbene",
                "params": ['metal', 'substituents', 'oxidation'],
                "description": "Generate N-heterocyclic carbene complex."
            },
            "generate_peapod": {
                "module": "generators.molecule.carbon_nanostructures",
                "function": "generate_peapod",
                "params": ['configuration', 'n_fullerenes', 'spacing_factor'],
                "description": "Generate carbon nanotube peapod (fullerenes inside CNT)."
            },
            "generate_peptide": {
                "module": "generators.molecule.biomolecules",
                "function": "generate_peptide",
                "params": ['sequence', 'secondary_structure', 'termini'],
                "description": "Generate peptide structure."
            },
            "generate_phthalocyanine": {
                "module": "generators.molecule.porphyrins",
                "function": "generate_phthalocyanine",
                "params": ['name', 'peripheral_substituents', 'axial_ligand'],
                "description": "Generate phthalocyanine structure."
            },
            "generate_piano_stool": {
                "module": "generators.molecule.organometallics",
                "function": "generate_piano_stool",
                "params": ['name', 'add_hydrogens'],
                "description": "Generate piano-stool complex."
            },
            "generate_porphyrin": {
                "module": "generators.molecule.porphyrins",
                "function": "generate_porphyrin",
                "params": ['metalloporphyrin', 'meso_substituents', 'axial_ligands'],
                "description": "Generate porphyrin or metalloporphyrin structure."
            },
            "generate_subphthalocyanine": {
                "module": "generators.molecule.porphyrins",
                "function": "generate_subphthalocyanine",
                "params": ['metal', 'axial_ligand'],
                "description": "Generate subphthalocyanine (bowl-shaped with 3 isoindole units)."
            },
            "generate_tripeptide_library": {
                "module": "generators.molecule.biomolecules",
                "function": "generate_tripeptide_library",
                "params": ['output_format'],
                "description": "Generate tripeptide library statistics."
            },
            # New: Conformer and isomer generation
            "generate_conformers": {
                "module": "generators.molecule.conformers",
                "function": "generate_conformers",
                "params": ['smiles', 'n_conformers', 'method', 'prune_threshold', 'optimize', 'force_field'],
                "description": "Generate molecular conformers from SMILES."
            },
            "generate_tautomers": {
                "module": "generators.molecule.conformers",
                "function": "generate_tautomers",
                "params": ['smiles', 'max_tautomers'],
                "description": "Generate tautomers of a molecule."
            },
            "enumerate_stereoisomers": {
                "module": "generators.molecule.conformers",
                "function": "enumerate_stereoisomers",
                "params": ['smiles', 'include_undefined_only'],
                "description": "Enumerate all stereoisomers (R/S, E/Z)."
            },
            "predict_protonation_states": {
                "module": "generators.molecule.conformers",
                "function": "predict_protonation_states",
                "params": ['smiles', 'ph', 'ph_range'],
                "description": "Predict protonation states at given pH."
            },
            "generate_rotamers": {
                "module": "generators.molecule.conformers",
                "function": "generate_rotamers",
                "params": ['smiles', 'rotatable_bond_idx', 'n_steps'],
                "description": "Generate rotamers by rotating around a specific bond."
            },
            "generate_ligand_conformers": {
                "module": "generators.molecule.conformers",
                "function": "generate_ligand_conformers",
                "params": ['smiles', 'n_conformers', 'binding_constraints'],
                "description": "Generate ligand conformers for docking."
            },
        }
    },

    # -------------------------------------------------------------------------
    # NANOTUBE (Category)
    # -------------------------------------------------------------------------
    "nanotube": {
        "description": "Nanotube structures",
        "operations": {
            "generate_axial_heterostructure": {
                "module": "generators.nanotube.nanowires",
                "function": "generate_axial_heterostructure",
                "params": ['materials', 'segment_lengths_nm', 'diameter_nm'],
                "description": "Generate axial heterostructure nanowire."
            },
            "generate_cnt": {
                "module": "generators.nanotube.cnt",
                "function": "generate_cnt",
                "params": ['n', 'm', 'length_nm', 'capped', 'open_ended'],
                "description": "Generate single-walled carbon nanotube."
            },
            "generate_cnt_bundle": {
                "module": "generators.nanotube.cnt",
                "function": "generate_cnt_bundle",
                "params": ['chirality', 'n_tubes', 'length_nm', 'bundle_arrangement'],
                "description": "Generate CNT bundle."
            },
            "generate_coreshell_nanowire": {
                "module": "generators.nanotube.nanowires",
                "function": "generate_coreshell_nanowire",
                "params": ['system', 'core_diameter_nm', 'shell_thickness_nm', 'length_nm'],
                "description": "Generate core-shell nanowire."
            },
            "generate_doped_cnt": {
                "module": "generators.nanotube.cnt",
                "function": "generate_doped_cnt",
                "params": ['n', 'm', 'length_nm', 'dopant', 'doping_level'],
                "description": "Generate doped CNT."
            },
            "generate_functionalized_cnt": {
                "module": "generators.nanotube.cnt",
                "function": "generate_functionalized_cnt",
                "params": ['n', 'm', 'length_nm', 'functional_group', 'coverage'],
                "description": "Generate functionalized CNT."
            },
            "generate_mwcnt": {
                "module": "generators.nanotube.cnt",
                "function": "generate_mwcnt",
                "params": ['inner_chirality', 'n_walls', 'length_nm', 'wall_spacing'],
                "description": "Generate multi-walled carbon nanotube."
            },
            "generate_nanowire": {
                "module": "generators.nanotube.nanowires",
                "function": "generate_nanowire",
                "params": ['material', 'diameter_nm', 'length_nm', 'growth_direction', 'cross_section', 'passivated'],
                "description": "Generate semiconductor/metal nanowire."
            },
        }
    },
    
    # -------------------------------------------------------------------------
    # OUTPUT_FORMATS (Category)
    # -------------------------------------------------------------------------
    "output_formats": {
        "description": "Output Formats structures",
        "operations": {
            "export_cif": {
                "module": "generators.output_formats.converters",
                "function": "export_cif",
                "params": ['structure', 'data_name'],
                "description": "Export structure to CIF format."
            },
            "export_lammps": {
                "module": "generators.output_formats.converters",
                "function": "export_lammps",
                "params": ['structure', 'atom_style', 'units'],
                "description": "Export structure to LAMMPS data file format."
            },
            "export_pdb": {
                "module": "generators.output_formats.converters",
                "function": "export_pdb",
                "params": ['structure', 'title'],
                "description": "Export structure to PDB format."
            },
            "export_quantum_espresso": {
                "module": "generators.output_formats.converters",
                "function": "export_quantum_espresso",
                "params": ['structure', 'pseudopot_dir', 'ecutwfc', 'ecutrho'],
                "description": "Export structure to Quantum ESPRESSO pw.x input format."
            },
            "export_vasp": {
                "module": "generators.output_formats.converters",
                "function": "export_vasp",
                "params": ['structure', 'direct', 'comment'],
                "description": "Export structure to VASP POSCAR format."
            },
            "export_xyz": {
                "module": "generators.output_formats.converters",
                "function": "export_xyz",
                "params": ['structure', 'comment'],
                "description": "Export structure to XYZ format."
            },
            # New: Database adapters
            "export_ase_atoms": {
                "module": "generators.output_formats.database_adapters",
                "function": "export_ase_atoms",
                "params": ['structure', 'calculator', 'add_constraints'],
                "description": "Convert structure to ASE Atoms object."
            },
            "export_aiida_structuredata": {
                "module": "generators.output_formats.database_adapters",
                "function": "export_aiida_structuredata",
                "params": ['structure', 'node_label', 'node_description'],
                "description": "Convert structure to AiiDA StructureData format."
            },
            "export_aflow": {
                "module": "generators.output_formats.database_adapters",
                "function": "export_aflow",
                "params": ['structure', 'auid', 'keywords'],
                "description": "Export structure in AFLOW-compatible format."
            },
            "export_materials_project": {
                "module": "generators.output_formats.database_adapters",
                "function": "export_materials_project",
                "params": ['structure', 'material_id', 'task_id'],
                "description": "Export structure in Materials Project format."
            },
            "export_optimade": {
                "module": "generators.output_formats.database_adapters",
                "function": "export_optimade",
                "params": ['structure', 'entry_id', 'database_id'],
                "description": "Export structure in OPTIMADE JSON format."
            },
            "export_jarvis": {
                "module": "generators.output_formats.database_adapters",
                "function": "export_jarvis",
                "params": ['structure', 'jid', 'dataset'],
                "description": "Export structure in JARVIS-compatible format."
            },
        }
    },

    # -------------------------------------------------------------------------
    # PHOTONIC (Category)
    # -------------------------------------------------------------------------
    "photonic": {
        "description": "Photonic structures",
        "operations": {
            "generate_2d_photonic_crystal": {
                "module": "generators.photonic.crystals",
                "function": "generate_2d_photonic_crystal",
                "params": ['crystal_type', 'material', 'lattice_constant_nm', 'hole_radius_fraction', 'slab_thickness_nm', 'supercell'],
                "description": "Generate 2D photonic crystal slab."
            },
            "generate_3d_photonic_crystal": {
                "module": "generators.photonic.crystals",
                "function": "generate_3d_photonic_crystal",
                "params": ['crystal_type', 'material', 'lattice_constant_um', 'fill_fraction', 'n_periods'],
                "description": "Generate 3D photonic crystal."
            },
            "generate_broken_symmetry_photonic": {
                "module": "generators.photonic.metamaterials",
                "function": "generate_broken_symmetry_photonic",
                "params": ['base_type', 'symmetry_breaking', 'defect_position', 'material'],
                "description": "Generate photonic crystal with broken symmetry (defect modes)."
            },
            "generate_metasurface_unit_cell": {
                "module": "generators.photonic.metamaterials",
                "function": "generate_metasurface_unit_cell",
                "params": ['element_type', 'material', 'period_nm', 'element_size_nm', 'rotation_deg', 'substrate'],
                "description": "Generate metasurface unit cell."
            },
            "generate_photonic_crystal": {
                "module": "generators.photonic.metamaterials",
                "function": "generate_photonic_crystal",
                "params": ['crystal_type', 'material', 'lattice_constant_nm', 'fill_fraction', 'n_periods'],
                "description": "Generate photonic crystal structure."
            },
            "generate_photonic_defect": {
                "module": "generators.photonic.crystals",
                "function": "generate_photonic_defect",
                "params": ['base_crystal', 'defect_type', 'material', 'lattice_constant_nm'],
                "description": "Generate photonic crystal with defect for cavity/waveguide."
            },
        }
    },
    
    # -------------------------------------------------------------------------
    # QUALITY_CONTROL (Category)
    # -------------------------------------------------------------------------
    "quality_control": {
        "description": "Quality Control and Symmetry Analysis",
        "operations": {
            "analyze_symmetry": {
                "module": "generators.quality_control.symmetry",
                "function": "analyze_symmetry",
                "params": ['structure', 'symprec', 'angle_tolerance'],
                "description": "Analyze crystal symmetry - find space group, point group, etc."
            },
            "generate_kpath": {
                "module": "generators.quality_control.kpath",
                "function": "generate_kpath",
                "params": ['lattice_type', 'path', 'n_points'],
                "description": "Generate k-path for band structure."
            },
            "get_subgroups": {
                "module": "generators.quality_control.symmetry",
                "function": "get_subgroups",
                "params": ['space_group', 'strict'],
                "description": "Get maximal subgroups of a space group."
            },
            "get_symmetry_path": {
                "module": "generators.quality_control.symmetry",
                "function": "get_symmetry_path",
                "params": ['start_spg', 'end_spg', 'max_depth'],
                "description": "Find group-subgroup path between two space groups."
            },
            "get_all_subgroup_paths": {
                "module": "generators.quality_control.symmetry",
                "function": "get_all_subgroup_paths",
                "params": ['start_spg', 'end_spg', 'max_depth'],
                "description": "Find ALL possible group-subgroup paths between two space groups."
            },
            "transform_to_subgroup": {
                "module": "generators.quality_control.symmetry",
                "function": "transform_to_subgroup",
                "params": ['structure', 'target_spacegroup', 'eps', 'max_steps'],
                "description": "Transform structure to a subgroup (Bilbao TRANSTRU-like). "
                              "E.g., transform cubic BiFeO3 to rhombohedral for phase comparison."
            },
            "transform_by_path": {
                "module": "generators.quality_control.symmetry",
                "function": "transform_by_path",
                "params": ['structure', 'path', 'eps'],
                "description": "Transform structure following a specific group-subgroup path."
            },
            "validate_structure": {
                "module": "generators.quality_control.symmetry",
                "function": "validate_structure",
                "params": ['structure'],
                "description": "Validate crystal structure (check overlaps, etc)."
            },
            # New: Structure optimization
            "optimize_structure_ase": {
                "module": "generators.quality_control.optimizer",
                "function": "optimize_structure_ase",
                "params": ['structure', 'optimizer', 'fmax', 'steps', 'fix_cell', 'constraints'],
                "description": "Optimize structure using ASE with classical force fields."
            },
            "optimize_structure_mlff": {
                "module": "generators.quality_control.optimizer",
                "function": "optimize_structure_mlff",
                "params": ['structure', 'model', 'optimizer', 'fmax', 'steps', 'fix_cell', 'device'],
                "description": "Optimize structure using ML force field (CHGNet, M3GNet, MACE)."
            },
            "relax_cell_parameters": {
                "module": "generators.quality_control.optimizer",
                "function": "relax_cell_parameters",
                "params": ['structure', 'model', 'target_pressure_GPa', 'mask'],
                "description": "Relax cell parameters under optional pressure constraint."
            },
            "calculate_elastic_tensor": {
                "module": "generators.quality_control.optimizer",
                "function": "calculate_elastic_tensor",
                "params": ['structure', 'model', 'strain_magnitude'],
                "description": "Calculate elastic tensor from stress-strain relationships."
            },
            "generate_piezoelectric_template": {
                "module": "generators.quality_control.optimizer",
                "function": "generate_piezoelectric_template",
                "params": ['material', 'supercell', 'polarization_axis'],
                "description": "Generate template for piezoelectric constant calculation."
            },
        }
    },

    # -------------------------------------------------------------------------
    # QUANTUM (Category)
    # -------------------------------------------------------------------------
    "quantum": {
        "description": "Quantum structures",
        "operations": {
            "generate_coreshell_qd": {
                "module": "generators.quantum.dots",
                "function": "generate_coreshell_qd",
                "params": ['system', 'core_diameter_nm', 'shell_thickness_nm', 'gradient_shell'],
                "description": "Generate core-shell quantum dot."
            },
            "generate_nv_center": {
                "module": "generators.quantum.qubits",
                "function": "generate_nv_center",
                "params": ['supercell', 'n_nv_centers'],
                "description": "Generate diamond with NV center(s)."
            },
            "generate_quantum_dot": {
                "module": "generators.quantum.dots",
                "function": "generate_quantum_dot",
                "params": ['material', 'diameter_nm', 'shape', 'ligand', 'faceted'],
                "description": "Generate quantum dot structure."
            },
            "generate_qubit_structure": {
                "module": "generators.quantum.qubits",
                "function": "generate_qubit_structure",
                "params": ['platform', 'supercell'],
                "description": "Generate qubit material structure."
            },
            "generate_spin_qubit_array": {
                "module": "generators.quantum.qubits",
                "function": "generate_spin_qubit_array",
                "params": ['platform', 'n_qubits', 'spacing_nm'],
                "description": "Generate array of spin qubits."
            },
        }
    },
    
    # -------------------------------------------------------------------------
    # SURFACE (Category)
    # -------------------------------------------------------------------------
    "surface": {
        "description": "Surface structures",
        "operations": {
            "generate_adatom_superlattice": {
                "module": "generators.surface.adatoms",
                "function": "generate_adatom_superlattice",
                "params": ['system', 'pattern', 'n_layers', 'vacuum', 'site_type'],
                "description": "Generate adatom superlattice on metal surface."
            },
            "generate_advacancy_superlattice": {
                "module": "generators.surface.adatoms",
                "function": "generate_advacancy_superlattice",
                "params": ['element', 'facet', 'pattern', 'n_layers', 'vacuum'],
                "description": "Generate ordered vacancy pattern on surface."
            },
            "generate_dealloyed_surface": {
                "module": "generators.surface.alloys",
                "function": "generate_dealloyed_surface",
                "params": ['parent_alloy', 'composition', 'dealloying_depth', 'porosity', 'size', 'n_layers', 'vacuum', 'seed'],
                "description": "Generate de-alloyed porous surface."
            },
            "generate_interface_with_defects": {
                "module": "generators.surface.interfaces",
                "function": "generate_interface_with_defects",
                "params": ['interface', 'n_layers_top', 'n_layers_bottom', 'defect_type', 'defect_concentration', 'defect_location', 'vacuum', 'seed'],
                "description": "Generate interface with point defects."
            },
            "generate_kinked_surface": {
                "module": "generators.surface.slabs",
                "function": "generate_kinked_surface",
                "params": ['base_surface', 'kink_density', 'kink_type'],
                "description": "Generate surface with kink sites."
            },
            "generate_nanoparticle": {
                "module": "generators.surface.nanoparticles",
                "function": "generate_nanoparticle",
                "params": ['element', 'shape', 'size_nm', 'use_magic_number'],
                "description": "Generate nanoparticle structure."
            },
            "generate_oxide_superlattice": {
                "module": "generators.surface.interfaces",
                "function": "generate_oxide_superlattice",
                "params": ['materials', 'layer_thicknesses', 'n_periods', 'vacuum'],
                "description": "Generate oxide superlattice."
            },
            "generate_polar_interface": {
                "module": "generators.surface.interfaces",
                "function": "generate_polar_interface",
                "params": ['interface', 'n_layers_top', 'n_layers_bottom', 'top_termination', 'bottom_termination', 'vacuum', 'compensation'],
                "description": "Generate polar oxide interface."
            },
            "generate_reconstruction": {
                "module": "generators.surface.reconstructions",
                "function": "generate_reconstruction",
                "params": ['surface', 'n_layers', 'vacuum'],
                "description": "Generate reconstructed surface."
            },
            "generate_roughened_surface": {
                "module": "generators.surface.stepped",
                "function": "generate_roughened_surface",
                "params": ['element', 'facet', 'roughness_rms', 'n_layers', 'size', 'vacuum', 'seed'],
                "description": "Generate roughened surface with controlled roughness."
            },
            "generate_segregation_profile": {
                "module": "generators.surface.alloys",
                "function": "generate_segregation_profile",
                "params": ['alloy', 'bulk_composition', 'surface_enrichment', 'profile_type', 'decay_length', 'size', 'n_layers', 'vacuum', 'seed'],
                "description": "Generate alloy with segregation profile."
            },
            "generate_slab": {
                "module": "generators.surface.slabs",
                "function": "generate_slab",
                "params": ['bulk_structure', 'miller_index', 'min_slab_thickness', 'min_vacuum', 'in_unit_planes', 'center_slab', 'symmetric', 'primitive'],
                "description": "Generate slab from bulk structure."
            },
            "generate_stepped_surface_with_kinks": {
                "module": "generators.surface.stepped",
                "function": "generate_stepped_surface_with_kinks",
                "params": ['element', 'base_facet', 'terrace_width', 'kink_spacing', 'n_layers', 'vacuum'],
                "description": "Generate stepped surface with kinks."
            },
            "generate_surface_alloy": {
                "module": "generators.surface.alloys",
                "function": "generate_surface_alloy",
                "params": ['alloy_type', 'size', 'n_layers', 'vacuum'],
                "description": "Generate surface alloy structure."
            },
            "generate_surface_with_adterrace": {
                "module": "generators.surface.slabs",
                "function": "generate_surface_with_adterrace",
                "params": ['base_surface', 'adterrace_size', 'adterrace_height'],
                "description": "Generate surface with ad-terrace (island)."
            },
            "generate_vicinal_surface": {
                "module": "generators.surface.slabs",
                "function": "generate_vicinal_surface",
                "params": ['base_surface', 'miscut_angle_deg', 'miscut_direction', 'n_terraces'],
                "description": "Generate vicinal (stepped) surface."
            },
            "generate_wulff_nanoparticle": {
                "module": "generators.surface.nanoparticles",
                "function": "generate_wulff_nanoparticle",
                "params": ['element', 'size_nm', 'custom_ratios'],
                "description": "Generate Wulff-shape nanoparticle based on surface energies."
            },
        }
    },
    
    # -------------------------------------------------------------------------
    # THERMOELECTRIC (Category)
    # -------------------------------------------------------------------------
    "thermoelectric": {
        "description": "Thermoelectric structures",
        "operations": {
            "generate_thermoelectric": {
                "module": "generators.thermoelectric.materials",
                "function": "generate_thermoelectric",
                "params": ['material', 'supercell'],
                "description": "Generate thermoelectric material structure."
            },
        }
    },
    
    # -------------------------------------------------------------------------
    # TWIST (Category)
    # -------------------------------------------------------------------------
    "twist": {
        "description": "Twist structures",
        "operations": {
            "generate_ferroelectric_2d": {
                "module": "generators.twist.ferroelectric",
                "function": "generate_ferroelectric_2d",
                "params": ['material', 'polarization_state', 'n_layers', 'size', 'vacuum'],
                "description": "Generate 2D ferroelectric structure."
            },
            "generate_heterostructure_stack": {
                "module": "generators.twist.multilayers",
                "function": "generate_heterostructure_stack",
                "params": ['materials', 'twist_angles', 'supercell'],
                "description": "Generate custom heterostructure stack."
            },
            "generate_intercalated_bilayer": {
                "module": "generators.twist.intercalation",
                "function": "generate_intercalated_bilayer",
                "params": ['host_material', 'intercalant', 'stage', 'concentration', 'supercell_size', 'vacuum'],
                "description": "Generate intercalated bilayer structure."
            },
            "generate_magic_angle_graphene": {
                "module": "generators.twist.moire",
                "function": "generate_magic_angle_graphene",
                "params": ['magic_order', 'vacuum'],
                "description": "Generate magic-angle twisted bilayer graphene (TBG)."
            },
            "generate_moire_superlattice": {
                "module": "generators.twist.moire",
                "function": "generate_moire_superlattice",
                "params": ['material', 'twist_angle', 'n_layers', 'vacuum', 'use_commensurate'],
                "description": "Generate Moiré superlattice structure."
            },
            "generate_multilayer_stack": {
                "module": "generators.twist.multilayers",
                "function": "generate_multilayer_stack",
                "params": ['system', 'supercell'],
                "description": "Generate multilayer twisted structure."
            },
            "generate_sliding_ferroelectric": {
                "module": "generators.twist.ferroelectric",
                "function": "generate_sliding_ferroelectric",
                "params": ['material', 'stacking', 'twist_angle', 'vacuum'],
                "description": "Generate sliding ferroelectric bilayer."
            },
            "generate_stacked_layers": {
                "module": "generators.twist.stacking",
                "function": "generate_stacked_layers",
                "params": ['material', 'stacking', 'n_layers', 'vacuum'],
                "description": "Generate bulk stacked layered structure."
            },
            "generate_staged_graphite": {
                "module": "generators.twist.intercalation",
                "function": "generate_staged_graphite",
                "params": ['intercalant', 'stage', 'n_galleries', 'vacuum'],
                "description": "Generate staged intercalated graphite."
            },
            "generate_turbostratic_graphite": {
                "module": "generators.twist.stacking",
                "function": "generate_turbostratic_graphite",
                "params": ['n_layers', 'disorder_type', 'disorder_amount', 'vacuum', 'seed'],
                "description": "Generate turbostratic (disordered stacking) graphite."
            },
            "generate_twisted_bilayer": {
                "module": "generators.twist.bilayers",
                "function": "generate_twisted_bilayer",
                "params": ['material', 'twist_angle', 'stacking', 'supercell_size', 'vacuum'],
                "description": "Generate twisted homo-bilayer structure."
            },
            "generate_twisted_bilayer_simple": {
                "module": "generators.twist.base",
                "function": "generate_twisted_bilayer_simple",
                "params": ['material', 'twist_angle_deg', 'interlayer_distance_A'],
                "description": "Generate simplified twisted bilayer structure (small angle approximation)."
            },
            "generate_twisted_heterobilayer": {
                "module": "generators.twist.bilayers",
                "function": "generate_twisted_heterobilayer",
                "params": ['layer1', 'layer2', 'twist_angle', 'supercell_size', 'strain_distribution', 'vacuum'],
                "description": "Generate twisted hetero-bilayer (different materials)."
            },
            "generate_twisted_multilayer": {
                "module": "generators.twist.multilayers",
                "function": "generate_twisted_multilayer",
                "params": ['system', 'supercell'],
                "description": "Alias for generate_multilayer_stack."
            },
            "generate_twisted_trilayer": {
                "module": "generators.twist.multilayers",
                "function": "generate_twisted_trilayer",
                "params": ['stacking', 'twist_angle', 'supercell'],
                "description": "Generate twisted trilayer graphene (wrapper)."
            },
        }
    },
    
    # -------------------------------------------------------------------------
    # TWO_D (Category)
    # -------------------------------------------------------------------------
    "two_d": {
        "description": "Two D structures",
        "operations": {
            "generate_2d_grain_boundary": {
                "module": "generators.two_d.grain_boundaries",
                "function": "generate_2d_grain_boundary",
                "params": ['base_material', 'gb_type', 'length', 'width', 'vacuum'],
                "description": "Generate 2D grain boundary structure."
            },
            "generate_chiral_nanoribbon": {
                "module": "generators.two_d.ribbons",
                "function": "generate_chiral_nanoribbon",
                "params": ['base_material', 'chiral_angle', 'width', 'length', 'passivation', 'vacuum'],
                "description": "Generate chiral nanoribbon with arbitrary edge angle."
            },
            "generate_encapsulated_structure": {
                "module": "generators.two_d.he_2d",
                "function": "generate_encapsulated_structure",
                "params": ['active_layer', 'encapsulant', 'n_encap_layers'],
                "description": "Generate encapsulated 2D material (e.g., hBN/graphene/hBN)."
            },
            "generate_graphene": {
                "module": "generators.two_d.xenes",
                "function": "generate_graphene",
                "params": ['size', 'vacuum', 'stacking', 'n_layers'],
                "description": "Generate graphene using ASE for maximum accuracy."
            },
            "generate_holey_graphene": {
                "module": "generators.two_d.porous",
                "function": "generate_holey_graphene",
                "params": ['hole_size', 'hole_spacing', 'base_size', 'vacuum'],
                "description": "Generate holey graphene with custom hole pattern."
            },
            "generate_janus_2d": {
                "module": "generators.two_d.janus",
                "function": "generate_janus_2d",
                "params": ['formula', 'metal', 'top_chalcogen', 'bottom_chalcogen', 'a', 'vacuum', 'size'],
                "description": "Generate Janus 2D structure."
            },
            "generate_lateral_heterostructure": {
                "module": "generators.two_d.he_2d",
                "function": "generate_lateral_heterostructure",
                "params": ['material1', 'material2', 'interface_type', 'width_nm'],
                "description": "Generate lateral heterostructure with in-plane junction."
            },
            "generate_magnetic_2d": {
                "module": "generators.two_d.magnetic_2d",
                "function": "generate_magnetic_2d",
                "params": ['material', 'magnetic_state', 'size', 'vacuum'],
                "description": "Generate 2D magnetic material structure with magnetic moments."
            },
            "generate_mxene": {
                "module": "generators.two_d.mxenes",
                "function": "generate_mxene",
                "params": ['formula', 'termination', 'size', 'vacuum'],
                "description": "Generate MXene structure."
            },
            "generate_mxene_multilayer": {
                "module": "generators.two_d.mxenes",
                "function": "generate_mxene_multilayer",
                "params": ['formula', 'termination', 'n_layers', 'interlayer_spacing', 'vacuum'],
                "description": "Generate multi-layer MXene structure."
            },
            "generate_silicene": {
                "module": "generators.two_d.xenes",
                "function": "generate_silicene",
                "params": ['size', 'vacuum', 'passivate'],
                "description": "Generate silicene monolayer (buckled Si honeycomb)."
            },
            "generate_janus_tmd": {
                "module": "generators.two_d.tmds",
                "function": "generate_janus_tmd",
                "params": ['formula', 'phase', 'size', 'vacuum'],
                "description": "Generate Janus TMD with asymmetric chalcogens (e.g., MoSSe)."
            },
            "generate_nanoflake": {
                "module": "generators.two_d.flakes",
                "function": "generate_nanoflake",
                "params": ['base_material', 'shape', 'size', 'passivation', 'vacuum'],
                "description": "Generate 2D nanoflake (finite 2D fragment)."
            },
            "generate_nanomesh": {
                "module": "generators.two_d.flakes",
                "function": "generate_nanomesh",
                "params": ['base_material', 'supercell', 'hole_pattern', 'hole_fraction', 'vacuum'],
                "description": "Generate 2D nanomesh (periodic hole array)."
            },
            "generate_nanoribbon": {
                "module": "generators.two_d.ribbons",
                "function": "generate_nanoribbon",
                "params": ['base_material', 'edge_type', 'width', 'length', 'passivation', 'vacuum'],
                "description": "Generate nanoribbon structure."
            },
            "generate_polycrystalline_2d": {
                "module": "generators.two_d.grain_boundaries",
                "function": "generate_polycrystalline_2d",
                "params": ['base_material', 'n_grains', 'box_size', 'vacuum', 'seed'],
                "description": "Generate polycrystalline 2D material with multiple grains."
            },
            "generate_porous_2d": {
                "module": "generators.two_d.porous",
                "function": "generate_porous_2d",
                "params": ['material', 'size', 'vacuum'],
                "description": "Generate porous 2D material structure."
            },
            "generate_quantum_dot": {
                "module": "generators.two_d.flakes",
                "function": "generate_quantum_dot",
                "params": ['base_material', 'diameter', 'shape', 'passivation', 'vacuum'],
                "description": "Generate 2D quantum dot (nm-scale flake)."
            },
            "generate_tmd": {
                "module": "generators.two_d.tmds",
                "function": "generate_tmd",
                "params": ['formula', 'phase', 'size', 'vacuum'],
                "description": "Generate TMD monolayer structure."
            },
            "generate_tmd_heterostructure": {
                "module": "generators.two_d.tmds",
                "function": "generate_tmd_heterostructure",
                "params": ['tmd_1', 'tmd_2', 'alignment', 'vacuum'],
                "description": "Generate TMD vertical heterostructure."
            },
            "generate_vertical_heterostructure": {
                "module": "generators.two_d.he_2d",
                "function": "generate_vertical_heterostructure",
                "params": ['materials', 'supercell', 'twist_angles'],
                "description": "Generate vertical van der Waals heterostructure."
            },
            "generate_xene": {
                "module": "generators.two_d.xenes",
                "function": "generate_xene",
                "params": ['xene_type', 'size', 'vacuum', 'passivate'],
                "description": "Generate X-ene monolayer structure."
            },
            # New: 2D High-Entropy Materials
            "generate_2d_hea": {
                "module": "generators.two_d.hea_2d",
                "function": "generate_2d_hea",
                "params": ['elements', 'anion', 'structure_type', 'composition', 'size', 'vacuum', 'ordering', 'seed'],
                "description": "Generate 2D high-entropy material (oxide, carbide, nitride, TMD, MXene)."
            },
            "generate_2d_heo": {
                "module": "generators.two_d.hea_2d",
                "function": "generate_2d_heo",
                "params": ['elements', 'composition', 'size', 'vacuum', 'ordering', 'seed'],
                "description": "Generate 2D high-entropy oxide."
            },
            "generate_2d_hec": {
                "module": "generators.two_d.hea_2d",
                "function": "generate_2d_hec",
                "params": ['elements', 'composition', 'size', 'vacuum', 'ordering', 'seed'],
                "description": "Generate 2D high-entropy carbide (MXene-like)."
            },
            "generate_2d_he_tmd": {
                "module": "generators.two_d.hea_2d",
                "function": "generate_2d_he_tmd",
                "params": ['metals', 'chalcogen', 'composition', 'size', 'vacuum', 'ordering', 'seed'],
                "description": "Generate 2D high-entropy TMD."
            },
            "generate_2d_he_mxene": {
                "module": "generators.two_d.hea_2d",
                "function": "generate_2d_he_mxene",
                "params": ['metals', 'composition', 'size', 'vacuum', 'termination', 'ordering', 'seed'],
                "description": "Generate 2D high-entropy MXene."
            },
        }
    },

    # -------------------------------------------------------------------------
    # META_STRUCTURES (Category) - NEW
    # -------------------------------------------------------------------------
    "meta_structures": {
        "description": "Meta-structures: Quasicrystals and Artificial Lattices",
        "operations": {
            "generate_quasicrystal": {
                "module": "generators.meta_structures.quasicrystals",
                "function": "generate_quasicrystal",
                "params": ['quasicrystal_type', 'approximant_order', 'element', 'vacuum'],
                "description": "Generate quasicrystal structure or approximant."
            },
            "generate_penrose_tiling": {
                "module": "generators.meta_structures.quasicrystals",
                "function": "generate_penrose_tiling",
                "params": ['order', 'element', 'bond_length', 'vacuum', 'vertex_decoration'],
                "description": "Generate Penrose tiling (5-fold quasicrystal)."
            },
            "generate_fibonacci_approximant": {
                "module": "generators.meta_structures.quasicrystals",
                "function": "generate_fibonacci_approximant",
                "params": ['order', 'element', 'a_short', 'a_long'],
                "description": "Generate 1D Fibonacci chain approximant."
            },
            "generate_icosahedral_quasicrystal": {
                "module": "generators.meta_structures.quasicrystals",
                "function": "generate_icosahedral_quasicrystal",
                "params": ['approximant_order', 'elements', 'a'],
                "description": "Generate icosahedral quasicrystal approximant."
            },
            "generate_kagome_lattice": {
                "module": "generators.meta_structures.artificial_lattices",
                "function": "generate_kagome_lattice",
                "params": ['element', 'a', 'size', 'vacuum', 'magnetic_order', 'buckling'],
                "description": "Generate Kagome lattice (flat band system)."
            },
            "generate_lieb_lattice": {
                "module": "generators.meta_structures.artificial_lattices",
                "function": "generate_lieb_lattice",
                "params": ['elements', 'a', 'size', 'vacuum', 'dimerization'],
                "description": "Generate Lieb lattice (flat band at E=0)."
            },
            "generate_checkerboard_lattice": {
                "module": "generators.meta_structures.artificial_lattices",
                "function": "generate_checkerboard_lattice",
                "params": ['elements', 'a', 'size', 'vacuum', 'crossing_type'],
                "description": "Generate checkerboard lattice."
            },
            "generate_dice_lattice": {
                "module": "generators.meta_structures.artificial_lattices",
                "function": "generate_dice_lattice",
                "params": ['element', 'a', 'size', 'vacuum', 'alpha'],
                "description": "Generate dice (T3) lattice with pseudospin-1 fermions."
            },
            "generate_honeycomb_variants": {
                "module": "generators.meta_structures.artificial_lattices",
                "function": "generate_honeycomb_variants",
                "params": ['variant', 'element', 'a', 'size', 'vacuum', 'alpha'],
                "description": "Generate honeycomb lattice variants (alpha-T3, Kekulé, Haldane, etc.)."
            },
        }
    },

}


def get_all_operations() -> Dict[str, List[str]]:
    """Get all available operations organized by category."""
    result = {}
    for category, data in GENERATOR_REGISTRY.items():
        result[category] = list(data["operations"].keys())
    return result


def get_operation_info(operation: str) -> Dict[str, Any]:
    """Get information about a specific operation."""
    for category, data in GENERATOR_REGISTRY.items():
        if operation in data["operations"]:
            op_info = data["operations"][operation].copy()
            op_info["category"] = category
            return {"success": True, "operation": operation, "info": op_info}
    
    return {
        "success": False,
        "error": {
            "code": "UNKNOWN_OPERATION",
            "message": f"Operation '{operation}' not found",
            "available_operations": get_all_operations()
        }
    }


def list_category_operations(category: str) -> Dict[str, Any]:
    """List all operations in a category."""
    if category not in GENERATOR_REGISTRY:
        return {
            "success": False,
            "error": {
                "code": "UNKNOWN_CATEGORY",
                "message": f"Category '{category}' not found",
                "available_categories": list(GENERATOR_REGISTRY.keys())
            }
        }
    
    data = GENERATOR_REGISTRY[category]
    return {
        "success": True,
        "category": category,
        "description": data["description"],
        "operations": {
            name: info["description"] 
            for name, info in data["operations"].items()
        }
    }