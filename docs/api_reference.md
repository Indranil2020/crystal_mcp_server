# Crystal MCP Generator API Reference

Complete catalog of all **228** available generator operations, organized by category.

## Categories
- [**ADSORPTION**](#adsorption) (14 ops): Adsorption structures
- [**BATTERY**](#battery) (5 ops): Battery structures
- [**BIOMOLECULAR**](#biomolecular) (5 ops): Biomolecular structures
- [**BULK**](#bulk) (51 ops): Bulk structures
- [**CATALYST**](#catalyst) (8 ops): Catalyst structures
- [**DEFECT**](#defect) (18 ops): Defect structures
- [**ELECTRONIC**](#electronic) (8 ops): Electronic structures
- [**EXTERNAL_FIELDS**](#external_fields) (6 ops): External Fields structures
- [**HIGH_PRESSURE**](#high_pressure) (4 ops): High Pressure structures
- [**MAGNETIC**](#magnetic) (6 ops): Magnetic structures
- [**MOLECULE**](#molecule) (24 ops): Molecule structures
- [**NANOTUBE**](#nanotube) (8 ops): Nanotube structures
- [**OUTPUT_FORMATS**](#output_formats) (6 ops): Output Formats structures
- [**PHOTONIC**](#photonic) (6 ops): Photonic structures
- [**QUALITY_CONTROL**](#quality_control) (1 ops): Quality Control structures
- [**QUANTUM**](#quantum) (5 ops): Quantum structures
- [**SURFACE**](#surface) (17 ops): Surface structures
- [**THERMOELECTRIC**](#thermoelectric) (1 ops): Thermoelectric structures
- [**TWIST**](#twist) (15 ops): Twist structures
- [**TWO_D**](#two_d) (20 ops): Two D structures

---

## adsorption
_Adsorption structures_

| Operation | Description | Parameters |
|-----------|-------------|------------|
| `generate_adsorbate_on_surface` | Place adsorbate molecule on surface. | `surface_structure`, `adsorbate`, `site`, `coverage`, `height_offset`, `rotation_deg` |
| `generate_coadsorption` | Generate co-adsorption configuration. | `surface_structure`, `adsorbates`, `sites`, `positions` |
| `generate_confined_il` | Generate confined ionic liquid. | `ionic_liquid`, `confinement`, `pore_diameter_A` |
| `generate_coverage_pattern` | Generate ordered adsorbate coverage pattern. | `surface`, `adsorbate`, `pattern`, `height` |
| `generate_electrode_il_interface` | Generate electrode-ionic liquid interface. | `electrode`, `ionic_liquid`, `electrode_layers`, `il_thickness_A`, `n_ion_pairs` |
| `generate_fragment_adsorbate` | Add dissociated fragment to surface. | `surface`, `fragment`, `site`, `coverage`, `site_index` |
| `generate_ice_layer` | Generate ice layer on surface. | `surface`, `ice_type`, `n_bilayers`, `proton_ordering` |
| `generate_ionic_liquid_bulk` | Generate bulk ionic liquid structure. | `ionic_liquid`, `box_size_A`, `n_ion_pairs` |
| `generate_mixed_sam` | Generate mixed SAM with two components. | `sam1`, `sam2`, `ratio`, `mixing_pattern`, `size` |
| `generate_sac_graphene` | Generate SAC on graphene-based support. | `metal`, `coordination`, `supercell`, `n_sac_sites` |
| `generate_sac_oxide` | Generate SAC on oxide support. | `metal`, `support`, `supercell`, `anchoring_site` |
| `generate_sac_tmd` | Generate SAC on 2D TMD support. | `metal`, `support`, `defect_type`, `supercell` |
| `generate_sam` | Generate self-assembled monolayer. | `sam_type`, `size`, `defect_density` |
| `generate_water_bilayer` | Generate water bilayer on surface. | `surface`, `orientation`, `coverage`, `dissociation_fraction` |

[Back to Top](#categories)

---

## battery
_Battery structures_

| Operation | Description | Parameters |
|-----------|-------------|------------|
| `generate_cathode` | Generate cathode structure. | `cathode_type`, `lithiation`, `supercell` |
| `generate_cathode_electrolyte_interface` | Generate cathode-electrolyte interface. | `cathode`, `electrolyte`, `interface_orientation` |
| `generate_electrolyte` | Generate solid electrolyte structure. | `electrolyte_type`, `supercell`, `li_site_disorder` |
| `generate_lithiated_series` | Generate series of lithiated structures. | `cathode_type`, `n_points` |
| `generate_surface_coated_cathode` | Generate surface-coated cathode. | `cathode_type`, `coating`, `coating_thickness_nm`, `supercell` |

[Back to Top](#categories)

---

## biomolecular
_Biomolecular structures_

| Operation | Description | Parameters |
|-----------|-------------|------------|
| `generate_dna_surface` | Generate DNA origami on 2D surface. | `surface`, `dna_type`, `surface_size` |
| `generate_lipid_bilayer` | Generate lipid bilayer structure. | `lipid`, `size`, `cholesterol_fraction`, `water_buffer`, `seed` |
| `generate_membrane_with_protein` | Generate membrane with embedded protein. | `lipid`, `membrane_size`, `protein_type`, `seed` |
| `generate_peptide_surface` | Generate peptide-semiconductor interface. | `peptide`, `surface`, `surface_size`, `coverage_nmol_nm2` |
| `generate_protein_surface` | Generate protein-surface interface. | `protein`, `surface`, `surface_size`, `orientation`, `coverage` |

[Back to Top](#categories)

---

## bulk
_Bulk structures_

| Operation | Description | Parameters |
|-----------|-------------|------------|
| `apply_strain` | Apply strain to a structure. | `structure_dict`, `strain`, `strain_type` |
| `apply_strain_tensor` | Apply full 3x3 strain tensor to structure. | `structure_dict`, `strain_tensor` |
| `generate_all_cubic_prototypes` | Generate all common cubic structure types for an element. | `element` |
| `generate_antiphase_boundary` | Generate an antiphase boundary (APB). | `structure_dict`, `apb_plane`, `shift_vector` |
| `generate_antisite` | Generate antisite defect by swapping two atoms. | `structure_dict`, `site_a`, `site_b` |
| `generate_bain_path` | Generate structures along Bain transformation path (BCC ↔ FCC). | `element`, `n_steps`, `a_bcc` |
| `generate_breathing_mode` | Generate breathing mode distortion (e.g., breathing Kagome). | `structure_dict`, `amplitude`, `pattern` |
| `generate_cdw_phase` | Generate charge density wave structure. | `material`, `phase`, `temperature` |
| `generate_clathrate` | Generate clathrate structure. | `clathrate`, `empty_cages`, `supercell` |
| `generate_cuprate` | Generate cuprate superconductor structure. | `cuprate_type`, `n_cuo2_planes`, `oxygen_content`, `supercell` |
| `generate_custom_cuprate` | Generate custom cuprate with specified layer sequence. | `n_cuo2_planes`, `charge_reservoir`, `spacer`, `a`, `c_per_layer` |
| `generate_custom_stacking` | Generate a structure with custom stacking sequence. | `element`, `sequence`, `a`, `c_layer`, `bilayer_element` |
| `generate_dislocation_dipole` | Generate dislocation dipole structure. | `structure_dict`, `dislocation_type`, `separation`, `supercell` |
| `generate_dopant_distribution` | Generate a doped structure with specified dopant distribution. | `structure_dict`, `host_element`, `dopant_element`, `concentration`, `distribution`, `seed` |
| `generate_double_perovskite` | Generate double perovskite A2BB'X6 structure. | `a_site`, `b_site_1`, `b_site_2`, `x_site`, `ordering`, `a` |
| `generate_elastic_dipole` | Generate elastic dipole structure. | `structure_dict`, `dipole_type`, `dipole_tensor` |
| `generate_empty_clathrate` | Generate empty clathrate framework (superconducting Si-46). | `framework`, `clathrate_type` |
| `generate_energetic_cocrystal` | Generate energetic cocrystal structure. | `component_a`, `component_b`, `ratio` |
| `generate_energetic_material` | Generate energetic material crystal structure. | `material`, `polymorph` |
| `generate_epitaxial_strain` | Apply epitaxial strain from substrate. | `structure_dict`, `substrate_a`, `substrate_b`, `relax_c`, `poisson_ratio` |
| `generate_frank_loop` | Generate Frank partial dislocation loop. | `structure_dict`, `loop_type`, `radius`, `habit_plane` |
| `generate_from_spacegroup` | Generate structure in any of 230 space groups using PyXtal. | `spacegroup`, `elements`, `composition`, `a`, `b`, `c`, `alpha`, `beta`, `gamma`, `factor` |
| `generate_from_spacegroup_symbol` | Generate structure from space group symbol (e.g., 'Fm-3m', 'P6_3/mmc'). | `symbol`, `elements`, `composition`, `factor` |
| `generate_graded_supercell` | Generate compositionally graded supercell. | `base_structure`, `element_from`, `element_to`, `gradient_axis`, `supercell` |
| `generate_high_pressure_phase` | Generate a known high-pressure phase. | `phase_name` |
| `generate_interstitial` | Generate interstitial defect structure. | `structure_dict`, `element`, `site_type`, `crystal_type`, `position` |
| `generate_jahn_teller` | Generate Jahn-Teller distorted structure. | `structure_dict`, `jt_metal`, `jt_type`, `elongation` |
| `generate_layered_perovskite` | Generate layered perovskite structure. | `a_site`, `a_prime_site`, `b_site`, `x_site`, `n_layers`, `family` |
| `generate_magnetic_ordering` | Generate magnetic ordering on a structure. | `structure_dict`, `magnetic_elements`, `ordering`, `moments`, `axis`, `q_vector` |
| `generate_nte_material` | Generate negative thermal expansion material. | `material` |
| `generate_pbx_model` | Generate plastic-bonded explosive (PBX) model. | `energetic`, `binder`, `energetic_fraction` |
| `generate_peierls_distortion` | Generate Peierls distorted 1D chain structure. | `structure_dict`, `chain_direction`, `dimerization` |
| `generate_perovskite` | Generate simple perovskite ABX3 structure. | `a_site`, `b_site`, `x_site`, `a`, `phase`, `tilting` |
| `generate_polytype` | Generate a polytype structure with specified stacking sequence. | `element`, `polytype`, `a`, `c_layer`, `bilayer_element`, `vacuum` |
| `generate_porous_energetic` | Generate porous energetic material (e.g., for PBX modeling). | `material`, `void_fraction`, `void_shape`, `void_radius` |
| `generate_prismatic_loop` | Generate prismatic dislocation loop. | `structure_dict`, `loop_type`, `radius`, `orientation` |
| `generate_prototype` | Generate common prototype structure. | `prototype`, `elements`, `a`, `c` |
| `generate_prussian_blue` | Generate Prussian blue analog structure. | `metal_a`, `metal_b`, `a` |
| `generate_random_alloy` | Generate random alloy by random substitution. | `structure_dict`, `substitutions`, `seed` |
| `generate_random_alloy_supercell` | Generate random alloy supercell. | `base_structure`, `substitutions`, `supercell`, `seed` |
| `generate_sft` | Generate stacking fault tetrahedron (SFT). | `structure_dict`, `n_vacancies` |
| `generate_slab_supercell` | Generate slab supercell for surface calculations. | `base_structure`, `miller_index`, `slab_thickness`, `vacuum` |
| `generate_spin_peierls` | Generate spin-Peierls distorted structure. | `structure_dict`, `magnetic_element`, `dimerization` |
| `generate_sqs` | Generate Special Quasirandom Structure (SQS) using simulated annealing. | `structure_dict`, `target_composition`, `n_iterations`, `correlation_cutoff`, `seed` |
| `generate_stacking_fault` | Generate a stacking fault structure. | `structure_dict`, `fault_type`, `fault_plane`, `position` |
| `generate_superhydride` | Generate superhydride structure (MHn). | `metal`, `h_content`, `a` |
| `generate_twin_boundary` | Generate a twin boundary structure. | `structure_dict`, `twin_plane`, `n_layers_below`, `n_layers_above`, `twin_type` |
| `generate_vacancy` | Generate vacancy defect structure. | `structure_dict`, `site_index`, `element`, `n_vacancies`, `vacancy_type` |
| `generate_zeolite` | Generate zeolite structure. | `framework`, `cation`, `Si_Al_ratio`, `supercell` |
| `generate_zeolite_with_guest` | Generate zeolite with guest molecules. | `framework`, `guest_molecule`, `loading` |
| `make_supercell` | Create supercell from structure. | `structure`, `scaling` |

[Back to Top](#categories)

---

## catalyst
_Catalyst structures_

| Operation | Description | Parameters |
|-----------|-------------|------------|
| `generate_cluster` | Generate isolated metal cluster. | `cluster_type`, `return_atoms_only` |
| `generate_hea_surface` | Generate High-Entropy Alloy surface. | `hea_type`, `miller_index`, `supercell`, `termination`, `surface_segregation` |
| `generate_photocatalyst_heterojunction` | Generate photocatalyst heterojunction. | `junction_type`, `materials`, `interface_orientation` |
| `generate_saa_surface` | Generate Single Atom Alloy surface. | `saa_type`, `miller_index`, `supercell`, `n_dopant`, `dopant_layer`, `vacuum` |
| `generate_sac_2d_support` | Generate 2D SAC with specific coordination motif. | `metal`, `coordination_motif`, `support`, `supercell`, `n_sac_sites` |
| `generate_strained_catalyst` | Generate strain-engineered catalyst surface. | `catalyst_type`, `strain_type`, `strain_magnitude`, `supercell` |
| `generate_supported_cluster` | Generate oxide-supported metal cluster. | `cluster`, `support`, `surface`, `supercell`, `adsorption_site` |
| `generate_zeolite_cluster` | Generate zeolite-confined metal cluster. | `metal`, `zeolite`, `cluster_size`, `location` |

[Back to Top](#categories)

---

## defect
_Defect structures_

| Operation | Description | Parameters |
|-----------|-------------|------------|
| `generate_amorphous_interface` | Generate crystalline-amorphous interface. | `crystalline`, `amorphous`, `crystal_layers`, `amorphous_thickness_A` |
| `generate_amorphous_structure` | Generate amorphous structure using random placement with constraints. | `material`, `n_atoms`, `density_gcm3`, `seed` |
| `generate_antisite` | Generate antisite defect (swap two atoms). | `host_structure`, `site1`, `site2` |
| `generate_crack_tip` | Generate structure with crack tip. | `element`, `crack_mode`, `crack_length`, `supercell_size` |
| `generate_dislocation` | Generate structure with dislocation. | `element`, `dislocation_type`, `supercell` |
| `generate_f_center` | Generate F-center (color center) by removing anion. | `host_structure`, `anion_site` |
| `generate_grain_boundary` | Generate grain boundary structure. | `element`, `boundary_type`, `supercell` |
| `generate_hea_sqs` | Generate High-Entropy Alloy SQS structure. | `elements`, `composition`, `supercell_size`, `seed` |
| `generate_interstitial` | Generate interstitial defect. | `host_structure`, `interstitial_species`, `interstitial_type`, `position` |
| `generate_melt_quench` | Generate amorphous structure using simulated melt-quench. | `elements`, `composition`, `n_atoms`, `quench_rate` |
| `generate_precipitate` | Generate structure with precipitate. | `matrix_element`, `precipitate_element`, `precipitate_type`, `precipitate_radius`, `supercell_size` |
| `generate_radiation_cascade` | Generate radiation damage cascade structure. | `element`, `pka_energy_keV`, `supercell_size`, `seed` |
| `generate_sqs` | Generate Special Quasirandom Structure. | `alloy`, `composition`, `supercell_size`, `n_iterations`, `seed` |
| `generate_sro_structure` | Generate structure with specified short-range order. | `alloy`, `sro_parameter`, `supercell_size`, `seed` |
| `generate_stacking_fault` | Generate structure with stacking fault. | `element`, `fault_type`, `fault_position`, `supercell` |
| `generate_substitution` | Generate substitutional defect (dopant). | `host_structure`, `site`, `dopant`, `dopant_type` |
| `generate_vacancy` | Generate vacancy defect. | `host_structure`, `vacancy_site`, `vacancy_type` |
| `generate_void` | Generate structure with void. | `element`, `void_shape`, `void_radius`, `supercell_size`, `void_position` |

[Back to Top](#categories)

---

## electronic
_Electronic structures_

| Operation | Description | Parameters |
|-----------|-------------|------------|
| `generate_2deg_interface` | Generate 2DEG heterostructure interface. | `system`, `spacer_nm`, `doping_nm`, `doping_density_cm2` |
| `generate_cuprate_superconductor` | Generate cuprate superconductor with doping control. | `base`, `doping`, `n_CuO_planes`, `supercell` |
| `generate_majorana_nanowire` | Generate Majorana nanowire structure. | `platform`, `wire_length_nm`, `wire_diameter_nm`, `al_shell_thickness_nm`, `tunnel_barrier` |
| `generate_quantum_well` | Generate quantum well heterostructure. | `well_material`, `barrier_material`, `well_thickness_nm`, `barrier_thickness_nm`, `n_wells`, `alloy_fraction` |
| `generate_superconductor` | Generate superconductor structure. | `material`, `supercell` |
| `generate_superlattice` | Generate superlattice structure. | `material1`, `material2`, `period1_nm`, `period2_nm`, `n_periods` |
| `generate_topological_insulator` | Generate topological insulator structure. | `material`, `thickness_QL`, `supercell` |
| `generate_weyl_semimetal` | Generate Weyl or Dirac semimetal structure. | `material`, `supercell` |

[Back to Top](#categories)

---

## external_fields
_External Fields structures_

| Operation | Description | Parameters |
|-----------|-------------|------------|
| `apply_strain_field` | Apply strain field to existing structure. | `structure_dict`, `strain_profile`, `max_strain`, `direction` |
| `generate_antiskyrmion` | Generate antiskyrmion structure. | `material`, `size_nm`, `supercell_size` |
| `generate_magnetic_bobber` | Generate magnetic bobber (truncated skyrmion tube). | `material`, `depth_nm`, `supercell_size` |
| `generate_rippled_structure` | Generate rippled 2D structure. | `material`, `amplitude_angstrom`, `wavelength_angstrom`, `direction`, `size` |
| `generate_skyrmion` | Generate structure with skyrmion spin texture. | `material`, `skyrmion_radius_nm`, `supercell_size`, `helicity`, `polarity` |
| `generate_strained_structure` | Generate uniformly strained structure. | `element`, `strain_xx`, `strain_yy`, `strain_zz`, `strain_type`, `supercell` |

[Back to Top](#categories)

---

## high_pressure
_High Pressure structures_

| Operation | Description | Parameters |
|-----------|-------------|------------|
| `generate_high_pressure_phase` | Generate high-pressure phase structure. | `phase`, `supercell` |
| `generate_polymeric_nitrogen` | Generate polymeric nitrogen structure. | `phase`, `supercell` |
| `generate_super_ionic_phase` | Generate super-ionic phase structure. | `compound`, `phase`, `temperature_K` |
| `generate_superhydride` | Generate superhydride structure. | `compound`, `pressure_GPa`, `supercell` |

[Back to Top](#categories)

---

## magnetic
_Magnetic structures_

| Operation | Description | Parameters |
|-----------|-------------|------------|
| `generate_antiferromagnet` | Generate antiferromagnetic structure. | `material`, `ordering`, `supercell` |
| `generate_full_heusler` | Generate full Heusler alloy (X2YZ, L21 structure). | `formula`, `supercell` |
| `generate_half_heusler` | Generate half Heusler alloy (XYZ, C1b structure). | `formula`, `supercell` |
| `generate_inverse_heusler` | Generate inverse Heusler alloy (Xa structure). | `formula`, `supercell` |
| `generate_magnetic_structure` | Generate magnetic material with specified ordering. | `material`, `ordering`, `supercell` |
| `generate_quaternary_heusler` | Generate quaternary Heusler alloy (XX'YZ, LiMgPdSn-type). | `formula`, `supercell` |

[Back to Top](#categories)

---

## molecule
_Molecule structures_

| Operation | Description | Parameters |
|-----------|-------------|------------|
| `generate_cage_compound` | Wrapper to generate any cage compound. | `cage_type` |
| `generate_cof` | Generate COF structure. | `cof_name`, `n_layers`, `stacking` |
| `generate_corrole` | Generate corrole (contracted porphyrin with 3 meso positions). | `name`, `meso_substituents` |
| `generate_cyclopentadienyl_ring` | Generate Cp ring atom positions. | `radius`, `z`, `rotation` |
| `generate_dipeptide_library` | Generate all 400 dipeptide combinations. | `output_format` |
| `generate_dna_strand` | Generate DNA structure. | `sequence`, `form`, `single_stranded` |
| `generate_endohedral_fullerene` | Generate endohedral fullerene. | `system` |
| `generate_fullerene` | Generate fullerene structure. | `fullerene`, `optimize` |
| `generate_g_quadruplex` | Generate G-quadruplex structure. | `sequence`, `topology`, `cation` |
| `generate_graphene_quantum_dot` | Generate graphene quantum dot. | `shape`, `size_nm`, `edge_type`, `edge_passivation` |
| `generate_metal_cluster` | Generate metal cluster. | `cluster`, `ligand` |
| `generate_metallocene` | Generate metallocene structure. | `name`, `conformation`, `add_hydrogens` |
| `generate_mof` | Generate MOF structure. | `mof_name`, `supercell`, `defective`, `defect_fraction` |
| `generate_molecule` | Generate small molecule structure. | `formula`, `optimize` |
| `generate_nanobud` | Generate nanobud (CNT with fullerene attachments). | `tube_chirality`, `n_buds`, `bud_type`, `attachment_pattern` |
| `generate_nanohorn` | Generate carbon nanohorn structure. | `horn_type`, `n_atoms`, `closed_tip` |
| `generate_nhc_carbene` | Generate N-heterocyclic carbene complex. | `metal`, `substituents`, `oxidation` |
| `generate_peapod` | Generate carbon nanotube peapod (fullerenes inside CNT). | `configuration`, `n_fullerenes`, `spacing_factor` |
| `generate_peptide` | Generate peptide structure. | `sequence`, `secondary_structure`, `termini` |
| `generate_phthalocyanine` | Generate phthalocyanine structure. | `name`, `peripheral_substituents`, `axial_ligand` |
| `generate_piano_stool` | Generate piano-stool complex. | `name`, `add_hydrogens` |
| `generate_porphyrin` | Generate porphyrin or metalloporphyrin structure. | `metalloporphyrin`, `meso_substituents`, `axial_ligands` |
| `generate_subphthalocyanine` | Generate subphthalocyanine (bowl-shaped with 3 isoindole units). | `metal`, `axial_ligand` |
| `generate_tripeptide_library` | Generate tripeptide library statistics. | `output_format` |

[Back to Top](#categories)

---

## nanotube
_Nanotube structures_

| Operation | Description | Parameters |
|-----------|-------------|------------|
| `generate_axial_heterostructure` | Generate axial heterostructure nanowire. | `materials`, `segment_lengths_nm`, `diameter_nm` |
| `generate_cnt` | Generate single-walled carbon nanotube. | `n`, `m`, `length_nm`, `capped`, `open_ended` |
| `generate_cnt_bundle` | Generate CNT bundle. | `chirality`, `n_tubes`, `length_nm`, `bundle_arrangement` |
| `generate_coreshell_nanowire` | Generate core-shell nanowire. | `system`, `core_diameter_nm`, `shell_thickness_nm`, `length_nm` |
| `generate_doped_cnt` | Generate doped CNT. | `n`, `m`, `length_nm`, `dopant`, `doping_level` |
| `generate_functionalized_cnt` | Generate functionalized CNT. | `n`, `m`, `length_nm`, `functional_group`, `coverage` |
| `generate_mwcnt` | Generate multi-walled carbon nanotube. | `inner_chirality`, `n_walls`, `length_nm`, `wall_spacing` |
| `generate_nanowire` | Generate semiconductor/metal nanowire. | `material`, `diameter_nm`, `length_nm`, `growth_direction`, `cross_section`, `passivated` |

[Back to Top](#categories)

---

## output_formats
_Output Formats structures_

| Operation | Description | Parameters |
|-----------|-------------|------------|
| `export_cif` | Export structure to CIF format. | `structure`, `data_name` |
| `export_lammps` | Export structure to LAMMPS data file format. | `structure`, `atom_style`, `units` |
| `export_pdb` | Export structure to PDB format. | `structure`, `title` |
| `export_quantum_espresso` | Export structure to Quantum ESPRESSO pw.x input format. | `structure`, `pseudopot_dir`, `ecutwfc`, `ecutrho` |
| `export_vasp` | Export structure to VASP POSCAR format. | `structure`, `direct`, `comment` |
| `export_xyz` | Export structure to XYZ format. | `structure`, `comment` |

[Back to Top](#categories)

---

## photonic
_Photonic structures_

| Operation | Description | Parameters |
|-----------|-------------|------------|
| `generate_2d_photonic_crystal` | Generate 2D photonic crystal slab. | `crystal_type`, `material`, `lattice_constant_nm`, `hole_radius_fraction`, `slab_thickness_nm`, `supercell` |
| `generate_3d_photonic_crystal` | Generate 3D photonic crystal. | `crystal_type`, `material`, `lattice_constant_um`, `fill_fraction`, `n_periods` |
| `generate_broken_symmetry_photonic` | Generate photonic crystal with broken symmetry (defect modes). | `base_type`, `symmetry_breaking`, `defect_position`, `material` |
| `generate_metasurface_unit_cell` | Generate metasurface unit cell. | `element_type`, `material`, `period_nm`, `element_size_nm`, `rotation_deg`, `substrate` |
| `generate_photonic_crystal` | Generate photonic crystal structure. | `crystal_type`, `material`, `lattice_constant_nm`, `fill_fraction`, `n_periods` |
| `generate_photonic_defect` | Generate photonic crystal with defect for cavity/waveguide. | `base_crystal`, `defect_type`, `material`, `lattice_constant_nm` |

[Back to Top](#categories)

---

## quality_control
_Quality Control structures_

| Operation | Description | Parameters |
|-----------|-------------|------------|
| `generate_kpath` | Generate k-path for band structure. | `lattice_type`, `path`, `n_points` |

[Back to Top](#categories)

---

## quantum
_Quantum structures_

| Operation | Description | Parameters |
|-----------|-------------|------------|
| `generate_coreshell_qd` | Generate core-shell quantum dot. | `system`, `core_diameter_nm`, `shell_thickness_nm`, `gradient_shell` |
| `generate_nv_center` | Generate diamond with NV center(s). | `supercell`, `n_nv_centers` |
| `generate_quantum_dot` | Generate quantum dot structure. | `material`, `diameter_nm`, `shape`, `ligand`, `faceted` |
| `generate_qubit_structure` | Generate qubit material structure. | `platform`, `supercell` |
| `generate_spin_qubit_array` | Generate array of spin qubits. | `platform`, `n_qubits`, `spacing_nm` |

[Back to Top](#categories)

---

## surface
_Surface structures_

| Operation | Description | Parameters |
|-----------|-------------|------------|
| `generate_adatom_superlattice` | Generate adatom superlattice on metal surface. | `system`, `pattern`, `n_layers`, `vacuum`, `site_type` |
| `generate_advacancy_superlattice` | Generate ordered vacancy pattern on surface. | `element`, `facet`, `pattern`, `n_layers`, `vacuum` |
| `generate_dealloyed_surface` | Generate de-alloyed porous surface. | `parent_alloy`, `composition`, `dealloying_depth`, `porosity`, `size`, `n_layers`, `vacuum`, `seed` |
| `generate_interface_with_defects` | Generate interface with point defects. | `interface`, `n_layers_top`, `n_layers_bottom`, `defect_type`, `defect_concentration`, `defect_location`, `vacuum`, `seed` |
| `generate_kinked_surface` | Generate surface with kink sites. | `base_surface`, `kink_density`, `kink_type` |
| `generate_nanoparticle` | Generate nanoparticle structure. | `element`, `shape`, `size_nm`, `use_magic_number` |
| `generate_oxide_superlattice` | Generate oxide superlattice. | `materials`, `layer_thicknesses`, `n_periods`, `vacuum` |
| `generate_polar_interface` | Generate polar oxide interface. | `interface`, `n_layers_top`, `n_layers_bottom`, `top_termination`, `bottom_termination`, `vacuum`, `compensation` |
| `generate_reconstruction` | Generate reconstructed surface. | `surface`, `n_layers`, `vacuum` |
| `generate_roughened_surface` | Generate roughened surface with controlled roughness. | `element`, `facet`, `roughness_rms`, `n_layers`, `size`, `vacuum`, `seed` |
| `generate_segregation_profile` | Generate alloy with segregation profile. | `alloy`, `bulk_composition`, `surface_enrichment`, `profile_type`, `decay_length`, `size`, `n_layers`, `vacuum`, `seed` |
| `generate_slab` | Generate slab from bulk structure. | `bulk_structure`, `miller_index`, `min_slab_thickness`, `min_vacuum`, `in_unit_planes`, `center_slab`, `symmetric`, `primitive` |
| `generate_stepped_surface_with_kinks` | Generate stepped surface with kinks. | `element`, `base_facet`, `terrace_width`, `kink_spacing`, `n_layers`, `vacuum` |
| `generate_surface_alloy` | Generate surface alloy structure. | `alloy_type`, `size`, `n_layers`, `vacuum` |
| `generate_surface_with_adterrace` | Generate surface with ad-terrace (island). | `base_surface`, `adterrace_size`, `adterrace_height` |
| `generate_vicinal_surface` | Generate vicinal (stepped) surface. | `base_surface`, `miscut_angle_deg`, `miscut_direction`, `n_terraces` |
| `generate_wulff_nanoparticle` | Generate Wulff-shape nanoparticle based on surface energies. | `element`, `size_nm`, `custom_ratios` |

[Back to Top](#categories)

---

## thermoelectric
_Thermoelectric structures_

| Operation | Description | Parameters |
|-----------|-------------|------------|
| `generate_thermoelectric` | Generate thermoelectric material structure. | `material`, `supercell` |

[Back to Top](#categories)

---

## twist
_Twist structures_

| Operation | Description | Parameters |
|-----------|-------------|------------|
| `generate_ferroelectric_2d` | Generate 2D ferroelectric structure. | `material`, `polarization_state`, `n_layers`, `size`, `vacuum` |
| `generate_heterostructure_stack` | Generate custom heterostructure stack. | `materials`, `twist_angles`, `supercell` |
| `generate_intercalated_bilayer` | Generate intercalated bilayer structure. | `host_material`, `intercalant`, `stage`, `concentration`, `supercell_size`, `vacuum` |
| `generate_magic_angle_graphene` | Generate magic-angle twisted bilayer graphene (TBG). | `magic_order`, `vacuum` |
| `generate_moire_superlattice` | Generate Moiré superlattice structure. | `material`, `twist_angle`, `n_layers`, `vacuum`, `use_commensurate` |
| `generate_multilayer_stack` | Generate multilayer twisted structure. | `system`, `supercell` |
| `generate_sliding_ferroelectric` | Generate sliding ferroelectric bilayer. | `material`, `stacking`, `twist_angle`, `vacuum` |
| `generate_stacked_layers` | Generate bulk stacked layered structure. | `material`, `stacking`, `n_layers`, `vacuum` |
| `generate_staged_graphite` | Generate staged intercalated graphite. | `intercalant`, `stage`, `n_galleries`, `vacuum` |
| `generate_turbostratic_graphite` | Generate turbostratic (disordered stacking) graphite. | `n_layers`, `disorder_type`, `disorder_amount`, `vacuum`, `seed` |
| `generate_twisted_bilayer` | Generate twisted homo-bilayer structure. | `material`, `twist_angle`, `stacking`, `supercell_size`, `vacuum` |
| `generate_twisted_bilayer_simple` | Generate simplified twisted bilayer structure (small angle approximation). | `material`, `twist_angle_deg`, `interlayer_distance_A` |
| `generate_twisted_heterobilayer` | Generate twisted hetero-bilayer (different materials). | `layer1`, `layer2`, `twist_angle`, `supercell_size`, `strain_distribution`, `vacuum` |
| `generate_twisted_multilayer` | Alias for generate_multilayer_stack. | `system`, `supercell` |
| `generate_twisted_trilayer` | Generate twisted trilayer graphene (wrapper). | `stacking`, `twist_angle`, `supercell` |

[Back to Top](#categories)

---

## two_d
_Two D structures_

| Operation | Description | Parameters |
|-----------|-------------|------------|
| `generate_2d_grain_boundary` | Generate 2D grain boundary structure. | `base_material`, `gb_type`, `length`, `width`, `vacuum` |
| `generate_chiral_nanoribbon` | Generate chiral nanoribbon with arbitrary edge angle. | `base_material`, `chiral_angle`, `width`, `length`, `passivation`, `vacuum` |
| `generate_encapsulated_structure` | Generate encapsulated 2D material (e.g., hBN/graphene/hBN). | `active_layer`, `encapsulant`, `n_encap_layers` |
| `generate_graphene` | Generate graphene using ASE for maximum accuracy. | `size`, `vacuum`, `stacking`, `n_layers` |
| `generate_holey_graphene` | Generate holey graphene with custom hole pattern. | `hole_size`, `hole_spacing`, `base_size`, `vacuum` |
| `generate_janus_2d` | Generate Janus 2D structure. | `formula`, `metal`, `top_chalcogen`, `bottom_chalcogen`, `a`, `vacuum`, `size` |
| `generate_lateral_heterostructure` | Generate lateral heterostructure with in-plane junction. | `material1`, `material2`, `interface_type`, `width_nm` |
| `generate_magnetic_2d` | Generate 2D magnetic material structure with magnetic moments. | `material`, `magnetic_state`, `size`, `vacuum` |
| `generate_mxene` | Generate MXene structure. | `formula`, `termination`, `size`, `vacuum` |
| `generate_mxene_multilayer` | Generate multi-layer MXene structure. | `formula`, `termination`, `n_layers`, `interlayer_spacing`, `vacuum` |
| `generate_nanoflake` | Generate 2D nanoflake (finite 2D fragment). | `base_material`, `shape`, `size`, `passivation`, `vacuum` |
| `generate_nanomesh` | Generate 2D nanomesh (periodic hole array). | `base_material`, `supercell`, `hole_pattern`, `hole_fraction`, `vacuum` |
| `generate_nanoribbon` | Generate nanoribbon structure. | `base_material`, `edge_type`, `width`, `length`, `passivation`, `vacuum` |
| `generate_polycrystalline_2d` | Generate polycrystalline 2D material with multiple grains. | `base_material`, `n_grains`, `box_size`, `vacuum`, `seed` |
| `generate_porous_2d` | Generate porous 2D material structure. | `material`, `size`, `vacuum` |
| `generate_quantum_dot` | Generate 2D quantum dot (nm-scale flake). | `base_material`, `diameter`, `shape`, `passivation`, `vacuum` |
| `generate_tmd` | Generate TMD monolayer structure. | `formula`, `phase`, `size`, `vacuum` |
| `generate_tmd_heterostructure` | Generate TMD vertical heterostructure. | `tmd_1`, `tmd_2`, `alignment`, `vacuum` |
| `generate_vertical_heterostructure` | Generate vertical van der Waals heterostructure. | `materials`, `supercell`, `twist_angles` |
| `generate_xene` | Generate X-ene monolayer structure. | `xene_type`, `size`, `vacuum`, `passivate` |

[Back to Top](#categories)

---
