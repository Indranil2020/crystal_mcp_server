# Crystal MCP Server - Complete Capability Audit Report

**Date:** December 22, 2024
**Objective:** Rigorous verification of listing items against code capabilities.
**Legend:**
- **VERIFIED**: Code exists, is tested, and capability confirmed.
- **IMPLEMENTED**: Code exists but relies on simplified logic or missing dependencies.
- **PARTIAL**: Code exists but misses specific variations listed in Roman numerals.
- **MISSING**: No code found for the specific requirement.

---

### 1. Crystalline Bulk – Periodic 3-D
| ID | Catalogue Item | Status | Source (File:Line) | Capabilities & Gaps |
|:---|:---|:---|:---|:---|
| (i) | 230 Space Groups | **PARTIAL** | `bulk/spacegroups.py:97` | **Supported:** Generates all 230 space groups via PyXtal.<br>**Missing:** Explicit control for Primitive / Conventional / Niggli settings. |
| (ii) | Supercells | **VERIFIED** | `bulk/supercells.py:57` | **Supported:** Integer & rational matrix scaling supported.<br>**Missing:** Advanced rational matrix generation utilities. |
| (iii) | Non-diagonal | **VERIFIED** | `bulk/supercells.py:57` | **Supported:** [make_supercell](file:///home/niel/git/crystal-mcp-server/tests/testsuit/test_mcp_comprehensive.py#557-598) accepts arbitrary 3x3 matrices.<br>**Missing:** Automated generation of triclinic/sheared supercells. |
| (iv) | Long-period | **IMPLEMENTED**| `bulk/supercells.py:169`| **Supported:** `generate_graded_supercell` for long-period structures.<br>**Missing:** Specific commensurate/incommensurate logic. |
| (v) | Boundaries | **IMPLEMENTED**| `bulk/boundaries.py:43` | **Supported:** Twin boundaries, Antiphase boundaries.<br>**Missing:** Crystallographic shear planes. |
| (vi) | Polytypes | **IMPLEMENTED**| `bulk/polytypes.py:49` | **Supported:** Explicit stacking sequences (ABC, ABAC).<br>**Missing:** Automated 6H/15R notation parsing. |
| (vii) | Bain paths | **IMPLEMENTED**| `bulk/strain.py:196` | **Supported:** FCC-BCC Bain path generation.<br>**Missing:** Full strain tensor sampling. |
| (viii)| HEAs | **IMPLEMENTED**| `bulk/disorder.py:23` | **Supported:** Random substitution (equi/non-equi), SQS generation.<br>**Missing:** Short-range-order (SRO) metric optimization. |
| (ix) | Magnetic | **PARTIAL** | [magnetic/materials.py](file:///home/niel/git/crystal-mcp-server/src/python/generators/magnetic/materials.py) | **Supported:** FM, AFM (A/C/G types), FiM, Helical, Skyrmions.<br>**Missing:** 3-q, 4-q, DMI cycloids explicit generation. |
| (x) | Spin-orbit | **MISSING** | - | **Missing:** Magnetic space groups / Shubnikov groups support. |
| (xi) | Defects | **PARTIAL** | [defect/point_defects.py](file:///home/niel/git/crystal-mcp-server/src/python/generators/defect/point_defects.py)| **Supported:** Vacancy (Mono, Di, Schottky, Frenkel), Interstitial (Tet, Oct, Split, Crowdion), Antisite.<br>**Missing:** Trivacancy, Gradient defect profiles. |
| (xii) | Dopants | **IMPLEMENTED**| `bulk/defects.py:232` | **Supported:** Dopant distributions (Gaussian/Error-function).<br>**Missing:** Poisson profiles. |
| (xiii)| Elastic loops | **IMPLEMENTED**| `bulk/dislocations.py` | **Supported:** Dislocation dipoles.<br>**Missing:** Frank loops, Stacking-fault tetrahedra. |
| (xiv) | Quantum phases | **IMPLEMENTED**| `bulk/quantum_phases.py`| **Supported:** CDW, Peierls, Breathing Kagome, Jahn-Teller.<br>**Missing:** Spin-Peierls. |
| (xv) | Cuprates | **IMPLEMENTED**| `bulk/cuprates.py:86` | **Supported:** General cuprate generation logic.<br>**Missing:** Specific infinite-layer / 123 / 214 presets. |
| (xvi) | Clathrates | **IMPLEMENTED**| [bulk/clathrates.py](file:///home/niel/git/crystal-mcp-server/src/python/generators/bulk/clathrates.py) | **Supported:** Clathrate I/II generation with guests.<br>**Missing:** Types III / H / HS-I. |
| (xvii)| Zeolites | **IMPLEMENTED**| `bulk/zeolites.py:94` | **Supported:** Zeolite generation from IZA codes.<br>**Missing:** AlPO frameworks. |
| (xviii)| Prussian-blue | **IMPLEMENTED**| `bulk/special_phases.py`| **Supported:** Prussian Blue Analogues, Double Perovskites.<br>**Missing:** Cyanide perovskites. |
| (xix) | Ferroelectrics | **PARTIAL** | [external_fields/electric_fields.py](file:///home/niel/git/crystal-mcp-server/src/python/generators/external_fields/electric_fields.py)| **Supported:** Ferroelectric domains (180/90), Polarized cells.<br>**Missing:** Hybrid improper ferroelectrics (RP, DJ, Aurivillius). |
| (xxi) | NTE Materials | **IMPLEMENTED**| `bulk/special_phases.py`| **Supported:** Negative Thermal Expansion material generation.<br>**Missing:** Zero thermal expansion engineered cells. |
| (xxii)| Super-hydrides | **IMPLEMENTED**| `bulk/special_phases.py`| **Supported:** LaH10, YH9 generation.<br>**Missing:** - |
| (xxiii)| Energetic | **IMPLEMENTED**| `bulk/energetics.py` | **Supported:** Energetic material crystal generation.<br>**Missing:** Plastic-bonded explosives with voids. |

### 2. 2-D Monolayers & Slabs
| ID | Catalogue Item | Status | Source (File:Line) | Capabilities & Gaps |
|:---|:---|:---|:---|:---|
| (i) | Layer groups | **VERIFIED** | `surface/slabs.py:58` | **Supported:** General slab generation from any bulk.<br>**Missing:** Explicit layer group (80 types) symmetry analysis. |
| (ii) | Magnetic 2D | **IMPLEMENTED**| `two_d/magnetic_2d.py` | **Supported:** CrI3, Fe3GeTe2, CrGeTe3 presets.<br>**Missing:** 2-D Kagome magnets. |
| (iii) | Janus 2D | **IMPLEMENTED**| `two_d/tmds.py:272` | **Supported:** MoSSe and generalized Janus generation (dipole support).<br>**Missing:** - |
| (iv) | X-enes | **VERIFIED** | `two_d/xenes.py:109` | **Supported:** Silicene, Germanene, Stanene, Plumbene, Borophene (beta12), Phosphorene, Arsenene, Antimonene, Bismuthene.<br>**Missing:** - |
| (v) | 2D Oxides | **PARTIAL** | `two_d/mxenes.py:59` | **Supported:** MXenes (Ti3C2, Mo2TiC2), HEOs via [hea_2d.py](file:///home/niel/git/crystal-mcp-server/src/python/generators/two_d/hea_2d.py).<br>**Missing:** Hydroxides (NiFe-LDH) explicit presets. |
| (vi) | 2D HEAs | **IMPLEMENTED**| [two_d/hea_2d.py](file:///home/niel/git/crystal-mcp-server/src/python/generators/two_d/hea_2d.py) | **Supported:** 2-D High-entropy oxides, carbides, nitrides, TMDs, MXenes.<br>**Missing:** - |
| (vii) | Phases | **PARTIAL** | `two_d/tmds.py:58` | **Supported:** 1H, 2H, 1T, 1T', Td, 3R phases.<br>**Missing:** 1T''' phase. |
| (viii)| GBs | **IMPLEMENTED**| `two_d/grain_bound...` | **Supported:** 2D Grain boundaries.<br>**Missing:** Explicit Sigma (Σ3, Σ5, etc.) generator. |
| (ix) | Nanoribbons | **VERIFIED** | `two_d/ribbons.py:75` | **Supported:** Armchair & Zigzag ribbons.<br>**Missing:** Chiral/Bearded edges, specific edge passivation (NO2). |
| (x) | Nanomeshes | **IMPLEMENTED**| `two_d/flakes.py:214` | **Supported:** Nanomesh & Nanoflake generation.<br>**Missing:** Exact edge termination control for QDs. |
| (xi) | Porous 2D | **IMPLEMENTED**| `two_d/porous.py:79` | **Supported:** g-C3N4, C2N, Holey graphene, Graphyne.<br>**Missing:** Biphenylene (Implemented in porous.py). |

### 3. Multilayer Stacks & Twistronics
| ID | Catalogue Item | Status | Source (File:Line) | Capabilities & Gaps |
|:---|:---|:---|:---|:---|
| (i) | Twisted | **VERIFIED** | `twist/bilayers.py:34` | **Supported:** Twisted Homobilayers, Commensurate cell finding.<br>**Missing:** Incommensurate / Lattice-matched modes. |
| (ii) | Multilayers | **IMPLEMENTED**| [twist/multilayers.py](file:///home/niel/git/crystal-mcp-server/src/python/generators/twist/multilayers.py) | **Supported:** Arbitrary layer stacking (Trilayer, Quadlayer).<br>**Missing:** - |
| (iii) | Heterostructures| **VERIFIED** | `twist/bilayers.py:152`| **Supported:** Graphene/hBN, TMD/TMD, etc. twisted heterostructures.<br>**Missing:** - |
| (iv) | Bulk Stacking | **IMPLEMENTED**| `twist/stacking.py:71` | **Supported:** ABA, ABC, AA stacking generation.<br>**Missing:** Rhombohedral / Staggered explicit naming. |
| (v) | Ferroelectrics | **IMPLEMENTED**| `twist/ferroelectric.py`| **Supported:** Sliding ferroelectricity (In2Se3).<br>**Missing:** CuInP2S6 specific support. |
| (vi) | Intercalation | **IMPLEMENTED**| `twist/intercalation.py`| **Supported:** Intercalated layers (Li, K, etc.).<br>**Missing:** Organic molecule intercalation. |

### 4. Surfaces, Slabs, Polar
| ID | Catalogue Item | Status | Source (File:Line) | Capabilities & Gaps |
|:---|:---|:---|:---|:---|
| (i) | Miller Surfaces | **VERIFIED** | `surface/slabs.py:58` | **Supported:** Generic hkl surface generation.<br>**Missing:** Automatic dipole correction logic. |
| (ii) | Reconstructions | **IMPLEMENTED**| [surface/reconstructions.py](file:///home/niel/git/crystal-mcp-server/src/python/generators/surface/reconstructions.py)| **Supported:** Standard reconstructions (2x1, 7x7 etc).<br>**Missing:** Fe3O4(001)-√2×√2R45°. |
| (iii) | Adatom Lattice | **VERIFIED** | `surface/adatoms.py:99` | **Supported:** Ordered adatom/advacancy superlattices.<br>**Missing:** - |
| (iv) | Wulff Shapes | **IMPLEMENTED**| [surface/nanoparticles.py](file:///home/niel/git/crystal-mcp-server/src/python/generators/surface/nanoparticles.py)| **Supported:** Wulff construction based on surface energies.<br>**Missing:** Kinetics-controlled facets. |
| (v) | Core-Shell | **PARTIAL** | `quantum/dots.py:347` | **Supported:** Spherical Core-shell particles.<br>**Missing:** Cubic, Rod, Wire, Octahedron shell shapes. |
| (vi) | Alloying | **IMPLEMENTED**| `surface/alloys.py:63` | **Supported:** Pt-skin / Surface segregation profiles.<br>**Missing:** - |
| (vii) | Roughened | **IMPLEMENTED**| `surface/stepped.py:214` | **Supported:** Kinks, steps, vicinal surfaces.<br>**Missing:** - |
| (viii)| Polar | **IMPLEMENTED**| [external_fields/electric_fields.py](file:///home/niel/git/crystal-mcp-server/src/python/generators/external_fields/electric_fields.py)| **Supported:** Polar surface generation with compensation, Electric fields.<br>**Missing:** Specific LaAlO3/SrTiO3 preset (generic supported). |

### 5. Molecules
| ID | Catalogue Item | Status | Source (File:Line) | Capabilities & Gaps |
|:---|:---|:---|:---|:---|
| (i) | Small Molecules | **IMPLEMENTED**| `molecule/small_mole...`| **Supported:** Basic molecule generation.<br>**Missing:** GDB-17 enumeration / Database. |
| (ii) | Conformers | **IMPLEMENTED** | [molecule/conformers.py](file:///home/niel/git/crystal-mcp-server/src/python/generators/molecule/conformers.py) | **Supported:** RDKit-based generation.<br>**Note:** **UNTESTED** - RDKit not in env. |
| (iii) | Isomers | **IMPLEMENTED** | [molecule/conformers.py](file:///home/niel/git/crystal-mcp-server/src/python/generators/molecule/conformers.py) | **Supported:** Stereoisomers, Tautomers, Rotamers.<br>**Missing:** Atropisomers explicit check. |
| (iv) | States | **IMPLEMENTED** | [molecule/conformers.py](file:///home/niel/git/crystal-mcp-server/src/python/generators/molecule/conformers.py) | **Supported:** Protonation state prediction.<br>**Missing:** Redox state enumerator. |
| (v) | Binding Modes | **IMPLEMENTED** | [molecule/conformers.py](file:///home/niel/git/crystal-mcp-server/src/python/generators/molecule/conformers.py)| **Supported:** Ligand conformer generation for docking.<br>**Missing:** Explicit metal-binder denticity control. |
| (vi) | Organometallics | **IMPLEMENTED**| `molecule/organomet...` | **Supported:** Metallocenes, Piano-stool, Carbene complexes.<br>**Missing:** Grubbs catalysts. |
| (vii) | COFs | **IMPLEMENTED**| [molecule/frameworks.py](file:///home/niel/git/crystal-mcp-server/src/python/generators/molecule/frameworks.py)| **Supported:** 2D/3D COF generation.<br>**Missing:** Specific topology selection (hcb, sql etc). |
| (viii)| MOFs | **IMPLEMENTED**| [molecule/frameworks.py](file:///home/niel/git/crystal-mcp-server/src/python/generators/molecule/frameworks.py)| **Supported:** Common MOF generation.<br>**Missing:** Large topology library (12k+ codes). |
| (ix) | Cages | **IMPLEMENTED**| `molecule/cages.py:325` | **Supported:** M12L24, M24L48 cages.<br>**Missing:** Cucurbiturils, Cryptands. |
| (x) | Fullerenes | **IMPLEMENTED**| `molecule/carbon_nano...`| **Supported:** C60, C70, M@C60.<br>**Missing:** Heterofullerenes (C59N). |
| (xi) | Nanohorns | **IMPLEMENTED**| `molecule/carbon_nano...`| **Supported:** Carbon nanohorns, Peapods, Nanobuds.<br>**Missing:** Nanofunnels. |
| (xii) | Porphyrins | **IMPLEMENTED**| `molecule/porphyrins.py`| **Supported:** Porphyrins, Phthalocyanines, Corroles.<br>**Missing:** Sub-phthalocyanines with axial ligands. |
| (xiii)| Nucleic Acids | **IMPLEMENTED**| [molecule/biomolecules.py](file:///home/niel/git/crystal-mcp-server/src/python/generators/molecule/biomolecules.py)| **Supported:** DNA/RNA single strands.<br>**Missing:** G-quadruplexes (FOUND), i-motifs, PNA. |
| (xiv) | Peptides | **IMPLEMENTED**| [molecule/biomolecules.py](file:///home/niel/git/crystal-mcp-server/src/python/generators/molecule/biomolecules.py)| **Supported:** Alpha-helix, Beta-sheet.<br>**Missing:** 310-helix, PPII (FOUND), Random coil. |
| (xv) | Peptide Libs | **IMPLEMENTED** | [molecule/biomolecules.py](file:///home/niel/git/crystal-mcp-server/src/python/generators/molecule/biomolecules.py)| **Supported:** Dipeptide generator, Tripeptide summary/samples.<br>**Missing:** Full tripeptide library object. |

### 6. Adsorption
| ID | Catalogue Item | Status | Source (File:Line) | Capabilities & Gaps |
|:---|:---|:---|:---|:---|
| (i) | SACs | **IMPLEMENTED**| `adsorption/sac.py:77` | **Supported:** SACs on Graphene, Oxides, TMDs.<br>**Missing:** C2N, g-C3N4, Mo2TiC2Tx specific hosts. |
| (ii) | Ad-molecules | **VERIFIED** | [adsorption/adsorbates.py](file:///home/niel/git/crystal-mcp-server/src/python/generators/adsorption/adsorbates.py)| **Supported:** Standard adsorbates (CO, H2O etc) on surfaces.<br>**Missing:** Complex molecules (Guaiacol, Aniline). |
| (iii) | Fragments | **IMPLEMENTED**| `adsorption/fragments.py`| **Supported:** Dissociated fragments (OH, H, O).<br>**Missing:** Complex fragments (CxHyOz, NHx). |
| (iv) | Patterns | **IMPLEMENTED**| `adsorption/fragments.py`| **Supported:** Coverage patterns.<br>**Missing:** Specific p(2x2), c(4x2) notation parsers. |
| (v) | Co-adsorption | **IMPLEMENTED**| [adsorption/adsorbates.py](file:///home/niel/git/crystal-mcp-server/src/python/generators/adsorption/adsorbates.py)| **Supported:** Multiple species co-adsorption.<br>**Missing:** Lateral interaction modeling. |
| (vi) | Water/Ice | **IMPLEMENTED**| `adsorption/water_layers.py`| **Supported:** Bilayer ice, water adlayers.<br>**Missing:** H-up/H-down explicit control. |
| (vii) | Ionic Liquids | **IMPLEMENTED**| `adsorption/ionic_li...`| **Supported:** IL Monolayers on Graphene.<br>**Missing:** Specific [BMIM][PF6] molecule models. |
| (viii)| SAMs | **IMPLEMENTED**| `adsorption/sams.py:68` | **Supported:** Alkanethiols on Au.<br>**Missing:** Silanes/Phosphonates on Oxides. |
| (ix) | Polymer Graft | **MISSING** | - | **Missing:** Polymer surface grafting (PMMA, PEG). |

### 7. Defects & Disorder
| ID | Catalogue Item | Status | Source (File:Line) | Capabilities & Gaps |
|:---|:---|:---|:---|:---|
| (i) | Point Defects | **PARTIAL** | [defect/point_defects.py](file:///home/niel/git/crystal-mcp-server/src/python/generators/defect/point_defects.py)| **Supported:** Vacancy (Mono, Di, Schottky, Frenkel), Interstitial, Antisite.<br>**Missing:** Trivacancy explicit support. |
| (ii) | Line Defects | **IMPLEMENTED**| `defect/extended_...` | **Supported:** Dislocation loops, dipoles.<br>**Missing:** Explicit Burgers vector library. |
| (iii) | Planar Defects | **IMPLEMENTED**| `defect/extended_...` | **Supported:** Grain boundaries, Stacking faults, Twin boundaries, APBs.<br>**Missing:** Inversion-domain boundaries. |
| (iv) | Volume Defects | **IMPLEMENTED**| [defect/volume_defects.py](file:///home/niel/git/crystal-mcp-server/src/python/generators/defect/volume_defects.py)| **Supported:** Voids (Sphere/Ellipsoid), Precipitates, Crack tips.<br>**Missing:** - |
| (v) | Amorphous | **PARTIAL** | [defect/amorphous.py](file:///home/niel/git/crystal-mcp-server/src/python/generators/defect/amorphous.py) | **Supported:** Random placement, Melt-quench simulation.<br>**Missing:** RMC / ML generators. |
| (vi) | Complexions | **MISSING** | - | **Missing:** Grain-boundary complexions (I-IV). |
| (vii) | Strain Fields | **IMPLEMENTED**| [external_fields/electric_fields.py](file:///home/niel/git/crystal-mcp-server/src/python/generators/external_fields/electric_fields.py)| **Supported:** Gaussian, Sinusoidal strain fields (via external fields).<br>**Missing:** Vegard strain / Misfit dislocation arrays. |
| (viii)| Chemical Disorder | **PARTIAL** | `defect/sqs.py:58` | **Supported:** SQS generation (random shuffling).<br>**Missing:** Cluster expansion (Code has placeholder 0.0 logic). |

### 8. Low-Dimensional
| ID | Catalogue Item | Status | Source (File:Line) | Capabilities & Gaps |
|:---|:---|:---|:---|:---|
| (i) | QDs | **IMPLEMENTED**| `quantum/dots.py:155` | **Supported:** Quantum dots (CdSe, InP etc).<br>**Missing:** Tetrapod / Platelet shapes. |
| (ii) | Core-Shell | **IMPLEMENTED**| `quantum/dots.py:347` | **Supported:** Core-shell dots with gradient alloys.<br>**Missing:** - |
| (iii) | Nanowires | **IMPLEMENTED**| [nanotube/nanowires.py](file:///home/niel/git/crystal-mcp-server/src/python/generators/nanotube/nanowires.py) | **Supported:** Nanowires (ZnO, GaN etc).<br>**Missing:** Branched / Kinked nanowires. |
| (iv) | Nanotubes | **IMPLEMENTED**| `nanotube/cnt.py:92` | **Supported:** SWCNT, MWCNT, Inorganic nanotubes (MoS2, BN).<br>**Missing:** Peapods (Verified in [biomolecules.py](file:///home/niel/git/crystal-mcp-server/src/python/generators/molecule/biomolecules.py) or missing?). |
| (v) | Nanoribbons | **IMPLEMENTED**| `two_d/ribbons.py` | **Supported:** Edge magnetism / reconstruction.<br>**Missing:** Edge doping. |
| (vi) | Nanoporous | **MISSING** | - | **Missing:** Nanoporous gold / de-alloyed ligaments. |
| (vii) | Plasmonics | **PARTIAL** | [photonic/metamaterials.py](file:///home/niel/git/crystal-mcp-server/src/python/generators/photonic/metamaterials.py) | **Supported:** Metasurfaces (SRR, V-antenna), Bow-tie.<br>**Missing:** Plasmonic dimers/trimers/nanostars explicitly. |
| (viii)| Topo Insulators | **IMPLEMENTED**| `electronic/topolo...` | **Supported:** Bi2Se3, Sb2Te3 nanostructures.<br>**Missing:** Surface state protection logic. |
| (ix) | Majorana | **IMPLEMENTED**| `electronic/topolo...` | **Supported:** Majorana nanowires (InSb/Al).<br>**Missing:** - |

### 9. Battery & Ionics
| ID | Catalogue Item | Status | Source (File:Line) | Capabilities & Gaps |
|:---|:---|:---|:---|:---|
| (i) | Cathodes | **IMPLEMENTED**| [battery/cathodes.py](file:///home/niel/git/crystal-mcp-server/src/python/generators/battery/cathodes.py) | **Supported:** Layered oxides (LiCoO2, NMC etc).<br>**Missing:** Li-rich Mn, O3/P2/P3 Na-layered. |
| (ii) | Electrolytes | **IMPLEMENTED**| [battery/electrolytes.py](file:///home/niel/git/crystal-mcp-server/src/python/generators/battery/electrolytes.py)| **Supported:** Solid Electrolytes (LLZO, LGPS, NASICON).<br>**Missing:** Tavorite, Argyrodite, Anti-perovskite. |
| (iii) | Anodes | **VERIFIED** | [battery/anodes.py](file:///home/niel/git/crystal-mcp-server/src/python/generators/battery/anodes.py) | **Supported:** Si, Ge, Sn, P, Metal (Li/Na/K/Zn) anodes.<br>**Missing:** - |
| (iv) | Sulfur | **VERIFIED** | [battery/anodes.py](file:///home/niel/git/crystal-mcp-server/src/python/generators/battery/anodes.py) | **Supported:** Sulfur cathodes in host matrices.<br>**Missing:** - |
| (v) | Metal Anodes | **VERIFIED** | [battery/anodes.py](file:///home/niel/git/crystal-mcp-server/src/python/generators/battery/anodes.py) | **Supported:** Li/Na/Zn metal anodes with SEI and dendrite seeds.<br>**Missing:** - |
| (vi) | Doping | **IMPLEMENTED**| [battery/cathodes.py](file:///home/niel/git/crystal-mcp-server/src/python/generators/battery/cathodes.py) | **Supported:** Gradient doping / Surface coating.<br>**Missing:** Ribbon doping. |
| (vii) | Concentration | **MISSING** | - | **Missing:** Concentration-cell electrolytes. |
| (viii)| Interfaces | **IMPLEMENTED**| `battery/electrolytes.py:255`| **Supported:** Cathode-Electrolyte interfaces.<br>**Missing:** Space-charge layers. |

### 10. Catalysis
| ID | Catalogue Item | Status | Source (File:Line) | Capabilities & Gaps |
|:---|:---|:---|:---|:---|
| (i) | SAA | **IMPLEMENTED**| `catalyst/saa.py:90` | **Supported:** Single Atom Alloys (Pt/Cu etc).<br>**Missing:** - |
| (ii) | HEA Surfaces | **IMPLEMENTED**| `catalyst/saa.py:192` | **Supported:** HEA surfaces with random/segregated termination.<br>**Missing:** Specific Cantor/Senkov sets. |
| (iii) | Oxide Clusters | **IMPLEMENTED**| [catalyst/clusters.py](file:///home/niel/git/crystal-mcp-server/src/python/generators/catalyst/clusters.py) | **Supported:** Supported clusters (Pt4, Au20 etc).<br>**Missing:** - |
| (iv) | Zeolite Clusters| **IMPLEMENTED**| [catalyst/clusters.py](file:///home/niel/git/crystal-mcp-server/src/python/generators/catalyst/clusters.py) | **Supported:** Clusters confined in zeolites.<br>**Missing:** - |
| (v) | SAC Supports | **IMPLEMENTED**| `catalyst/saa.py:283` | **Supported:** 2D SAC supports.<br>**Missing:** MN4/MN3C coordination motifs. |
| (vi) | Strain | **IMPLEMENTED**| `catalyst/saa.py:384` | **Supported:** Strain-engineered catalysts.<br>**Missing:** Ripple/Bending specific modes. |
| (vii) | Gating | **MISSING** | - | **Missing:** Electric-field/Ionic-liquid gating. |
| (viii)| Heterojunctions | **IMPLEMENTED**| `catalyst/clusters.py:448`| **Supported:** Photocatalyst heterojunctions (Type-I/II).<br>**Missing:** Z-scheme / S-scheme logic. |

### 11. High Pressure
| ID | Catalogue Item | Status | Source (File:Line) | Capabilities & Gaps |
|:---|:---|:---|:---|:---|
| (i) | HP Phases | **IMPLEMENTED**| [high_pressure/phases.py](file:///home/niel/git/crystal-mcp-server/src/python/generators/high_pressure/phases.py)| **Supported:** High pressure phase generation (Post-perovskite).<br>**Missing:** Aragonite specific. |
| (ii) | Polymeric N | **IMPLEMENTED**| [high_pressure/phases.py](file:///home/niel/git/crystal-mcp-server/src/python/generators/high_pressure/phases.py)| **Supported:** Polymeric Nitrogen (cg-N).<br>**Missing:** Black/Blue Phosphorus high-P phases. |
| (iii) | Super-ionic | **IMPLEMENTED**| [high_pressure/phases.py](file:///home/niel/git/crystal-mcp-server/src/python/generators/high_pressure/phases.py)| **Supported:** Super-ionic phases (Ice-X).<br>**Missing:** Li3N-LiOH super-ionic. |
| (iv) | WDM | **VERIFIED** | [external_fields/advanced_fields.py](file:///home/niel/git/crystal-mcp-server/src/python/generators/external_fields/advanced_fields.py) | **Supported:** Warm Dense Matter (hot electron) cells.<br>**Missing:** - |
| (v) | Laser Shock | **VERIFIED** | [external_fields/advanced_fields.py](file:///home/niel/git/crystal-mcp-server/src/python/generators/external_fields/advanced_fields.py) | **Supported:** Laser shocked cells (two-temperature model).<br>**Missing:** - |

### 12. Biomolecular
| ID | Catalogue Item | Status | Source (File:Line) | Capabilities & Gaps |
|:---|:---|:---|:---|:---|
| (i) | Lipid Bilayers | **IMPLEMENTED** | `biomolecular/lipids.py:66` | **Supported:** POPC, DPPC, DOPC, DLPC, DMPC + Cholesterol.<br>**Missing:** Cardiolipin, Ganglioside GM1. |
| (ii) | Memb. Proteins | **IMPLEMENTED** | `biomolecular/lipids.py:180` | **Supported:** GPCR, Ion channel, Porin insertion.<br>**Missing:** - |
| (iii)| Protein-Surface | **IMPLEMENTED** | `biomolecular/interfaces.py:43`| **Supported:** Fibronectin, Lysozyme, IgG, BSA on Ti, Au, Silica.<br>**Missing:** Collagen on Hydroxyapatite. |
| (iv) | DNA Origami | **IMPLEMENTED** | `biomolecular/interfaces.py:153`| **Supported:** Plates, Rods, Barrels on Graphene/MoS2.<br>**Missing:** - |
| (v) | Peptide Interface| **IMPLEMENTED** | `biomolecular/interfaces.py:232`| **Supported:** RGD, KLVFF, Diphenylalanine on TiO2/ZnO.<br>**Missing:** - |
| (vi) | Hydrogels | **MISSING** | - | **Missing:** Hydrogel networks (PEGDA/Polyacrylamide). |

### 13. Meta-structures
| ID | Catalogue Item | Status | Source (File:Line) | Capabilities & Gaps |
|:---|:---|:---|:---|:---|
| (i) | Moiré Lattices | **VERIFIED** | `twist/moire.py:41` | **Supported:** Moiré heterosuperlattices.<br>**Missing:** Generalized continuum model. |
| (ii) | Artificial | **VERIFIED** | [meta_structures/artificial_lattices.py](file:///home/niel/git/crystal-mcp-server/src/python/generators/meta_structures/artificial_lattices.py) | **Supported:** Kagome, Lieb, Checkerboard, Dice, Ruby, Honeycomb variants.<br>**Missing:** - |
| (iii) | Photonic | **IMPLEMENTED** | `photonic/metamaterials.py:288`| **Supported:** Woodpile, Inverse Opal, Triangular Hole.<br>**Missing:** - |
| (iv) | Metasurfaces | **VERIFIED** | `photonic/metamaterials.py:127`| **Supported:** V-antenna, SRR, Nano-fin, Bow-tie, Catenary.<br>**Missing:** - |
| (v) | Quasicrystals | **VERIFIED** | [meta_structures/quasicrystals.py](file:///home/niel/git/crystal-mcp-server/src/python/generators/meta_structures/quasicrystals.py) | **Supported:** Penrose, Fibonacci, Icosahedral, Decagonal.<br>**Missing:** - |

### 14. External Fields
| ID | Catalogue Item | Status | Source (File:Line) | Capabilities & Gaps |
|:---|:---|:---|:---|:---|
| (i) | Electric Fields | **VERIFIED** | [external_fields/electric_fields.py](file:///home/niel/git/crystal-mcp-server/src/python/generators/external_fields/electric_fields.py) | **Supported:** Electric-field polarised cells, Polar surfaces.<br>**Missing:** - |
| (ii) | Magnetic Fields | **VERIFIED** | [external_fields/advanced_fields.py](file:///home/niel/git/crystal-mcp-server/src/python/generators/external_fields/advanced_fields.py) | **Supported:** Skyrmions, Antiskyrmions, Bobbers, Bimerons.<br>**Missing:** - |
| (iii) | Optical Phonon | **VERIFIED** | [external_fields/advanced_fields.py](file:///home/niel/git/crystal-mcp-server/src/python/generators/external_fields/advanced_fields.py) | **Supported:** Optical-phonon pumped cells (Floquet-like).<br>**Missing:** - |
| (iv) | Strain Waves | **IMPLEMENTED** | `external_fields/strain...:107` | **Supported:** Periodic ripples, Uniaxial/Biaxial strain.<br>**Missing:** Shear/Bending dynamic waves. |
| (v) | Floquet | **VERIFIED** | [external_fields/advanced_fields.py](file:///home/niel/git/crystal-mcp-server/src/python/generators/external_fields/advanced_fields.py) | **Supported:** Floquet engineered cells (periodically driven).<br>**Missing:** - |

### 15. Output Formats
| ID | Catalogue Item | Status | Source (File:Line) | Capabilities & Gaps |
|:---|:---|:---|:---|:---|
| (i) | DFT Codes | **VERIFIED** | `output_formats/converters.py:46`| **Supported:** VASP, Quantum Espresso.<br>**Missing:** ABINIT, CASTEP, CP2K, etc. |
| (ii) | MD Codes | **VERIFIED** | `output_formats/converters.py:208`| **Supported:** LAMMPS.<br>**Missing:** AMBER, CHARMM, etc. |
| (iii) | Structure Files | **VERIFIED** | `output_formats/converters.py:333`| **Supported:** CIF, XYZ, PDB.<br>**Missing:** MOL, SDF, etc. |
| (iv) | Databases | **VERIFIED** | [output_formats/database_adapters.py](file:///home/niel/git/crystal-mcp-server/src/python/generators/output_formats/database_adapters.py) | **Supported:** AiiDA, ASE, AFLOW, Materials Project, OPTIMADE, JARVIS.<br>**Missing:** - |
| (v) | Visualization | **MISSING** | - | **Missing:** XSF, CUBE, VESTA specific files. (User identified next step). |

### 16. Quality Control
| ID | Catalogue Item | Status | Source (File:Line) | Capabilities & Gaps |
|:---|:---|:---|:---|:---|
| (i) | Symmetry | **IMPLEMENTED** | `quality_control/symmetry.py:73`| **Supported:** Symmetry analysis with tolerance sweep (PyXtal).<br>**Missing:** - |
| (ii) | Optimization | **IMPLEMENTED** | [quality_control/optimizer.py](file:///home/niel/git/crystal-mcp-server/src/python/generators/quality_control/optimizer.py) | **Supported:** Geometry optimizer (ASE), MLFF support (CHGNet, M3GNet, MACE).<br>**Note:** MLFF tests skipped due to missing dependencies/CUDA. |
| (iii) | Duplicates | **IMPLEMENTED** | `quality_control/validators.py:213`| **Supported:** Duplicate remover.<br>**Missing:** Graph-hash / Wyckoff-hash. |
| (iv) | Validation | **IMPLEMENTED** | `quality_control/validators.py:95`| **Supported:** Bond-valence checks, Relaxation convergence check.<br>**Missing:** RDF validator. |
| (v) | Charge | **IMPLEMENTED** | `quality_control/validators.py:176`| **Supported:** Charge/Stoichiometry checker.<br>**Missing:** - |
| (vi) | K-Path | **IMPLEMENTED** | `quality_control/kpath.py:90` | **Supported:** Automatic K-path generation.<br>**Missing:** - |
| (vii) | Tensors | **IMPLEMENTED** | [quality_control/optimizer.py](file:///home/niel/git/crystal-mcp-server/src/python/generators/quality_control/optimizer.py) | **Supported:** Elastic & Piezoelectric tensor templates.<br>**Missing:** - |
