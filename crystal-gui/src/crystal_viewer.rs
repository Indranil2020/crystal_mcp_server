//! 3D Crystal Structure Viewer
//! 
//! Renders crystal structures using three-d for 3D visualization.
//! Similar to Materials Project visualization style.

use serde::{Deserialize, Serialize};

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct CrystalStructure {
    pub lattice: Lattice,
    #[serde(default)]
    pub atoms: Vec<Atom>,
    #[serde(default)]
    pub sites: Option<Vec<Atom>>,
    pub space_group: SpaceGroup,
    #[serde(default)]
    pub metadata: Option<Metadata>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Lattice {
    pub a: f64,
    pub b: f64,
    pub c: f64,
    pub alpha: f64,
    pub beta: f64,
    pub gamma: f64,
    #[serde(default)]
    pub matrix: Option<Vec<Vec<f64>>>,
    #[serde(default)]
    pub volume: Option<f64>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Atom {
    pub element: String,
    pub coords: [f64; 3],
    pub cartesian: [f64; 3],
    #[serde(default)]
    pub wyckoff: Option<String>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SpaceGroup {
    pub number: u32,
    #[serde(default)]
    pub symbol: Option<String>,
    #[serde(default)]
    pub crystal_system: Option<String>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Metadata {
    #[serde(default)]
    pub formula: Option<String>,
    #[serde(default)]
    pub natoms: Option<u32>,
    #[serde(default)]
    pub volume: Option<f64>,
    #[serde(default)]
    pub density: Option<f64>,
}

pub struct ElementColors;

impl ElementColors {
    pub fn get_color(element: &str) -> [f32; 3] {
        match element {
            "H" => [1.0, 1.0, 1.0],
            "C" => [0.3, 0.3, 0.3],
            "N" => [0.0, 0.0, 1.0],
            "O" => [1.0, 0.0, 0.0],
            "F" => [0.0, 1.0, 0.0],
            "Si" => [0.5, 0.5, 0.8],
            "P" => [1.0, 0.5, 0.0],
            "S" => [1.0, 1.0, 0.0],
            "Cl" => [0.0, 1.0, 0.0],
            "Fe" => [0.8, 0.4, 0.0],
            "Cu" => [0.8, 0.5, 0.2],
            "Zn" => [0.5, 0.5, 0.5],
            "Al" => [0.7, 0.7, 0.8],
            "Na" => [0.7, 0.0, 0.7],
            "K" => [0.5, 0.0, 0.5],
            "Ca" => [0.0, 0.8, 0.0],
            "Ti" => [0.6, 0.6, 0.6],
            "Mg" => [0.0, 0.6, 0.0],
            "Au" => [1.0, 0.84, 0.0],
            "Ag" => [0.75, 0.75, 0.75],
            _ => [0.5, 0.5, 0.5],
        }
    }

    pub fn get_radius(element: &str) -> f32 {
        match element {
            "H" => 0.31,
            "C" => 0.77,
            "N" => 0.71,
            "O" => 0.66,
            "F" => 0.57,
            "Si" => 1.11,
            "P" => 1.07,
            "S" => 1.05,
            "Cl" => 1.02,
            "Fe" => 1.26,
            "Cu" => 1.28,
            "Zn" => 1.22,
            "Al" => 1.21,
            "Na" => 1.66,
            "K" => 2.03,
            "Ca" => 1.76,
            "Ti" => 1.47,
            "Mg" => 1.41,
            "Au" => 1.44,
            "Ag" => 1.45,
            _ => 1.0,
        }
    }
}

pub struct CrystalViewer {
    pub structure: Option<CrystalStructure>,
    pub rotation: [f32; 3],
    pub zoom: f32,
    pub show_unit_cell: bool,
    pub show_bonds: bool,
    pub atom_scale: f32,
}

impl Default for CrystalViewer {
    fn default() -> Self {
        Self {
            structure: None,
            rotation: [0.0, 0.0, 0.0],
            zoom: 1.0,
            show_unit_cell: true,
            show_bonds: true,
            atom_scale: 0.5,
        }
    }
}

impl CrystalViewer {
    pub fn new() -> Self {
        Self::default()
    }

    pub fn set_structure(&mut self, structure: CrystalStructure) {
        self.structure = Some(structure);
        self.reset_view();
    }

    pub fn clear_structure(&mut self) {
        self.structure = None;
    }

    pub fn reset_view(&mut self) {
        self.rotation = [0.0, 0.0, 0.0];
        self.zoom = 1.0;
    }

    pub fn rotate(&mut self, dx: f32, dy: f32) {
        self.rotation[0] += dy * 0.01;
        self.rotation[1] += dx * 0.01;
    }

    pub fn zoom_in(&mut self) {
        self.zoom *= 1.1;
    }

    pub fn zoom_out(&mut self) {
        self.zoom /= 1.1;
    }

    pub fn get_atom_positions(&self) -> Vec<([f32; 3], String, [f32; 3], f32)> {
        let mut positions = Vec::new();
        
        if let Some(ref structure) = self.structure {
            // Use atoms if available, otherwise try sites
            let atoms = if !structure.atoms.is_empty() {
                &structure.atoms
            } else if let Some(ref sites) = structure.sites {
                sites
            } else {
                return positions;
            };
            
            for atom in atoms {
                let pos = [
                    atom.cartesian[0] as f32,
                    atom.cartesian[1] as f32,
                    atom.cartesian[2] as f32,
                ];
                let color = ElementColors::get_color(&atom.element);
                let radius = ElementColors::get_radius(&atom.element) * self.atom_scale;
                positions.push((pos, atom.element.clone(), color, radius));
            }
        }
        
        positions
    }

    pub fn get_unit_cell_lines(&self) -> Vec<([f32; 3], [f32; 3])> {
        let mut lines = Vec::new();
        
        if let Some(ref structure) = self.structure {
            // Get lattice vectors from matrix or calculate from parameters
            let (a, b, c) = if let Some(ref m) = structure.lattice.matrix {
                if m.len() >= 3 && m[0].len() >= 3 && m[1].len() >= 3 && m[2].len() >= 3 {
                    (
                        [m[0][0] as f32, m[0][1] as f32, m[0][2] as f32],
                        [m[1][0] as f32, m[1][1] as f32, m[1][2] as f32],
                        [m[2][0] as f32, m[2][1] as f32, m[2][2] as f32],
                    )
                } else {
                    // Fallback to simple cubic
                    let la = structure.lattice.a as f32;
                    let lb = structure.lattice.b as f32;
                    let lc = structure.lattice.c as f32;
                    ([la, 0.0, 0.0], [0.0, lb, 0.0], [0.0, 0.0, lc])
                }
            } else {
                // No matrix, use lattice parameters (assume orthorhombic for simplicity)
                let la = structure.lattice.a as f32;
                let lb = structure.lattice.b as f32;
                let lc = structure.lattice.c as f32;
                ([la, 0.0, 0.0], [0.0, lb, 0.0], [0.0, 0.0, lc])
            };
            let o = [0.0f32, 0.0, 0.0];
            
            fn add(p1: [f32; 3], p2: [f32; 3]) -> [f32; 3] {
                [p1[0] + p2[0], p1[1] + p2[1], p1[2] + p2[2]]
            }
            
            lines.push((o, a));
            lines.push((o, b));
            lines.push((o, c));
            lines.push((a, add(a, b)));
            lines.push((a, add(a, c)));
            lines.push((b, add(b, a)));
            lines.push((b, add(b, c)));
            lines.push((c, add(c, a)));
            lines.push((c, add(c, b)));
            lines.push((add(a, b), add(add(a, b), c)));
            lines.push((add(a, c), add(add(a, c), b)));
            lines.push((add(b, c), add(add(b, c), a)));
        }
        
        lines
    }

    pub fn get_structure_info(&self) -> Option<String> {
        self.structure.as_ref().map(|s| {
            let formula = s.metadata.as_ref()
                .and_then(|m| m.formula.as_ref())
                .map(|f| f.as_str())
                .unwrap_or("Unknown");
            let symbol = s.space_group.symbol.as_deref().unwrap_or("?");
            let natoms = s.metadata.as_ref()
                .and_then(|m| m.natoms)
                .unwrap_or(0);
            let volume = s.lattice.volume.unwrap_or(0.0);
            
            format!(
                "Formula: {}\nSpace Group: {} ({})\nAtoms: {}\nVolume: {:.2} Å³\nLattice: a={:.3} b={:.3} c={:.3}\nAngles: α={:.1}° β={:.1}° γ={:.1}°",
                formula,
                symbol,
                s.space_group.number,
                natoms,
                volume,
                s.lattice.a, s.lattice.b, s.lattice.c,
                s.lattice.alpha, s.lattice.beta, s.lattice.gamma
            )
        })
    }
}
