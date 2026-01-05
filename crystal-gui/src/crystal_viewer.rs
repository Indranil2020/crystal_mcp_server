use serde::{Deserialize, Serialize};
use eframe::egui;
use std::f32::consts::PI;

// --- Data Structures ---

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
    #[serde(default)]
    pub cartesian: [f64; 3],
    #[serde(default)]
    pub wyckoff: Option<String>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SpaceGroup {
    pub number: u32,
    #[serde(default)]
    pub symbol: Option<String>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Metadata {
    #[serde(default)]
    pub formula: Option<String>,
    #[serde(default)]
    pub natoms: Option<u32>,
}

// --- Colors & Radii ---

pub struct ElementInfo;

impl ElementInfo {
    pub fn get_color(element: &str) -> egui::Color32 {
        let (r, g, b) = match element {
            "H" => (255, 255, 255),
            "C" => (105, 105, 105), // Dark grey
            "N" => (48, 80, 248),
            "O" => (255, 13, 13),
            "F" => (144, 224, 80),
            "Si" => (240, 200, 160),
            "P" => (255, 128, 0),
            "S" => (255, 255, 48),
            "Cl" => (31, 240, 31),
            "Br" => (166, 41, 41),
            "I" => (148, 0, 148),
            "Fe" => (224, 102, 51),
            "Cu" => (200, 115, 51),
            "Zn" => (125, 128, 176),
            "Au" => (255, 209, 35),
            "Ag" => (192, 192, 192),
            "Mg" => (0, 100, 0),
            "Ca" => (61, 255, 0),
            "Li" => (204, 128, 255),
            "Na" => (171, 92, 242),
            "K" => (143, 64, 212),
            _ => (255, 192, 203), // Pink default
        };
        egui::Color32::from_rgb(r, g, b)
    }

    pub fn get_radius(element: &str) -> f32 {
        match element {
            "H" => 0.31,
            "C" => 0.76,
            "N" => 0.71,
            "O" => 0.66,
            "F" => 0.57,
            "Si" => 1.11,
            "P" => 1.07,
            "S" => 1.05,
            "Cl" => 0.99,
            "Br" => 1.14,
            "I" => 1.33,
            "Fe" => 1.25,
            "Cu" => 1.28,
            "Zn" => 1.22,
            "Au" => 1.44,
            "Ag" => 1.45,
            "Li" => 1.28,
            "Na" => 1.66,
            "K" => 2.03,
            "Mg" => 1.41,
            "Ca" => 1.76,
            _ => 1.0,
        }
    }
}

// --- Viewer ---

pub struct CrystalViewer {
    pub structure: Option<CrystalStructure>,
    
    // Camera state
    rotation: [f32; 2], // X and Y rotation
    zoom: f32,
    pan: [f32; 2],
    
    // Settings
    pub show_bonds: bool,
    pub show_unit_cell: bool,
    pub atom_scale: f32,
}

impl Default for CrystalViewer {
    fn default() -> Self {
        Self {
            structure: None,
            rotation: [0.0, 0.0],
            zoom: 1.0,
            pan: [0.0, 0.0],
            show_bonds: true,
            show_unit_cell: true,
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
        self.rotation = [0.0, 0.0];
        self.zoom = 1.0;
        self.pan = [0.0, 0.0];
    }
    
    pub fn rotate(&mut self, dx: f32, dy: f32) {
        self.rotation[0] += dx * 0.01;
        self.rotation[1] += dy * 0.01;
    }

    pub fn zoom_in(&mut self) {
        self.zoom *= 1.1;
    }

    pub fn zoom_out(&mut self) {
        self.zoom /= 1.1;
    }
    
    pub fn get_structure_info(&self) -> Option<String> {
        self.structure.as_ref().map(|s| {
            let formula = s.metadata.as_ref()
                .and_then(|m| m.formula.clone())
                .unwrap_or_else(|| "Unknown".to_string());
            let natoms = s.metadata.as_ref()
                .and_then(|m| m.natoms)
                .unwrap_or(0);
            
            format!("Formula: {}\nAtoms: {}", formula, natoms)
        })
    }
    
    /// Detect bonds based on distance
    fn detect_bonds(atoms: &[Atom]) -> Vec<(usize, usize)> {
        let mut bonds = Vec::new();
        // Covalent radii sum * tolerance
        let tolerance = 1.25; 
        
        for i in 0..atoms.len() {
            for j in (i+1)..atoms.len() {
                let p1 = atoms[i].cartesian;
                let p2 = atoms[j].cartesian;
                let elem1 = &atoms[i].element;
                let elem2 = &atoms[j].element;
                
                let dist_sq = (p1[0]-p2[0]).powi(2) + (p1[1]-p2[1]).powi(2) + (p1[2]-p2[2]).powi(2);
                let dist = dist_sq.sqrt();
                
                let r1 = ElementInfo::get_radius(elem1);
                let r2 = ElementInfo::get_radius(elem2);
                let max_dist = (r1 + r2) as f64 * tolerance;
                
                if dist < max_dist {
                    bonds.push((i, j));
                }
            }
        }
        bonds
    }

    pub fn render(&mut self, ui: &mut egui::Ui) {
        let available_size = ui.available_size();
        let (response, painter) = ui.allocate_painter(available_size, egui::Sense::drag());
        
        let rect = response.rect;
        let center = rect.center();
        
        // Background
        painter.rect_filled(rect, 0.0, egui::Color32::from_rgb(20, 20, 30));
        
        // Input Handling
        if response.dragged() {
            let delta = response.drag_delta();
            // Left click rotate, Right/Middle click pan (simplified: just rotate for now)
            self.rotate(delta.x, delta.y);
        }
        
        if response.hovered() {
            let scroll = ui.input(|i| i.raw_scroll_delta.y);
            if scroll > 0.0 { self.zoom_in(); }
            else if scroll < 0.0 { self.zoom_out(); }
        }
        
        let Some(structure) = &self.structure else {
            painter.text(
                center,
                egui::Align2::CENTER_CENTER,
                "No structure loaded",
                egui::FontId::proportional(20.0),
                egui::Color32::GRAY
            );
            return;
        };
        
        let atoms = structure.atoms.as_slice();
        let atoms = if atoms.is_empty() {
             structure.sites.as_deref().unwrap_or(&[])
        } else {
             atoms
        };
        
        if atoms.is_empty() { return; }

        // 3D Projection Logic
        // Rotate points
        let rot_x = self.rotation[1]; // Pitch
        let rot_y = self.rotation[0]; // Yaw
        
        let cos_x = rot_x.cos() as f64;
        let sin_x = rot_x.sin() as f64;
        let cos_y = rot_y.cos() as f64;
        let sin_y = rot_y.sin() as f64;
        
        // View Scale
        let scale = 30.0 * self.zoom;
        
        struct RenderItem {
            z_depth: f64,
            draw_fn: Box<dyn FnOnce(&egui::Painter)>,
        }
        
        let mut render_queue: Vec<RenderItem> = Vec::new();
        
        let project = |p: [f64; 3]| -> ([f64; 3], egui::Pos2) {
            // Rotate Y
            let x1 = p[0] * cos_y - p[2] * sin_y;
            let z1 = p[0] * sin_y + p[2] * cos_y;
            // Rotate X
            let y2 = p[1] * cos_x - z1 * sin_x;
            let z2 = p[1] * sin_x + z1 * cos_x;
            
            // Perspective (optional, using Ortho for simplicity/clarity often preferred for schematics)
            // Let's stick to weak perspective / orthographic for correct sizing
            let screen_x = center.x as f64 + x1 * scale as f64 + self.pan[0] as f64;
            let screen_y = center.y as f64 - y2 * scale as f64 + self.pan[1] as f64;
            
            ([x1, y2, z2], egui::Pos2::new(screen_x as f32, screen_y as f32))
        };
        
        // 1. Bonds
        if self.show_bonds {
            let bonds = Self::detect_bonds(atoms);
            for (idx1, idx2) in bonds {
                let p1_raw = atoms[idx1].cartesian;
                let p2_raw = atoms[idx2].cartesian;
                
                let (p1_rot, pos1) = project(p1_raw);
                let (p2_rot, pos2) = project(p2_raw);
                
                // Average Z for sorting
                let z_depth = (p1_rot[2] + p2_rot[2]) / 2.0;
                
                // Bond is a line
                let bond_width = 4.0 * self.zoom * self.atom_scale;
                let bond_color = egui::Color32::from_gray(180);
                
                render_queue.push(RenderItem {
                    z_depth,
                    draw_fn: Box::new(move |painter| {
                        painter.line_segment([pos1, pos2], egui::Stroke::new(bond_width, bond_color));
                        // Add highlights for cylinder effect? Too expensive/complex for 2D lines.
                        // Simple shading: darker if further back?
                    }),
                });
            }
        }
        
        // 2. Atoms (Glossy Spheres)
        for atom in atoms {
            let (p_rot, pos) = project(atom.cartesian);
            let radius = ElementInfo::get_radius(&atom.element) * 10.0 * self.atom_scale as f64 * self.zoom as f64;
            let color = ElementInfo::get_color(&atom.element);
            let z_depth = p_rot[2];
            
            render_queue.push(RenderItem {
                z_depth,
                draw_fn: Box::new(move |painter| {
                    // Draw Glossy Sphere
                    
                    // 1. Base circle (darker edge)
                    painter.circle_filled(pos, radius as f32, color);
                    
                    // 2. Outline/Shading (slight border)
                    painter.circle_stroke(pos, radius as f32, egui::Stroke::new(1.0, egui::Color32::BLACK));
                    
                    // 3. Highlight (Specular reflection)
                    // Offset highlight to top-left
                    let highlight_offset = radius as f32 * 0.3;
                    let highlight_pos = egui::Pos2::new(pos.x - highlight_offset, pos.y - highlight_offset);
                    let highlight_radius = radius as f32 * 0.25;
                    let highlight_color = egui::Color32::from_white_alpha(180); // Semi-transparent white
                    
                    painter.circle_filled(highlight_pos, highlight_radius, highlight_color);
                    
                    // 4. Element Symbol Text
                    if radius > 8.0 {
                        painter.text(
                            pos, 
                            egui::Align2::CENTER_CENTER, 
                            &atom.element, 
                            egui::FontId::proportional(radius as f32), 
                            egui::Color32::BLACK
                        );
                    }
                })
            });
        }
        
        // Sort by Z (painters algorithm: draw furthest first)
        // Z increases towards viewer? p_rot calculated: 
        // z2 = p[1]*sin + z1*cos. 
        // Standard handedness: Z comes out of screen.
        // So smaller Z (more negative) is further away.
        render_queue.sort_by(|a, b| a.z_depth.partial_cmp(&b.z_depth).unwrap_or(std::cmp::Ordering::Equal));
        
        // Render
        for item in render_queue {
            (item.draw_fn)(&painter);
        }
    }
}
