//! Main Application State and UI
//! 
//! Implements the egui-based GUI with chat interface and crystal visualization.

use eframe::egui;
use serde_json::{json, Value};
use std::sync::{Arc, Mutex};

use crate::crystal_viewer::{CrystalStructure, CrystalViewer};
use crate::llm_client::{ChatMessage, LlmClient};
use crate::mcp_client::{McpClient, Tool};
use crate::structure_export::save_structure;
use rfd::FileDialog;

pub struct CrystalApp {
    mcp_client: Arc<Mutex<Option<McpClient>>>,
    llm_client: Arc<Mutex<LlmClient>>,
    crystal_viewer: CrystalViewer,
    
    chat_history: Vec<ChatMessage>,
    chat_input: String,
    
    status_message: String,
    mcp_connected: bool,
    llm_connected: bool,
    
    available_tools: Vec<Tool>,
    selected_tool: Option<String>,
    tool_params: String,
    
    show_settings: bool,
    ollama_url: String,
    ollama_model: String,
    mcp_server_path: String,
    current_structure_json: Option<Value>,
}

impl CrystalApp {
    pub fn new(_cc: &eframe::CreationContext<'_>) -> Self {
        // Use absolute path - can be overridden in settings
        let mcp_server_path = "/home/niel/git/crystal-mcp-server/dist/index.js".to_string();

        Self {
            mcp_client: Arc::new(Mutex::new(None)),
            llm_client: Arc::new(Mutex::new(LlmClient::default())),
            crystal_viewer: CrystalViewer::new(),
            
            chat_history: vec![ChatMessage {
                role: "system".to_string(),
                content: "You are a helpful assistant for crystal structure generation. Help users create and analyze crystal structures.".to_string(),
            }],
            chat_input: String::new(),
            
            status_message: "Ready. Connect to MCP server and Ollama to begin.".to_string(),
            mcp_connected: false,
            llm_connected: false,
            
            available_tools: Vec::new(),
            selected_tool: None,
            tool_params: String::new(),
            
            show_settings: false,
            ollama_url: "http://localhost:11434".to_string(),
            ollama_model: "qwen2.5:7b".to_string(),
            mcp_server_path,
            current_structure_json: None,
        }
    }

    fn connect_mcp(&mut self) {
        let path = self.mcp_server_path.clone();
        let mut client = McpClient::new(path);
        
        match client.start() {
            Ok(()) => {
                self.available_tools = client.get_tools().to_vec();
                self.mcp_connected = true;
                self.status_message = format!("MCP connected. {} tools available.", self.available_tools.len());
                *self.mcp_client.lock().unwrap() = Some(client);
            }
            Err(e) => {
                self.status_message = format!("MCP connection failed: {}", e);
                self.mcp_connected = false;
            }
        }
    }

    fn disconnect_mcp(&mut self) {
        if let Some(mut client) = self.mcp_client.lock().unwrap().take() {
            client.stop();
        }
        self.mcp_connected = false;
        self.available_tools.clear();
        self.status_message = "MCP disconnected.".to_string();
    }

    fn connect_llm(&mut self) {
        let mut client = self.llm_client.lock().unwrap();
        client.set_model(&self.ollama_model);
        
        match client.check_connection() {
            Ok(true) => {
                self.llm_connected = true;
                let models = client.get_available_models().len();
                self.status_message = format!("Ollama connected. {} models available.", models);
            }
            Ok(false) | Err(_) => {
                self.llm_connected = false;
                self.status_message = "Ollama connection failed. Is Ollama running?".to_string();
            }
        }
    }

    fn format_tools_for_llm(&self) -> String {
        let mut output = String::new();
        for tool in &self.available_tools {
            output.push_str(&format!("## {}\n", tool.name));
            if let Some(ref desc) = tool.description {
                output.push_str(&format!("{}\n", desc));
            }
            if let Some(ref schema) = tool.input_schema {
                if let Some(props) = schema.get("properties") {
                    output.push_str("Parameters:\n");
                    if let Some(obj) = props.as_object() {
                        for (key, val) in obj {
                            let type_str = val.get("type").and_then(|t| t.as_str()).unwrap_or("any");
                            let desc_str = val.get("description").and_then(|d| d.as_str()).unwrap_or("");
                            output.push_str(&format!("  - {}: {} - {}\n", key, type_str, desc_str));
                        }
                    }
                }
                if let Some(required) = schema.get("required").and_then(|r| r.as_array()) {
                    let req_strs: Vec<&str> = required.iter().filter_map(|v| v.as_str()).collect();
                    if !req_strs.is_empty() {
                        output.push_str(&format!("Required: {}\n", req_strs.join(", ")));
                    }
                }
            }
            output.push('\n');
        }
        output
    }

    fn get_tool_names(&self) -> Vec<String> {
        self.available_tools.iter().map(|t| t.name.clone()).collect()
    }

    fn send_chat_message(&mut self) {
        if self.chat_input.trim().is_empty() {
            return;
        }

        let user_message = ChatMessage {
            role: "user".to_string(),
            content: self.chat_input.clone(),
        };
        self.chat_history.push(user_message);
        self.chat_input.clear();

        if !self.llm_connected {
            self.connect_llm();
        }
        if !self.llm_connected {
            self.chat_history.push(ChatMessage {
                role: "assistant".to_string(),
                content: "LLM not connected. Is Ollama running?".to_string(),
            });
            return;
        }

        if !self.mcp_connected {
            self.connect_mcp();
        }
        if !self.mcp_connected {
            self.chat_history.push(ChatMessage {
                role: "assistant".to_string(),
                content: "MCP not connected. Please check MCP server path in Settings.".to_string(),
            });
            return;
        }

        let tools_desc = self.format_tools_for_llm();
        let tool_names = self.get_tool_names();
        let client = self.llm_client.lock().unwrap();

        match client.chat_with_tools(&self.chat_history, &tools_desc) {
            Ok(response) => {
                // Try to parse tool call from response
                let tool_call = self.parse_tool_call(&response);
                
                if let Some((tool_name, params)) = tool_call {
                    if !tool_names.iter().any(|t| t == &tool_name) {
                        self.chat_history.push(ChatMessage {
                            role: "assistant".to_string(),
                            content: format!(
                                "Tool `{}` is not available. Available tools:\n{}",
                                tool_name,
                                tool_names.join("\n")
                            ),
                        });
                        return;
                    }

                    // Auto-fill tool parameters
                    self.selected_tool = Some(tool_name.clone());
                    self.tool_params = serde_json::to_string_pretty(&params).unwrap_or_default();
                    
                    self.chat_history.push(ChatMessage {
                        role: "assistant".to_string(),
                        content: format!("Executing `{}`...", tool_name),
                    });
                    
                    // Auto-execute if MCP is connected
                    drop(client); // Release lock before calling tool
                    if self.mcp_connected {
                        self.call_tool();
                    }
                } else {
                    self.chat_history.push(ChatMessage {
                        role: "assistant".to_string(),
                        content: response,
                    });
                }
            }
            Err(e) => {
                self.chat_history.push(ChatMessage {
                    role: "assistant".to_string(),
                    content: format!("Error: {}", e),
                });
            }
        }
    }
    
    fn parse_tool_call(&self, response: &str) -> Option<(String, Value)> {
        // Try to extract JSON from response
        let json_str = if let Some(start) = response.find("```json") {
            let start = start + 7;
            if let Some(end) = response[start..].find("```") {
                response[start..start + end].trim()
            } else {
                return None;
            }
        } else if let Some(start) = response.find('{') {
            if let Some(end) = response.rfind('}') {
                &response[start..=end]
            } else {
                return None;
            }
        } else {
            return None;
        };
        
        // Parse JSON
        let parsed: Value = serde_json::from_str(json_str).ok()?;
        let tool_name = parsed.get("tool")?.as_str()?.to_string();
        let params = parsed.get("params")?.clone();
        
        Some((tool_name, params))
    }

    fn call_tool(&mut self) {
        let Some(ref tool_name) = self.selected_tool else {
            self.status_message = "No tool selected. Please select a tool first.".to_string();
            return;
        };

        let params: Value = match serde_json::from_str(&self.tool_params) {
            Ok(p) => p,
            Err(e) => {
                self.status_message = format!("Invalid JSON params: {}", e);
                return;
            }
        };
        
        self.status_message = format!("Executing {}...", tool_name);
        
        // Call MCP tool and get result
        let tool_result = {
            let mut client_guard = self.mcp_client.lock().unwrap();
            let Some(ref mut client) = *client_guard else {
                self.status_message = "MCP not connected. Click 'Connect MCP' first.".to_string();
                return;
            };
            client.call_tool(tool_name, params)
        }; // client_guard dropped here
        
        match tool_result {
            Ok(result) => {
                self.status_message = format!("Tool {} executed.", tool_name);
                
                // Check for errors first
                if result.is_error {
                    if let Some(content) = result.content.first() {
                        self.status_message = format!("Tool error: {}", content.text);
                        self.chat_history.push(ChatMessage {
                            role: "assistant".to_string(),
                            content: self.status_message.clone(),
                        });
                    }
                    return;
                }

                for content in &result.content {
                    if content.content_type != "text" {
                        continue;
                    }
                    let trimmed = content.text.trim();
                    if trimmed.starts_with("<json-data>") && trimmed.ends_with("</json-data>") {
                        continue;
                    }
                    if trimmed.contains("<json-data>") && trimmed.contains("</json-data>") {
                        continue;
                    }

                    self.chat_history.push(ChatMessage {
                        role: "assistant".to_string(),
                        content: content.text.clone(),
                    });
                }
                
                // Try to extract structure from response
                let mut found_structure = false;
                let mut parse_error: Option<String> = None;
                
                for content in &result.content {
                    if content.content_type == "text" {
                        // Check for <json-data> tag format from MCP server
                        if content.text.contains("<json-data>") {
                            if let Some(start) = content.text.find("<json-data>") {
                                let json_start = start + "<json-data>".len();
                                if let Some(end) = content.text[json_start..].find("</json-data>") {
                                    let json_str = content.text[json_start..json_start + end].trim();
                                    if let Ok(wrapper) = serde_json::from_str::<Value>(json_str) {
                                        if let Some(structure_val) = wrapper.get("structure").or_else(|| wrapper.get("data").and_then(|d| d.get("structure"))) {
                                            self.current_structure_json = Some(structure_val.clone());
                                            match serde_json::from_value::<CrystalStructure>(structure_val.clone()) {
                                                Ok(structure) => {
                                                    self.crystal_viewer.set_structure(structure);
                                                    found_structure = true;
                                                    break;
                                                }
                                                Err(e) => {
                                                    parse_error = Some(e.to_string());
                                                    break;
                                                }
                                            }
                                        }
                                        if let Some(slab_val) = wrapper.get("slab").or_else(|| wrapper.get("data").and_then(|d| d.get("slab"))) {
                                            self.current_structure_json = Some(slab_val.clone());
                                            match serde_json::from_value::<CrystalStructure>(slab_val.clone()) {
                                                Ok(structure) => {
                                                    self.crystal_viewer.set_structure(structure);
                                                    found_structure = true;
                                                    break;
                                                }
                                                Err(e) => {
                                                    parse_error = Some(e.to_string());
                                                    break;
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        }
                        else if let Ok(structure) = serde_json::from_str::<CrystalStructure>(&content.text) {
                            self.current_structure_json = Some(serde_json::to_value(&structure).unwrap_or_default());
                            self.crystal_viewer.set_structure(structure);
                            found_structure = true;
                            break;
                        }
                        else if content.text.contains("\"structure\"") {
                            if let Ok(wrapper) = serde_json::from_str::<Value>(&content.text) {
                                if let Some(structure_val) = wrapper.get("structure").or_else(|| wrapper.get("data").and_then(|d| d.get("structure"))) {
                                    self.current_structure_json = Some(structure_val.clone());
                                    match serde_json::from_value::<CrystalStructure>(structure_val.clone()) {
                                        Ok(structure) => {
                                            self.crystal_viewer.set_structure(structure);
                                            found_structure = true;
                                            break;
                                        }
                                        Err(e) => {
                                            parse_error = Some(e.to_string());
                                            break;
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
                
                if found_structure {
                    self.status_message = "Structure generated. View in Crystal Viewer panel.".to_string();
                } else if let Some(err) = parse_error {
                    self.status_message = format!("Structure payload received but failed to load in viewer: {}", err);
                    self.chat_history.push(ChatMessage {
                        role: "assistant".to_string(),
                        content: self.status_message.clone(),
                    });
                } else if !result.content.is_empty() {
                    self.status_message = format!("Tool {} completed (no structure to display).", tool_name);
                }
            }
            Err(e) => {
                self.status_message = format!("Tool error: {}", e);
                self.chat_history.push(ChatMessage {
                    role: "assistant".to_string(),
                    content: self.status_message.clone(),
                });
            }
        }
    }
    
    fn generate_and_show_visualization(&mut self) {
        let Some(ref structure) = self.current_structure_json else {
            self.status_message = "No structure to visualize.".to_string();
            return;
        };
        
        // Save structure to temp file
        let input_path = "/tmp/crystal_viz_input.json";
        let output_path = "/tmp/crystal_structure.png";
        
        let json_str = match serde_json::to_string(structure) {
            Ok(s) => s,
            Err(e) => {
                self.status_message = format!("Failed to serialize: {}", e);
                return;
            }
        };
        
        if let Err(e) = std::fs::write(input_path, &json_str) {
            self.status_message = format!("Failed to write input file: {}", e);
            return;
        }
        
        // Generate PNG using inline Python script (ASE)
        let python_script = format!(r#"
import json
import numpy as np
from ase import Atoms
from ase.visualize.plot import plot_atoms
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

with open('{}') as f:
    structure = json.load(f)

lattice = structure.get('lattice', {{}})
atoms_data = structure.get('atoms', structure.get('sites', []))

if not lattice or not atoms_data:
    print("Invalid structure")
    exit(1)

cell = np.array(lattice.get('matrix', [[10,0,0],[0,10,0],[0,0,10]]))
symbols = [a.get('element', 'X') for a in atoms_data]
positions = [a.get('cartesian', a.get('coords', [0,0,0])) for a in atoms_data]

atoms = Atoms(symbols=symbols, positions=positions, cell=cell, pbc=True)

fig, ax = plt.subplots(figsize=(10, 10))
ax.axis('off')
ax.set_title(f"{{len(symbols)}} atoms: {{', '.join(set(symbols))}}", fontsize=14)
plot_atoms(atoms, ax, radii=0.8, rotation='10x,10y,10z')
fig.savefig('{}', bbox_inches='tight', dpi=150, facecolor='white')
plt.close()
print("OK")
"#, input_path, output_path);
        
        let output = std::process::Command::new("python3")
            .arg("-c")
            .arg(&python_script)
            .output();
            
        match output {
            Ok(o) => {
                let stdout = String::from_utf8_lossy(&o.stdout);
                if !stdout.contains("OK") {
                    let stderr = String::from_utf8_lossy(&o.stderr);
                    self.status_message = format!("Visualization failed: {}", stderr);
                    return;
                }
            }
            Err(e) => {
                self.status_message = format!("Failed to run visualization: {}", e);
                return;
            }
        }
        
        // Open PNG with image viewer
        if std::path::Path::new(output_path).exists() {
            let _ = std::process::Command::new("xdg-open")
                .arg(output_path)
                .spawn();
            self.status_message = "Structure visualization opened.".to_string();
        } else {
            self.status_message = "Visualization file not created.".to_string();
        }
    }

    fn render_menu_bar(&mut self, ui: &mut egui::Ui) {
        egui::menu::bar(ui, |ui| {
            ui.menu_button("File", |ui| {
                if ui.button("Settings").clicked() {
                    self.show_settings = true;
                    ui.close_menu();
                }
                if ui.button("Save Structure...").clicked() {
                    if let Some(structure) = &self.crystal_viewer.structure {
                        if let Some(path) = FileDialog::new()
                            .add_filter("CIF", &["cif"])
                            .add_filter("POSCAR", &["vasp", "poscar"])
                            .add_filter("XYZ", &["xyz"])
                            .add_filter("JSON", &["json"])
                            .save_file() 
                        {
                            let ext = path.extension().and_then(|e| e.to_str()).unwrap_or("cif");
                            match save_structure(structure, path.to_str().unwrap(), ext) {
                                Ok(_) => self.status_message = format!("Saved to {:?}", path),
                                Err(e) => self.status_message = format!("Failed to save: {}", e),
                            }
                        }
                    } else {
                        self.status_message = "No structure to save.".to_string();
                    }
                    ui.close_menu();
                }
                ui.separator();
                if ui.button("Exit").clicked() {
                    std::process::exit(0);
                }
            });
            
            ui.menu_button("Connection", |ui| {
                if self.mcp_connected {
                    if ui.button("Disconnect MCP").clicked() {
                        self.disconnect_mcp();
                        ui.close_menu();
                    }
                } else if ui.button("Connect MCP").clicked() {
                    self.connect_mcp();
                    ui.close_menu();
                }
                
                ui.separator();
                
                if self.llm_connected {
                    ui.label("✓ Ollama Connected");
                } else if ui.button("Connect Ollama").clicked() {
                    self.connect_llm();
                    ui.close_menu();
                }
            });
            
            ui.menu_button("View", |ui| {
                if ui.button("Reset View").clicked() {
                    self.crystal_viewer.reset_view();
                    ui.close_menu();
                }
                ui.checkbox(&mut self.crystal_viewer.show_unit_cell, "Show Unit Cell");
                ui.checkbox(&mut self.crystal_viewer.show_bonds, "Show Bonds");
            });
        });
    }

    fn render_chat_panel(&mut self, ui: &mut egui::Ui) {
        ui.heading("Chat");
        
        egui::ScrollArea::vertical()
            .max_height(ui.available_height() - 60.0)
            .show(ui, |ui| {
                for msg in &self.chat_history {
                    if msg.role == "system" {
                        continue;
                    }
                    
                    let (prefix, color) = match msg.role.as_str() {
                        "user" => ("You: ", egui::Color32::from_rgb(100, 149, 237)),
                        "assistant" => ("AI: ", egui::Color32::from_rgb(50, 205, 50)),
                        _ => ("", egui::Color32::GRAY),
                    };
                    
                    ui.horizontal_wrapped(|ui| {
                        ui.label(egui::RichText::new(prefix).color(color).strong());
                        ui.label(&msg.content);
                    });
                    ui.add_space(4.0);
                }
            });
        
        ui.separator();
        
        ui.horizontal(|ui| {
            let response = ui.add(
                egui::TextEdit::singleline(&mut self.chat_input)
                    .hint_text("Type a message...")
                    .desired_width(ui.available_width() - 60.0)
            );
            
            if response.lost_focus() && ui.input(|i| i.key_pressed(egui::Key::Enter)) {
                self.send_chat_message();
            }
            
            if ui.button("Send").clicked() {
                self.send_chat_message();
            }
        });
    }

    fn render_tools_panel(&mut self, ui: &mut egui::Ui) {
        ui.heading("Tools");
        
        // Connection controls
        ui.horizontal(|ui| {
            if self.mcp_connected {
                ui.label(format!("✓ {} tools loaded", self.available_tools.len()));
                if ui.button("Disconnect").clicked() {
                    self.disconnect_mcp();
                }
            } else {
                ui.label("Not connected");
                if ui.button("Connect MCP").clicked() {
                    self.connect_mcp();
                }
            }
        });
        
        ui.separator();
        
        if !self.mcp_connected || self.available_tools.is_empty() {
            ui.label("Connect to MCP server to see tools.");
            
            // Still show quick actions even when not connected
            ui.add_space(16.0);
            ui.heading("Quick Actions (requires MCP)");
            
            if ui.button("Generate Si (Diamond)").clicked() {
                self.selected_tool = Some("generate_crystal".to_string());
                self.tool_params = r#"{
  "composition": ["Si", "Si", "Si", "Si", "Si", "Si", "Si", "Si"],
  "space_group": 227
}"#.to_string();
                if !self.mcp_connected {
                    self.connect_mcp();
                }
            }
            return;
        }
        
        // Tool selection dropdown
        ui.label("Select Tool:");
        egui::ComboBox::from_id_salt("tool_selector")
            .selected_text(self.selected_tool.as_deref().unwrap_or("-- Select --"))
            .width(250.0)
            .show_ui(ui, |ui| {
                ui.selectable_value(&mut self.selected_tool, None, "-- Select --");
                for tool in &self.available_tools {
                    ui.selectable_value(&mut self.selected_tool, Some(tool.name.clone()), &tool.name);
                }
            });
        
        ui.add_space(8.0);
        ui.label("Parameters (JSON):");
        ui.add(
            egui::TextEdit::multiline(&mut self.tool_params)
                .desired_rows(6)
                .code_editor()
        );
        
        ui.add_space(8.0);
        
        let mut should_clear = false;
        
        ui.horizontal(|ui| {
            if ui.button("Clear").clicked() {
                should_clear = true;
            }
        });
        
        if should_clear {
            self.tool_params.clear();
            self.selected_tool = None;
        }
    }

    fn render_viewer_panel(&mut self, ui: &mut egui::Ui) {
        ui.horizontal(|ui| {
            ui.heading("Crystal Viewer");
            ui.with_layout(egui::Layout::right_to_left(egui::Align::Center), |ui| {
                if self.current_structure_json.is_some() {
                    if ui.button("Show Structure").clicked() {
                        self.generate_and_show_visualization();
                    }
                }
            });
        });
        
        if let Some(info) = self.crystal_viewer.get_structure_info() {
            ui.label(info);
            ui.separator();
        }
        
        let available_size = ui.available_size();
        let (response, painter) = ui.allocate_painter(available_size, egui::Sense::drag());
        
        if response.dragged() {
            let delta = response.drag_delta();
            self.crystal_viewer.rotate(delta.x, delta.y);
        }
        
        if response.hovered() {
            let scroll = ui.input(|i| i.raw_scroll_delta.y);
            if scroll > 0.0 {
                self.crystal_viewer.zoom_in();
            } else if scroll < 0.0 {
                self.crystal_viewer.zoom_out();
            }
        }
        
        let rect = response.rect;
        let center = rect.center();
        
        painter.rect_filled(rect, 0.0, egui::Color32::from_rgb(30, 30, 40));
        
        if self.crystal_viewer.structure.is_some() {
            let scale = 20.0 * self.crystal_viewer.zoom;
            let rot = self.crystal_viewer.rotation;
            
            let project = |pos: [f32; 3]| -> egui::Pos2 {
                let x = pos[0];
                let y = pos[1];
                let z = pos[2];
                
                let cos_y = rot[1].cos();
                let sin_y = rot[1].sin();
                let x1 = x * cos_y - z * sin_y;
                let z1 = x * sin_y + z * cos_y;
                
                let cos_x = rot[0].cos();
                let sin_x = rot[0].sin();
                let y1 = y * cos_x - z1 * sin_x;
                
                egui::Pos2::new(
                    center.x + x1 * scale,
                    center.y - y1 * scale,
                )
            };
            
            if self.crystal_viewer.show_unit_cell {
                for (p1, p2) in self.crystal_viewer.get_unit_cell_lines() {
                    let pos1 = project(p1);
                    let pos2 = project(p2);
                    painter.line_segment([pos1, pos2], egui::Stroke::new(1.0, egui::Color32::from_rgb(100, 100, 100)));
                }
            }
            
            for (pos, element, color, radius) in self.crystal_viewer.get_atom_positions() {
                let screen_pos = project(pos);
                let screen_radius = radius * scale * 0.5;
                let egui_color = egui::Color32::from_rgb(
                    (color[0] * 255.0) as u8,
                    (color[1] * 255.0) as u8,
                    (color[2] * 255.0) as u8,
                );
                painter.circle_filled(screen_pos, screen_radius, egui_color);
                painter.circle_stroke(screen_pos, screen_radius, egui::Stroke::new(1.0, egui::Color32::BLACK));
            }
        } else {
            painter.text(
                center,
                egui::Align2::CENTER_CENTER,
                "No structure loaded.\nGenerate or load a crystal structure.",
                egui::FontId::proportional(16.0),
                egui::Color32::GRAY,
            );
        }
        
        ui.horizontal(|ui| {
            ui.label("Atom Scale:");
            ui.add(egui::Slider::new(&mut self.crystal_viewer.atom_scale, 0.1..=1.0));
        });
    }

    fn render_settings_window(&mut self, ctx: &egui::Context) {
        let mut show = self.show_settings;
        let mut should_close = false;
        
        egui::Window::new("Settings")
            .open(&mut show)
            .resizable(false)
            .show(ctx, |ui| {
                ui.heading("MCP Server");
                ui.horizontal(|ui| {
                    ui.label("Server Path:");
                    ui.text_edit_singleline(&mut self.mcp_server_path);
                });
                
                ui.add_space(16.0);
                ui.heading("Ollama");
                ui.horizontal(|ui| {
                    ui.label("URL:");
                    ui.text_edit_singleline(&mut self.ollama_url);
                });
                ui.horizontal(|ui| {
                    ui.label("Model:");
                    ui.text_edit_singleline(&mut self.ollama_model);
                });
                
                ui.add_space(16.0);
                if ui.button("Apply").clicked() {
                    should_close = true;
                }
            });
        
        self.show_settings = show && !should_close;
    }
}

impl eframe::App for CrystalApp {
    fn update(&mut self, ctx: &egui::Context, _frame: &mut eframe::Frame) {
        egui::TopBottomPanel::top("menu_bar").show(ctx, |ui| {
            self.render_menu_bar(ui);
        });
        
        egui::TopBottomPanel::bottom("status_bar").show(ctx, |ui| {
            ui.horizontal(|ui| {
                let mcp_status = if self.mcp_connected { "🟢 MCP" } else { "🔴 MCP" };
                let llm_status = if self.llm_connected { "🟢 LLM" } else { "🔴 LLM" };
                ui.label(mcp_status);
                ui.label(llm_status);
                ui.separator();
                ui.label(&self.status_message);
            });
        });
        
        egui::SidePanel::left("chat_panel")
            .default_width(350.0)
            .show(ctx, |ui| {
                self.render_chat_panel(ui);
            });
        
        egui::SidePanel::right("tools_panel")
            .default_width(300.0)
            .show(ctx, |ui| {
                self.render_tools_panel(ui);
            });
        
        egui::CentralPanel::default().show(ctx, |ui| {
            self.render_viewer_panel(ui);
        });
        
        if self.show_settings {
            self.render_settings_window(ctx);
        }
    }
}
