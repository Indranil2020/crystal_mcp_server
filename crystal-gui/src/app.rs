//! Main Application State and UI
//! 
//! Implements the egui-based GUI with chat interface and crystal visualization.

use eframe::egui;
use serde_json::{json, Value};
use std::sync::{Arc, Mutex};

use crate::crystal_viewer::{CrystalStructure, CrystalViewer};
use crate::llm_client::{ChatMessage, LlmClient};
use crate::mcp_client::McpClient;

pub struct CrystalApp {
    mcp_client: Arc<Mutex<Option<McpClient>>>,
    llm_client: Arc<Mutex<LlmClient>>,
    crystal_viewer: CrystalViewer,
    
    chat_history: Vec<ChatMessage>,
    chat_input: String,
    
    status_message: String,
    mcp_connected: bool,
    llm_connected: bool,
    
    available_tools: Vec<String>,
    selected_tool: Option<String>,
    tool_params: String,
    
    show_settings: bool,
    ollama_url: String,
    ollama_model: String,
    mcp_server_path: String,
}

impl CrystalApp {
    pub fn new(_cc: &eframe::CreationContext<'_>) -> Self {
        let mcp_server_path = std::env::current_dir()
            .map(|p| p.parent().unwrap_or(&p).join("dist/index.js").to_string_lossy().to_string())
            .unwrap_or_else(|_| "/home/niel/git/crystal-mcp-server/dist/index.js".to_string());

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
            ollama_model: "llama3.2:1b".to_string(),
            mcp_server_path,
        }
    }

    fn connect_mcp(&mut self) {
        let path = self.mcp_server_path.clone();
        let mut client = McpClient::new(path);
        
        match client.start() {
            Ok(()) => {
                self.available_tools = client.get_tools().iter().map(|t| t.name.clone()).collect();
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
            self.chat_history.push(ChatMessage {
                role: "assistant".to_string(),
                content: "LLM not connected. Please connect to Ollama first.".to_string(),
            });
            return;
        }

        let tools_desc = self.available_tools.join(", ");
        let client = self.llm_client.lock().unwrap();
        
        match client.chat_with_tools(&self.chat_history, &tools_desc) {
            Ok(response) => {
                self.chat_history.push(ChatMessage {
                    role: "assistant".to_string(),
                    content: response,
                });
            }
            Err(e) => {
                self.chat_history.push(ChatMessage {
                    role: "assistant".to_string(),
                    content: format!("Error: {}", e),
                });
            }
        }
    }

    fn call_tool(&mut self) {
        let Some(ref tool_name) = self.selected_tool else {
            self.status_message = "No tool selected.".to_string();
            return;
        };

        let params: Value = serde_json::from_str(&self.tool_params).unwrap_or(json!({}));
        
        let mut client_guard = self.mcp_client.lock().unwrap();
        let Some(ref mut client) = *client_guard else {
            self.status_message = "MCP not connected.".to_string();
            return;
        };

        match client.call_tool(tool_name, params) {
            Ok(result) => {
                self.status_message = format!("Tool {} executed.", tool_name);
                
                for content in &result.content {
                    if content.content_type == "text" {
                        if let Ok(structure) = serde_json::from_str::<CrystalStructure>(&content.text) {
                            self.crystal_viewer.set_structure(structure);
                            self.status_message = "Structure loaded into viewer.".to_string();
                        } else if content.text.contains("\"structure\"") {
                            if let Ok(wrapper) = serde_json::from_str::<Value>(&content.text) {
                                if let Some(structure_val) = wrapper.get("structure") {
                                    if let Ok(structure) = serde_json::from_value::<CrystalStructure>(structure_val.clone()) {
                                        self.crystal_viewer.set_structure(structure);
                                        self.status_message = "Structure loaded into viewer.".to_string();
                                    }
                                }
                            }
                        }
                    }
                }
            }
            Err(e) => {
                self.status_message = format!("Tool error: {}", e);
            }
        }
    }

    fn render_menu_bar(&mut self, ui: &mut egui::Ui) {
        egui::menu::bar(ui, |ui| {
            ui.menu_button("File", |ui| {
                if ui.button("Settings").clicked() {
                    self.show_settings = true;
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
        
        if !self.mcp_connected {
            ui.label("Connect to MCP server to see tools.");
            return;
        }
        
        egui::ComboBox::from_label("Select Tool")
            .selected_text(self.selected_tool.as_deref().unwrap_or("None"))
            .show_ui(ui, |ui| {
                for tool in &self.available_tools {
                    ui.selectable_value(&mut self.selected_tool, Some(tool.clone()), tool);
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
        if ui.button("Execute Tool").clicked() {
            self.call_tool();
        }
        
        ui.add_space(16.0);
        ui.heading("Quick Actions");
        
        if ui.button("Generate Si (Diamond)").clicked() {
            self.selected_tool = Some("generate_crystal".to_string());
            self.tool_params = r#"{
  "composition": ["Si", "Si", "Si", "Si", "Si", "Si", "Si", "Si"],
  "space_group": 227
}"#.to_string();
        }
        
        if ui.button("Generate NaCl (Rocksalt)").clicked() {
            self.selected_tool = Some("generate_prototype".to_string());
            self.tool_params = r#"{
  "prototype": "rocksalt",
  "elements": {"A": "Na", "X": "Cl"}
}"#.to_string();
        }
        
        if ui.button("Generate Perovskite").clicked() {
            self.selected_tool = Some("generate_prototype".to_string());
            self.tool_params = r#"{
  "prototype": "perovskite",
  "elements": {"A": "Ca", "B": "Ti", "X": "O"}
}"#.to_string();
        }
    }

    fn render_viewer_panel(&mut self, ui: &mut egui::Ui) {
        ui.heading("Crystal Viewer");
        
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
