//! Crystal Structure Generator GUI
//! 
//! A Rust-based GUI application that provides:
//! - Chat interface for interacting with crystal generation
//! - MCP client for communicating with crystal-mcp-server
//! - 3D visualization of crystal structures
//! - Structure editing capabilities

mod mcp_client;
mod llm_client;
mod crystal_viewer;
mod structure_export;
mod app;

use app::CrystalApp;
use eframe::egui;

fn main() -> eframe::Result<()> {
    tracing_subscriber::fmt::init();
    
    let options = eframe::NativeOptions {
        viewport: egui::ViewportBuilder::default()
            .with_inner_size([1400.0, 900.0])
            .with_min_inner_size([800.0, 600.0])
            .with_title("Crystal Structure Generator"),
        ..Default::default()
    };
    
    eframe::run_native(
        "Crystal Structure Generator",
        options,
        Box::new(|cc| Ok(Box::new(CrystalApp::new(cc)))),
    )
}
