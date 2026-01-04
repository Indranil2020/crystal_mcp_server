//! MCP Client for communicating with crystal-mcp-server
//! 
//! Implements JSON-RPC over stdio to communicate with the MCP server.

use std::io::{BufRead, BufReader, Write};
use std::process::{Child, Command, Stdio};
use std::sync::atomic::{AtomicU64, Ordering};

use serde::{Deserialize, Serialize};
use serde_json::{json, Value};
use thiserror::Error;

#[derive(Error, Debug)]
pub enum McpError {
    #[error("Failed to spawn MCP server: {0}")]
    SpawnError(String),
    #[error("Failed to send request: {0}")]
    SendError(String),
    #[error("Failed to read response: {0}")]
    ReadError(String),
    #[error("JSON parse error: {0}")]
    JsonError(#[from] serde_json::Error),
    #[error("MCP error: {code} - {message}")]
    McpError { code: i64, message: String },
    #[error("Server not initialized")]
    NotInitialized,
    #[error("IO error: {0}")]
    IoError(#[from] std::io::Error),
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Tool {
    pub name: String,
    pub description: String,
    pub input_schema: Value,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ToolResult {
    pub content: Vec<ContentItem>,
    #[serde(default)]
    pub is_error: bool,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ContentItem {
    #[serde(rename = "type")]
    pub content_type: String,
    pub text: String,
}

pub struct McpClient {
    process: Option<Child>,
    request_id: AtomicU64,
    initialized: bool,
    tools: Vec<Tool>,
    server_path: String,
}

impl McpClient {
    pub fn new(server_path: String) -> Self {
        Self {
            process: None,
            request_id: AtomicU64::new(1),
            initialized: false,
            tools: Vec::new(),
            server_path,
        }
    }

    fn next_id(&self) -> u64 {
        self.request_id.fetch_add(1, Ordering::SeqCst)
    }

    pub fn start(&mut self) -> Result<(), McpError> {
        let child = Command::new("node")
            .arg(&self.server_path)
            .stdin(Stdio::piped())
            .stdout(Stdio::piped())
            .stderr(Stdio::null())
            .spawn()
            .map_err(|e| McpError::SpawnError(e.to_string()))?;

        self.process = Some(child);
        self.initialize()?;
        self.list_tools()?;
        Ok(())
    }

    fn send_request(&mut self, method: &str, params: Value) -> Result<Value, McpError> {
        let id = self.next_id();
        let request = json!({
            "jsonrpc": "2.0",
            "id": id,
            "method": method,
            "params": params
        });

        let process = self.process.as_mut().ok_or(McpError::NotInitialized)?;
        
        {
            let stdin = process.stdin.as_mut().ok_or(McpError::SendError("No stdin".to_string()))?;
            let request_str = serde_json::to_string(&request)?;
            writeln!(stdin, "{}", request_str).map_err(|e| McpError::SendError(e.to_string()))?;
            stdin.flush().map_err(|e| McpError::SendError(e.to_string()))?;
        }

        let stdout = process.stdout.as_mut().ok_or(McpError::ReadError("No stdout".to_string()))?;
        let mut reader = BufReader::new(stdout);
        let mut line = String::new();
        reader.read_line(&mut line).map_err(|e| McpError::ReadError(e.to_string()))?;

        let response: Value = serde_json::from_str(&line)?;
        
        if let Some(error) = response.get("error") {
            let code = error.get("code").and_then(|c| c.as_i64()).unwrap_or(-1);
            let message = error.get("message").and_then(|m| m.as_str()).unwrap_or("Unknown error").to_string();
            return Err(McpError::McpError { code, message });
        }

        Ok(response.get("result").cloned().unwrap_or(Value::Null))
    }

    fn send_notification(&mut self, method: &str) -> Result<(), McpError> {
        let process = self.process.as_mut().ok_or(McpError::NotInitialized)?;
        
        let notification = json!({
            "jsonrpc": "2.0",
            "method": method
        });

        let stdin = process.stdin.as_mut().ok_or(McpError::SendError("No stdin".to_string()))?;
        let notification_str = serde_json::to_string(&notification)?;
        writeln!(stdin, "{}", notification_str).map_err(|e| McpError::SendError(e.to_string()))?;
        stdin.flush().map_err(|e| McpError::SendError(e.to_string()))?;
        Ok(())
    }

    fn initialize(&mut self) -> Result<(), McpError> {
        let params = json!({
            "protocolVersion": "2024-11-05",
            "capabilities": {},
            "clientInfo": {
                "name": "crystal-gui",
                "version": "0.1.0"
            }
        });

        let _result = self.send_request("initialize", params)?;
        self.send_notification("notifications/initialized")?;
        self.initialized = true;
        Ok(())
    }

    fn list_tools(&mut self) -> Result<(), McpError> {
        let result = self.send_request("tools/list", json!({}))?;
        
        if let Some(tools) = result.get("tools").and_then(|t| t.as_array()) {
            self.tools = tools
                .iter()
                .filter_map(|t| serde_json::from_value(t.clone()).ok())
                .collect();
        }
        Ok(())
    }

    pub fn get_tools(&self) -> &[Tool] {
        &self.tools
    }

    pub fn call_tool(&mut self, name: &str, arguments: Value) -> Result<ToolResult, McpError> {
        let params = json!({
            "name": name,
            "arguments": arguments
        });

        let result = self.send_request("tools/call", params)?;
        let tool_result: ToolResult = serde_json::from_value(result)?;
        Ok(tool_result)
    }

    pub fn is_initialized(&self) -> bool {
        self.initialized
    }

    pub fn stop(&mut self) {
        if let Some(mut process) = self.process.take() {
            let _ = process.kill();
        }
        self.initialized = false;
    }
}

impl Drop for McpClient {
    fn drop(&mut self) {
        self.stop();
    }
}
