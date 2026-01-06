//! LLM Client for communicating with Ollama
//!
//! Provides chat interface to local LLM models via Ollama API with native tool calling.

use serde::{Deserialize, Serialize};
use serde_json::Value;
use thiserror::Error;

#[derive(Error, Debug)]
pub enum LlmError {
    #[error("HTTP request failed: {0}")]
    RequestError(#[from] reqwest::Error),
    #[error("Ollama not available at {0}")]
    NotAvailable(String),
    #[error("Model not found: {0}")]
    ModelNotFound(String),
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ChatMessage {
    pub role: String,
    pub content: String,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub tool_calls: Option<Vec<ToolCall>>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ToolCall {
    pub function: ToolFunction,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ToolFunction {
    pub name: String,
    pub arguments: Value,
}

#[derive(Debug, Clone, Serialize)]
struct OllamaToolDef {
    #[serde(rename = "type")]
    tool_type: String,
    function: OllamaFunctionDef,
}

#[derive(Debug, Clone, Serialize)]
struct OllamaFunctionDef {
    name: String,
    description: String,
    parameters: Value,
}

#[derive(Debug, Clone, Serialize)]
struct OllamaChatRequest {
    model: String,
    messages: Vec<ChatMessage>,
    stream: bool,
    #[serde(skip_serializing_if = "Option::is_none")]
    tools: Option<Vec<OllamaToolDef>>,
    #[serde(skip_serializing_if = "Option::is_none")]
    format: Option<String>,
}

#[derive(Debug, Clone, Deserialize)]
struct OllamaChatResponse {
    message: ResponseMessage,
}

#[derive(Debug, Clone, Deserialize)]
struct ResponseMessage {
    #[serde(default)]
    content: String,
    #[serde(default)]
    tool_calls: Option<Vec<ToolCall>>,
}

#[derive(Debug, Clone, Deserialize)]
struct OllamaTagsResponse {
    models: Vec<OllamaModel>,
}

#[derive(Debug, Clone, Deserialize)]
struct OllamaModel {
    name: String,
}

/// Result from chat_with_tools - either tool calls or plain text
pub enum ChatResult {
    ToolCalls(Vec<ToolCall>),
    Text(String),
}

pub struct LlmClient {
    base_url: String,
    model: String,
    client: reqwest::blocking::Client,
    available_models: Vec<String>,
}

impl LlmClient {
    pub fn new(base_url: &str, model: &str) -> Self {
        Self {
            base_url: base_url.to_string(),
            model: model.to_string(),
            client: reqwest::blocking::Client::new(),
            available_models: Vec::new(),
        }
    }

    pub fn check_connection(&mut self) -> Result<bool, LlmError> {
        let url = format!("{}/api/tags", self.base_url);
        let response = self.client.get(&url).send()?;

        if response.status().is_success() {
            let tags: OllamaTagsResponse = response.json()?;
            self.available_models = tags.models.into_iter().map(|m| m.name).collect();
            Ok(true)
        } else {
            Err(LlmError::NotAvailable(self.base_url.clone()))
        }
    }

    pub fn get_available_models(&self) -> &[String] {
        &self.available_models
    }

    pub fn set_model(&mut self, model: &str) {
        self.model = model.to_string();
    }

    pub fn get_model(&self) -> &str {
        &self.model
    }

    /// Chat with native Ollama tool calling
    /// 
    /// tools_json: Array of MCP tool definitions (from get_tools)
    /// Returns either tool calls from the LLM or plain text
    pub fn chat_with_tools(
        &self,
        messages: &[ChatMessage],
        tools_json: &[Value],
    ) -> Result<ChatResult, LlmError> {
        let url = format!("{}/api/chat", self.base_url);
        
        // Convert MCP tools to Ollama format - NO FILTERING, send all tools
        let ollama_tools: Vec<OllamaToolDef> = tools_json
            .iter()
            .filter_map(|tool| {
                let name = tool.get("name")?.as_str()?.to_string();
                let description = tool
                    .get("description")
                    .and_then(|d| d.as_str())
                    .unwrap_or("")
                    .to_string();
                let parameters = tool
                    .get("inputSchema")
                    .cloned()
                    .unwrap_or(serde_json::json!({"type": "object", "properties": {}}));

                Some(OllamaToolDef {
                    tool_type: "function".to_string(),
                    function: OllamaFunctionDef {
                        name,
                        description,
                        parameters,
                    },
                })
            })
            .collect();

        // Use messages as-is - no hardcoded system prompts
        // Tool schemas are sent via the tools parameter (standard MCP/OpenAI format)
        let request = OllamaChatRequest {
            model: self.model.clone(),
            messages: messages.to_vec(),
            stream: false,
            tools: if ollama_tools.is_empty() {
                None
            } else {
                Some(ollama_tools)
            },
            // DO NOT use format: "json" - it DISABLES native tool calling!
            // Native tool calling returns tool_calls array, which is the proper MCP-compatible approach
            format: None,
        };
        
        // DEBUG: Print exact LLM request including tool schemas
        eprintln!("\n[DEBUG] ========== LLM REQUEST ==========");
        eprintln!("[DEBUG] Model: {}", request.model);
        eprintln!("[DEBUG] Messages ({}):", request.messages.len());
        for (i, msg) in request.messages.iter().enumerate() {
            eprintln!("[DEBUG]   [{}] {}: {}", i, msg.role, &msg.content[..msg.content.len().min(200)]);
        }
        if let Some(ref tools) = request.tools {
            eprintln!("[DEBUG] Tools ({}) with FULL SCHEMAS:", tools.len());
            for t in tools {
                eprintln!("[DEBUG]   - {} ({})", t.function.name, t.function.description.chars().take(60).collect::<String>());
                // Show schema parameters
                if let Some(props) = t.function.parameters.get("properties") {
                    let keys: Vec<_> = props.as_object().map(|o| o.keys().cloned().collect()).unwrap_or_default();
                    eprintln!("[DEBUG]     Parameters: {:?}", keys);
                }
            }
        }
        eprintln!("[DEBUG] Format: {:?} (None = native tool calling)", request.format);
        eprintln!("[DEBUG] =====================================\n");

        let response = self.client.post(&url).json(&request).send()?;

        if response.status().is_success() {
            let chat_response: OllamaChatResponse = response.json()?;

            // Check for tool calls first
            if let Some(tool_calls) = chat_response.message.tool_calls {
                if !tool_calls.is_empty() {
                    return Ok(ChatResult::ToolCalls(tool_calls));
                }
            }

            // Fall back to text content
            Ok(ChatResult::Text(chat_response.message.content))
        } else {
            Err(LlmError::ModelNotFound(self.model.clone()))
        }
    }

    /// Legacy method for plain chat without tools (for backward compatibility)
    pub fn chat(&self, messages: &[ChatMessage], format: Option<String>) -> Result<String, LlmError> {
        let url = format!("{}/api/chat", self.base_url);

        let request = OllamaChatRequest {
            model: self.model.clone(),
            messages: messages.to_vec(),
            stream: false,
            tools: None,
            format,
        };

        let response = self.client.post(&url).json(&request).send()?;

        if response.status().is_success() {
            let chat_response: OllamaChatResponse = response.json()?;
            Ok(chat_response.message.content)
        } else {
            Err(LlmError::ModelNotFound(self.model.clone()))
        }
    }
}

impl Default for LlmClient {
    fn default() -> Self {
        Self::new("http://localhost:11434", "qwen2.5:7b")
    }
}
