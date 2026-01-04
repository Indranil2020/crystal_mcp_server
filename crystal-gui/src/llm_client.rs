//! LLM Client for communicating with Ollama
//! 
//! Provides chat interface to local LLM models via Ollama API.

use serde::{Deserialize, Serialize};
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
}

#[derive(Debug, Clone, Serialize)]
struct OllamaChatRequest {
    model: String,
    messages: Vec<ChatMessage>,
    stream: bool,
}

#[derive(Debug, Clone, Deserialize)]
struct OllamaChatResponse {
    message: ChatMessage,
    done: bool,
}

#[derive(Debug, Clone, Deserialize)]
struct OllamaTagsResponse {
    models: Vec<OllamaModel>,
}

#[derive(Debug, Clone, Deserialize)]
struct OllamaModel {
    name: String,
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

    pub fn chat(&self, messages: &[ChatMessage]) -> Result<String, LlmError> {
        let url = format!("{}/api/chat", self.base_url);
        
        let request = OllamaChatRequest {
            model: self.model.clone(),
            messages: messages.to_vec(),
            stream: false,
        };

        let response = self.client
            .post(&url)
            .json(&request)
            .send()?;

        if response.status().is_success() {
            let chat_response: OllamaChatResponse = response.json()?;
            Ok(chat_response.message.content)
        } else {
            Err(LlmError::ModelNotFound(self.model.clone()))
        }
    }

    pub fn chat_with_tools(&self, messages: &[ChatMessage], tools_description: &str) -> Result<String, LlmError> {
        let mut enhanced_messages = messages.to_vec();
        
        let tool_instruction = format!(
            r#"You are a crystal structure assistant. Output ONLY a single valid JSON object.

AVAILABLE TOOLS:
{}

TOOL SELECTION:
- Bulk crystals (Si, Ge, NaCl, GaAs, perovskites) -> generate_prototype
- Slabs/surfaces -> generate_slab
- Graphene, nanotubes -> generate_nanostructure
- MXene, hBN, MoS2 -> generate_2d_material

LATTICE CONSTANTS (Angstroms):
Si=5.43, Ge=5.66, NaCl=5.64, GaAs=5.65, Fe=2.87, Cu=3.61

EXAMPLES:

Ge slab (100) 4 layers:
{{"tool": "generate_slab", "params": {{"structure": {{"prototype": "diamond", "elements": {{"A": "Ge"}}, "lattice_constant": 5.66}}, "miller_indices": [1, 0, 0], "thickness": 4, "vacuum": 15.0, "symmetric": false}}}}

Si diamond:
{{"tool": "generate_prototype", "params": {{"prototype": "diamond", "elements": {{"A": "Si"}}, "lattice_constant": 5.43}}}}

NaCl rocksalt:
{{"tool": "generate_prototype", "params": {{"prototype": "rocksalt", "elements": {{"A": "Na", "B": "Cl"}}, "lattice_constant": 5.64}}}}

Graphene:
{{"tool": "generate_nanostructure", "params": {{"type": "graphene", "params": {{"size": [4, 4, 1], "vacuum": 15.0}}}}}}

Output ONLY the JSON. No text before or after."#,
            tools_description
        );
        
        if let Some(first) = enhanced_messages.first_mut() {
            if first.role == "system" {
                first.content = tool_instruction;
            }
        } else {
            enhanced_messages.insert(0, ChatMessage {
                role: "system".to_string(),
                content: tool_instruction,
            });
        }

        self.chat(&enhanced_messages)
    }
}

impl Default for LlmClient {
    fn default() -> Self {
        Self::new("http://localhost:11434", "llama3.2:1b")
    }
}
