
# Crystal MCP Server Documentation

Welcome to the comprehensive documentation for the Crystal MCP Server. This documentation suite features extensive **visual diagrams, flowcharts, and decision trees** designed for scientific precision and accuracy.

## 📚 Documentation Suite

### Getting Started

#### 🎯 [Quick Start & Testing Guide](./testing_guide.md)
**Your First Steps.**
- E2E testing with MCP Inspector
- Integration with Claude Desktop
- Automated test suite walkthrough
- Installation verification

### Core Documentation

#### 📖 [Comprehensive Generator Catalog](./comprehensive_generator_catalog.md)
**Complete Visual Reference for All 228 Operations.**
- **Visual Category Hierarchy**: All 20 categories with operation trees
- **Detailed Examples**: JSON request/response for every major operation
- **Structure Diagrams**: Visual representations of generated structures
- **Application Index**: Find operations by scientific domain
- **Operation Counts**: Category-wise breakdown of capabilities

#### 🏗️ [System Architecture](./architecture.md)
**Complete System Design with Visual Diagrams.**
- **Component Dependency Graph**: TypeScript ↔ Python interaction
- **Request Lifecycle**: Timing diagrams with performance metrics
- **Schema Transformation Flow**: Data conversion pipelines
- **Generator Category Map**: Visual organization of all generators
- **Error Handling Flow**: Comprehensive error propagation chains
- **Memory & Performance Characteristics**: Resource usage patterns
- **Security Model**: Sandboxing and validation architecture

#### 🔄 [Data Flow Guide](./dataflow.md)
**Visual Flow Diagrams for All Operations.**
- **2D TMD Generation Pipeline**: Complete flow from request to structure
- **Defect Generation with Supercell**: Multi-stage processing
- **MLFF Optimization Workflow**: Machine learning force field integration
- **Heterostructure Stacking**: Non-orthogonal cell handling
- **Workflow Editing Chain**: Iterative refinement flows
- **Export Format Conversions**: Multi-format output pipelines
- **Error Recovery & Retry Logic**: Fault-tolerant operation flows

### Usage & Workflows

#### ⚡ [Workflows & Usage Patterns](./workflows.md)
**Common Scientific Workflows with Complete Examples.**
- **Defect Study Workflow**: Vacancy formation energy calculations
- **Surface Catalysis**: Slab generation → adsorbate placement → optimization
- **2D Heterostructures**: Twisted bilayer construction (magic angles)
- **High-Throughput Screening**: Space group scanning with MLFF filtering
- **Battery Materials**: Cathode delithiation and diffusion barriers
- **Iterative Refinement**: Structure optimization chains with version control
- **Multi-Material Comparison**: Parallel workflow execution
- **Best Practices**: Validation checkpoints and progressive refinement

#### 🔬 [Crystallography Reference Guide](./crystallography_guide.md)
**Visual Structure Library & Crystallographic Concepts.**
- **Common Structure Types**: FCC, BCC, diamond, wurtzite, perovskite (with ASCII diagrams)
- **Symmetry Operations**: Rotations, mirrors, inversion (visualized)
- **Miller Indices**: Plane visualization for (100), (110), (111)
- **Bravais Lattices**: Complete hierarchy of 14 lattice types
- **Reciprocal Lattice**: Brillouin zones and k-point sampling
- **Coordination Polyhedra**: Tetrahedral, octahedral, cubic environments
- **Point Groups vs Space Groups**: Symmetry relationships
- **Practical Tips**: Structure generation best practices

### Troubleshooting & Reference

#### 🔧 [Troubleshooting Guide](./troubleshooting_guide.md)
**Systematic Debugging with Decision Trees.**
- **Quick Diagnostic Flowchart**: Identify issue category in seconds
- **Installation Issues**: Decision tree for dependency problems
- **Connection Problems**: MCP integration debugging
- **Generation Failures**: Structure generation error diagnosis
- **Scientific Accuracy Issues**: Validation and correctness checks
- **Performance Problems**: Optimization and profiling guidance
- **Export Errors**: Format conversion troubleshooting
- **Complete Error Reference**: All error codes with solutions

#### 📋 [API Reference](./api_reference.md)
**Auto-Generated Complete Operation Catalog.**
- All **228** available generator functions
- Organized by 20 categories (Bulk, Defect, Twist, Surface, etc.)
- Parameter signatures and type information
- Generated from source code for 100% accuracy

### Development & Extension

#### 👨‍💻 [Developer Guide](./developer_guide.md)
**For Contributors & Extenders.**
- **Development Setup**: Environment configuration and dependencies
- **Adding New Generators**: Step-by-step integration guide
- **Schema & Data Flow**: Understanding internal data structures
- **Testing Strategies**: Unit, integration, and E2E testing
- **Code Style & Best Practices**: Maintaining scientific accuracy
- **Debugging Tips**: Common issues and solutions
- **Contributing Guidelines**: Pull request workflow

#### 🛠️ [Server Capabilities & FAQ](./server_details.md)
**Technical Specifications.**
- Operations, Resources, State management
- I/O handling and performance characteristics
- Frequently asked questions

---

## 📊 Visual Documentation Philosophy

This documentation emphasizes **graphical representations** over pure text:
- **Mermaid Diagrams**: Flowcharts, sequence diagrams, state machines
- **ASCII Structure Diagrams**: Crystal structures and coordination geometries
- **Decision Trees**: Systematic troubleshooting paths
- **Visual Hierarchies**: Category and operation organization
- **Flow Diagrams**: Data pipelines and processing stages

Designed for **scientific precision and accuracy** - suitable for researchers and developers who demand exact specifications.

---

## 🚀 Quick Navigation

### By User Type
- **First-Time Users** → [Testing Guide](./testing_guide.md) → [Workflows](./workflows.md)
- **Scientists/Researchers** → [Generator Catalog](./comprehensive_generator_catalog.md) → [Crystallography Guide](./crystallography_guide.md)
- **Developers** → [Architecture](./architecture.md) → [Developer Guide](./developer_guide.md) → [Data Flow](./dataflow.md)
- **Troubleshooting** → [Troubleshooting Guide](./troubleshooting_guide.md)

### By Task
- **Generate Structures** → [Generator Catalog](./comprehensive_generator_catalog.md)
- **Understand System** → [Architecture](./architecture.md) + [Data Flow](./dataflow.md)
- **Run Workflows** → [Workflows](./workflows.md)
- **Fix Problems** → [Troubleshooting Guide](./troubleshooting_guide.md)
- **Extend Server** → [Developer Guide](./developer_guide.md)
- **Learn Crystallography** → [Crystallography Guide](./crystallography_guide.md)

---

## 📂 Quick Links
- [Source Code](../src/)
- [Generator Registry](../src/python/generators/__init__.py)
- [Test Suite](../tests/)
