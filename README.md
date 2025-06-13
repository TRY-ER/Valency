# VALENCY
### Agentic Drug Discovery and Molecular Analysis Platform

<div align="center">
  <img src="./Resources/images/logo_rendered_main.png" alt="VALENCY Logo" width="200" style="margin: 20px 0;">
  
  [![License](https://img.shields.io/badge/license-MIT-blue.svg)](LICENSE)
  [![Python](https://img.shields.io/badge/python-3.12+-blue.svg)](https://python.org)
  [![React](https://img.shields.io/badge/react-18+-blue.svg)](https://reactjs.org)
  [![FastAPI](https://img.shields.io/badge/fastapi-latest-green.svg)](https://fastapi.tiangolo.com)
</div>

---

## 🧬 Overview

**VALENCY** is a comprehensive AI-powered platform for drug discovery and molecular analysis that combines cutting-edge artificial intelligence with robust scientific tools. Our platform offers both manual tool access and intelligent agentic automation, making complex molecular research accessible to scientists at all levels.

### 🎯 Mission
To accelerate drug discovery and molecular research through intelligent automation, comprehensive analysis tools, and seamless integration of multiple scientific databases and algorithms.

---

## ✨ Key Features

### 🤖 **Agentic Automation System**
- **Master Agent**: Natural language query processing with intelligent task delegation
- **Specialized Sub-Agents**: Domain-specific agents for proteins, drugs, and optimization
- **Multi-Agent Coordination**: Seamless orchestration of complex workflows
- **Conversational Interface**: Interactive chat-based research assistance

### 🔬 **Comprehensive Analysis Tools**
- **Molecular Structure Analysis**: SMILES-based exploration and visualization
- **Protein Structure Exploration**: UniProt & PDB integration with 3D visualization
- **Polymer Characterization**: Advanced polymer structure analysis
- **ADMET Property Prediction**: AI-powered drug safety and efficacy assessment

### 🎯 **Advanced Identification & Search**
- **Multi-Database Integration**: ChEMBL (1.6M+ compounds), PubChem (100M+ compounds), UniProt (150M+ entries)
- **Similarity Search Algorithms**: Chemical fingerprint comparison and substructure matching
- **Cross-Database Identification**: Unified search across multiple scientific databases

### 🧪 **AI-Powered Drug Optimization**
- **BRICS-Based Generation**: Fragment-based molecular generation
- **LSTM Neural Networks**: Advanced molecular candidate generation
- **Structure-Based Drug Design**: Intelligent lead compound optimization
- **Real-Time Property Prediction**: Instant ADMET analysis

### 📊 **Activity Management**
- **Comprehensive Tracking**: Monitor all research activities and agent interactions
- **Persistent Logs**: Fail-safe activity storage across browser sessions
- **Usage Analytics**: Tool usage statistics and resource optimization

---

## 🏗️ Architecture & Components

### **Frontend** (`/frontend/`)
Modern React-based user interface featuring:
- Interactive molecular visualizations
- Real-time chat interface with AI agents
- Responsive design with advanced animations
- Tool integration and workflow management
- 3D molecular structure viewers

### **Backend** (`/backend/`)
FastAPI-powered server infrastructure:
- RESTful API endpoints for all platform features
- Authentication and user management
- Database integration and data processing
- Tool orchestration and workflow management

### **Agents** (`/agents/`)
Intelligent agent system with specialized capabilities:
- **Master Agent**: Central orchestrator and query processor
- **Drug Agent**: Database searches and compound exploration
- **Protein Agent**: Protein analysis and structure prediction
- **Drug Optimization Agent**: ADMET prediction and molecular generation

### **Authentication Service** (`/auth_service/`)
Secure user authentication and authorization:
- JWT-based authentication
- User registration and management
- Email verification system
- Role-based access control

### **Engine** (`/engine/`)
Core computational algorithms:
- Vector storage and similarity search
- Molecular generation algorithms
- ADMET property prediction models
- Database integration utilities

### **Model Context Protocol Servers** (`/MCPs/`)
Specialized Model Context Protocol(MCP) servers for:
- ADMET property prediction
- PubChem Query
- ChEMBL Query
- Alphafold Data Retrieval
- UniProt Query
- RCSB PDB Query
- BRICS Generation Requests

---

## 🔧 Supported Databases

### **Scientific Databases**
| Database | Content | Records |
|----------|---------|---------|
| **ChEMBL** | Bioactivity data | 1.6M+ compounds |
| **PubChem** | Chemical information | 100M+ compounds |
| **UniProt** | Protein sequences | 150M+ entries |
| **PDB** | Protein structures | 200K+ structures |

---

## 🚀 Getting Started

### **Prerequisites**
- Python 3.12+
- Node.js 18+
- Modern web browser with WebGL support

### **Quick Overview**
1. **Explore the Platform**: Navigate through our interactive homepage
2. **Try Manual Tools**: Use individual analysis tools for specific tasks
3. **Engage with AI Agents**: Ask natural language questions for complex workflows
4. **Manage Activities**: Track your research progress and results

---

## 📱 User Interface

### **Homepage Features**
- **Interactive Platform Overview**: Comprehensive introduction to all features
- **Quick Action Buttons**: Direct access to popular tools
- **Animated Infographics**: Visual representation of platform capabilities
- **Progress Tracking**: Real-time system metrics and performance indicators

### **Tool Categories**
1. **Structure Analysis**: Molecular, protein, and polymer exploration
2. **Identification**: Database searches and similarity analysis
3. **Optimization**: AI-powered drug design and generation

---

## 🎯 Use Cases

### **Academic Research**
- Molecular property analysis for research publications
- Protein structure prediction and analysis
- Educational molecular visualization
- Cross-database literature searches

### **Pharmaceutical Development**
- Lead compound identification and optimization
- ADMET property screening
- Drug-target interaction analysis
- Candidate molecule generation

### **Biotechnology Applications**
- Protein engineering and design
- Molecular similarity screening
- Bioactivity understanding 

---