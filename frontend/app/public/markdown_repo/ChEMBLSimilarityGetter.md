# ChEMBL Similarity Getter Documentation
        docs: <DocRenderer filePath="/markdown_repo/UniProtExplorer.md"/>,
---

## Overview
---

The **ChEMBL Similarity Getter** is a powerful molecular similarity search tool that finds structurally similar molecules in the ChEMBL database using Tanimoto similarity calculations. This tool allows researchers to discover molecules with similar chemical structures by starting with either a SMILES string or a ChEMBL ID, making it invaluable for drug discovery, chemical space exploration, and structure-activity relationship studies.

## Features
---
- **Dual Search Methods**: Find similar molecules using SMILES strings or ChEMBL IDs
- **Adjustable Similarity Threshold**: Control search sensitivity from 0-100% similarity
- **Tanimoto Similarity**: Uses industry-standard molecular fingerprint comparison
- **Real-time Validation**: Automatic SMILES and ChEMBL ID format checking
- **Interactive Results**: Browse and compare multiple similar molecules
- **Comprehensive Data**: Detailed molecular information for each similar compound
- **Visual Structure Display**: 2D molecular structure visualization and analysis

### Tool Sections
---
1. **Search Configuration**
   - Search type selector (SMILES or ChEMBL ID)
   - Similarity threshold input box (0-100%)
   - Input validation with real-time feedback
   - Advanced search options

2. **Input Interface**
   - Dynamic input fields based on search type
   - SMILES string validation using backend services
   - ChEMBL ID format validation
   - Clear error messaging and guidance

3. **Results Display**
   - List of similar molecules with similarity scores
   - Interactive molecule selection
   - Detailed molecular information panels
   - Structure visualization and analysis tools

## Data Source
---
The tool accesses the ChEMBL database similarity service, providing:

- **Molecular Fingerprints**: Chemical structure representations for similarity calculations
- **Tanimoto Coefficients**: Quantitative similarity measurements between molecules
- **Structural Data**: SMILES, molecular properties, and structural information
- **Cross-References**: Links to ChEMBL entries and related database records
- **Similarity Scores**: Precise percentage similarity measurements
- **Comprehensive Profiles**: Complete molecular data for similar compounds

## Search Types
---

### **SMILES-Based Similarity Search**
Find molecules similar to a compound represented by its SMILES string.

**What you'll enter:**
- **SMILES String**: Chemical structure in SMILES notation (e.g., "CCO", "C1=CC=CC=C1")
- **Similarity Threshold**: Minimum similarity percentage (default 70%)

**What you'll get:**
- List of molecules with similar chemical structures
- Similarity scores for each match (Tanimoto coefficient)
- Complete molecular profiles for similar compounds
- Structure visualizations and property comparisons

### **ChEMBL ID-Based Similarity Search**
Find molecules similar to a known ChEMBL compound using its database identifier.

**What you'll enter:**
- **ChEMBL ID**: Unique ChEMBL identifier (e.g., "CHEMBL25", "CHEMBL192")
- **Similarity Threshold**: Minimum similarity percentage (default 70%)

**What you'll get:**
- Structurally similar molecules from the ChEMBL database
- Quantitative similarity measurements
- Comparative molecular data and properties
- Access to related compounds and analogs

## Usage
---
1. **Select Search Method**: Choose between SMILES string or ChEMBL ID search
2. **Set Similarity Threshold**: Enter the minimum similarity percentage in the input box (typically 70-90%)
3. **Enter Search Input**:
   - For SMILES: Input a valid chemical structure notation
   - For ChEMBL ID: Enter the complete ChEMBL identifier
4. **Validate Input**: The tool automatically validates your input format
5. **Execute Search**: Submit your query to find similar molecules
6. **Browse Results**: 
   - Review the list of similar molecules with similarity scores
   - Select individual molecules to view detailed information
   - Use arrow keys for quick navigation
7. **Analyze Data**: Examine molecular structures, properties, and relationships

## What You'll See
---
When you use the ChEMBL Similarity Getter, the interface displays several key sections:

### Search Interface
- **Search Type Dropdown**: Selection between SMILES and ChEMBL ID search methods
- **Similarity Threshold Input**: Numeric input box with percentage indicator
- **Input Field**: Dynamic form that changes based on search type selection
- **Validation Indicators**: Real-time feedback with success/error icons
- **Submit Button**: Enabled only when input is valid

### Results Display
The similarity search results show comprehensive information:

**Similarity List:**
- Ranked list of similar molecules by similarity score
- ChEMBL IDs and preferred names where available
- Percentage similarity scores (Tanimoto coefficients)
- Quick selection interface for detailed viewing

**Selected Molecule Details:**
- Molecular structure visualization (2D viewer)
- Chemical property information panel
- Similarity score highlighting
- Navigation controls for browsing results

**Comprehensive Data:**
- Complete ChEMBL molecular profiles
- Cross-references to related databases
- Structural alerts and property predictions
- Download options for further analysis

## Example Searches
---

### Find Molecules Similar to Ethanol
```
Search Type: SMILES String
SMILES: CCO
Similarity Threshold: 70%
```
*Returns alcohols and other molecules with similar functional groups*

### Find Aspirin Analogs
```
Search Type: ChEMBL ID
ChEMBL ID: CHEMBL25
Similarity Threshold: 80%
```
*Discovers aspirin-like compounds and structural analogs*

### High-Similarity Drug Discovery
```
Search Type: SMILES String
SMILES: CN1C=NC2=C1C(=O)N(C(=O)N2C)C
Similarity Threshold: 90%
```
*Finds very similar compounds to caffeine for medicinal chemistry*

### Broad Chemical Space Exploration
```
Search Type: ChEMBL ID
ChEMBL ID: CHEMBL192
Similarity Threshold: 60%
```
*Explores diverse compounds related to sildenafil with lower similarity*

## Data Fields Description
---
The ChEMBL Similarity Getter displays rich similarity and molecular information:

### Similarity Metrics
- **Similarity Score**: Tanimoto coefficient as percentage (0-100%)
- **Ranking**: Molecules ordered by decreasing similarity
- **Threshold Filtering**: Only molecules above similarity threshold shown
- **Score Precision**: High-precision similarity calculations

### Molecular Information
- **ChEMBL ID**: Unique database identifier for each similar molecule
- **Preferred Name**: Common or systematic name where available
- **Molecular Structure**: SMILES representation and 2D visualization
- **Chemical Properties**: Molecular weight, formula, and descriptors

### Structure Analysis
- **2D Visualization**: Interactive molecular structure display
- **Functional Groups**: Identification of key chemical features
- **Structural Alerts**: Potential reactivity or toxicity warnings
- **Property Predictions**: Calculated molecular descriptors

### Cross-References
- **Database Links**: Connections to PubChem, PDB, and other databases
- **Literature**: Associated publications and research
- **Activity Data**: Biological activity information where available

## Advanced Features
---

### Similarity Threshold Control
- **Adjustable Sensitivity**: Fine-tune search specificity from 0-100% using input box
- **Real-time Updates**: Immediate feedback on threshold changes
- **Optimal Ranges**: Guidance on typical threshold values for different applications
- **Score Distribution**: Understanding of similarity score ranges

### Input Validation
- **SMILES Validation**: Backend verification of chemical structure validity
- **ChEMBL ID Checking**: Format validation for database identifiers
- **Error Prevention**: Clear feedback prevents invalid searches
- **Format Guidance**: Examples and hints for proper input formatting

### Results Navigation
- **Keyboard Controls**: Arrow key navigation through similar molecules
- **Quick Selection**: Click-based molecule browsing
- **Score Sorting**: Results automatically ranked by similarity
- **Batch Analysis**: Efficient handling of multiple similar compounds

### Data Export and Integration
- **Structured Results**: Complete molecular data in organized format
- **Similarity Matrices**: Quantitative comparison data
- **Cross-Database Links**: Direct access to external molecular databases
- **Research Integration**: Compatible with computational chemistry workflows

## Use Cases
---
This tool is particularly valuable for:
- **Drug Discovery**: Finding lead compounds and structural analogs
- **Medicinal Chemistry**: Exploring structure-activity relationships
- **Chemical Space Analysis**: Mapping molecular similarity networks
- **Scaffold Hopping**: Discovering alternative chemical frameworks
- **Patent Research**: Finding structurally related compounds
- **Library Design**: Building focused compound collections
- **Academic Research**: Educational exploration of chemical similarity
- **Competitive Intelligence**: Analyzing competitor compound portfolios
- **Toxicity Prediction**: Finding compounds with known safety profiles
- **Synthesis Planning**: Identifying similar compounds with known synthetic routes
- **Bioactivity Prediction**: Leveraging similarity for activity forecasting
- **Chemical Database Mining**: Systematic exploration of molecular diversity
