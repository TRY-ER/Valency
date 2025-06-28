# RCSB Structure Similarity Search Documentation
---

## Overview
---

The **RCSB Structure Similarity Search** is a specialized protein structure analysis tool that finds protein structures with similar 3D shapes and conformations. This tool uses advanced structural comparison algorithms to identify proteins that share similar spatial arrangements, even when their sequences may differ significantly. It's particularly valuable for discovering functional relationships and structural homologs in the protein universe.

## Features
---
- **Multiple Input Methods**: Search using PDB entry IDs or custom structure files via URL
- **Shape Matching Algorithms**: Choose between strict and relaxed similarity matching criteria
- **Flexible Search Targets**: Compare against different structural representations
- **Multiple Return Types**: Get results as assemblies or individual polymer entities
- **File Format Support**: Works with various structural file formats (CIF, PDB, BCIF)
- **Real-time Results**: Interactive interface with immediate feedback
- **Similarity Scoring**: Results ranked by structural similarity

### Tool Sections
---
1. **Search Method Selector**
   - Dropdown to choose between PDB ID search or file URL search
   - Clear descriptions for each search type
   - Dynamic interface adaptation

2. **Input Configuration**
   - PDB entry ID input for known structures
   - File URL input for custom structures
   - Assembly ID specification for multi-chain complexes
   - File format selection for URL-based searches

3. **Search Parameters**
   - Similarity operator selection (strict vs relaxed matching)
   - Target search space configuration
   - Return type specification
   - Result limit controls

4. **Results Display**
   - Structured list of similar structures
   - Similarity rankings and scores
   - Expandable data viewer for detailed information

## Data Source
---
The tool searches the RCSB Protein Data Bank using structural comparison algorithms, providing:

- **3D Structure Comparisons**: Shape-based similarity matching
- **Assembly Information**: Complete protein complex structures
- **Polymer Entity Data**: Individual protein chain comparisons
- **Similarity Scores**: Quantitative measures of structural similarity
- **Cross-References**: Links to related structures and annotations
- **Geometric Analysis**: Spatial arrangement comparisons

## Search Types
---

### **PDB Entry ID Search**
Find structures similar to a known protein structure using its PDB identifier.

**What you'll enter:**
- **PDB Entry ID**: Four-character PDB code (e.g., "4HHB", "1TIM")
- **Assembly ID**: Specific assembly within the structure (default: "1")
- **Similarity Operator**: Matching strictness level
- **Target Search Space**: What to compare against
- **Return Type**: Format of results

**What you'll get:**
- List of structurally similar PDB entries
- Similarity rankings based on 3D shape matching
- Assembly or entity identifiers for further analysis

### **File URL Search**
Find structures similar to a custom structure file provided via web URL.

**What you'll enter:**
- **File URL**: Direct link to structure file (e.g., "https://files.rcsb.org/view/4HHB.cif")
- **File Format**: Format specification (CIF, PDB, BCIF, etc.)
- **Similarity Operator**: Matching algorithm selection
- **Target Search Space**: Comparison target type
- **Return Type**: Result identifier format

**What you'll get:**
- Structures similar to your uploaded file
- Ranked results based on structural similarity
- Cross-references to PDB database entries

## Search Parameters
---

### **Similarity Operators**
Choose how strictly structures should match:

**Strict Shape Match:**
- More stringent similarity requirements
- Finds very closely related structures
- Higher confidence in functional similarity
- Fewer but more relevant results

**Relaxed Shape Match:**
- More lenient similarity criteria
- Broader range of structural relatives
- Captures distant structural relationships
- More comprehensive result sets

### **Target Search Spaces**
Define what to compare your structure against:

**Assembly:**
- Compare against complete protein assemblies
- Includes multi-chain complexes
- Biological units and quaternary structures
- Full protein machinery comparisons

**Polymer Entity Instance:**
- Compare against individual protein chains
- Single polymer comparisons
- Domain-level structural analysis
- Chain-specific similarities

### **Return Types**
Specify the format of your results:

**Assembly Returns:**
- Complete assembly identifiers
- Multi-chain complex references
- Biological unit designations
- Quaternary structure identifiers

**Polymer Entity Returns:**
- Individual chain identifiers
- Single polymer references
- Domain-specific results
- Chain-level comparisons

## File Format Support
---
The tool accepts various structural file formats:

- **CIF**: Crystallographic Information File format
- **BCIF**: Binary CIF for faster processing
- **PDB**: Traditional Protein Data Bank format
- **CIF.GZ**: Compressed CIF files
- **PDB.GZ**: Compressed PDB files

## Usage
---
1. **Select Search Method**: Choose between PDB Entry ID or File URL search
2. **Enter Structure Information**:
   - For PDB ID: Enter the 4-character code and assembly ID
   - For File URL: Provide the complete URL and specify file format
3. **Configure Parameters**:
   - Select similarity operator (strict or relaxed)
   - Choose target search space (assembly or polymer entity)
   - Set return type preference
4. **Execute Search**: Click "Find Similar Structures" to start the analysis
5. **Review Results**: Explore the ranked list of similar structures
6. **Analyze Matches**: Use the data viewer to examine detailed similarity information

## What You'll See
---
When you use the Structure Similarity Search, the interface displays:

### Search Configuration
- **Method Selection**: Dropdown for PDB ID vs File URL search
- **Input Fields**: Dynamic forms based on selected search method
- **Parameter Controls**: Options for similarity matching and result formatting
- **Submit Button**: Execute button with loading feedback

### Results Display
The results show structures similar to your query:

**Structure Identifiers:**
- PDB entry codes for similar structures
- Assembly or entity IDs based on return type
- Cross-references to database entries

**Similarity Information:**
- Ranking based on structural similarity
- Quantitative similarity measures
- Geometric comparison data

**Detailed Metadata:**
- Complete structural information
- Experimental details for similar structures
- Functional annotations and classifications

## Example Searches
---

### Find Structures Similar to Hemoglobin
```
Search Type: PDB Entry ID Search
PDB Entry ID: 4HHB
Assembly ID: 1
Similarity Operator: Strict Shape Match
Target Search Space: Assembly
Return Type: Assembly
```
*Finds protein structures with similar overall shape to hemoglobin*

### Compare Custom Structure File
```
Search Type: File URL Search
File URL: https://files.rcsb.org/view/1TIM.cif
File Format: CIF
Similarity Operator: Relaxed Shape Match
Target Search Space: Polymer Entity Instance
Return Type: Polymer Entity
```
*Compares a custom structure file against individual protein chains*

### Broad Structural Survey
```
Search Type: PDB Entry ID Search
PDB Entry ID: 1TIM
Assembly ID: 1
Similarity Operator: Relaxed Shape Match
Target Search Space: Assembly
Return Type: Assembly
```
*Finds a wide range of structurally related proteins*

## Advanced Applications
---

### **Functional Annotation**
- Discover proteins with similar functions through structural similarity
- Identify functional domains and motifs
- Predict protein function based on structural homologs

### **Drug Discovery**
- Find alternative target proteins with similar binding sites
- Identify off-target effects through structural similarity
- Explore protein families for drug development

### **Evolutionary Studies**
- Trace structural evolution across protein families
- Identify conserved structural motifs
- Study divergent evolution of protein folds

### **Structural Classification**
- Classify novel structures into known fold families
- Identify structural relationships across species
- Organize protein structure space

## Use Cases
---
This tool is particularly valuable for:
- **Structural Biology Research**: Discovering related protein structures and folds
- **Protein Function Prediction**: Using structural similarity to infer function
- **Drug Target Analysis**: Finding structurally similar proteins for drug development
- **Evolutionary Studies**: Tracing structural relationships across species
- **Structural Classification**: Organizing and categorizing protein structures
- **Comparative Analysis**: Studying structural variations within protein families
- **Method Validation**: Comparing experimental and predicted structures
- **Database Mining**: Large-scale structural similarity searches
- **Educational Research**: Teaching structure-function relationships
- **Quality Assessment**: Evaluating structural models against known structures
