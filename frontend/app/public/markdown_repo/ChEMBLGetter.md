# ChEMBL Getter Documentation
---

## Overview
---

The **ChEMBL Getter** is a comprehensive molecular database search tool that provides access to the ChEMBL database, one of the world's largest repositories of bioactive molecules and their biological activities. This tool allows users to search for and retrieve detailed information about small molecules, drugs, and bioactive compounds using various search methods including ChEMBL IDs, preferred names, and synonyms.

## Features
---
- **Multiple Search Methods**: Four different ways to find molecules in the ChEMBL database
- **ChEMBL ID Search**: Direct lookup using unique ChEMBL identifiers
- **Multiple ChEMBL IDs**: Batch processing of multiple molecules simultaneously
- **Preferred Name Search**: Find molecules by their common scientific names
- **Synonym Search**: Discover molecules using brand names, trade names, or alternative names
- **Real-time Validation**: Automatic ChEMBL ID format checking
- **Comprehensive Results**: Detailed molecular information including structures, properties, and cross-references
- **Interactive Interface**: Dynamic forms that adapt based on search type selection

### Tool Sections
---
1. **Search Type Selector**
   - Dropdown menu to choose search method
   - Clear descriptions for each search option
   - Dynamic interface adaptation

2. **Search Input Forms**
   - Text fields with validation for ChEMBL IDs
   - Input areas for molecular names and synonyms
   - Multiple ID management for batch searches
   - Real-time format validation feedback

3. **Results Display**
   - Detailed molecular information panels
   - Multiple candidate selection for ambiguous searches
   - Batch result navigation for multiple molecules
   - Comprehensive data viewer with expandable sections

## Data Source
---
The tool accesses the ChEMBL database, providing comprehensive information about:

- **Molecular Structures**: Chemical structures, SMILES, InChI, and molecular properties
- **Drug Information**: Approved drugs, development phases, and regulatory information
- **Biological Activities**: Bioassay data, target interactions, and activity measurements
- **Cross-References**: Links to other databases (PubChem, PDB, UniProt, etc.)
- **Synonyms and Names**: Trade names, brand names, and alternative identifiers
- **Molecular Properties**: Physical and chemical properties, drug-likeness metrics

## Search Types
---

### **ChEMBL ID Search**
Direct lookup of molecules using their unique ChEMBL database identifiers.

**What you'll enter:**
- **ChEMBL ID**: Unique identifier starting with "CHEMBL" followed by numbers (e.g., "CHEMBL25", "CHEMBL192")

**What you'll get:**
- Complete molecular profile for the specific compound
- Chemical structure information (SMILES, InChI, molfile)
- Molecular properties and drug-likeness metrics
- Cross-references to other databases
- Regulatory and development information

### **Multiple ChEMBL IDs**
Batch processing to retrieve information for several molecules simultaneously.

**What you'll enter:**
- **Multiple ChEMBL IDs**: Add several ChEMBL identifiers to a list
- **ID Management**: Add, remove, and validate multiple identifiers

**What you'll get:**
- Combined results for all requested molecules
- Side-by-side comparison capabilities
- Individual molecule selection and navigation
- Batch data export and analysis options

### **Preferred Name Search**
Find molecules using their common scientific or pharmaceutical names.

**What you'll enter:**
- **Preferred Name**: Scientific name of the compound (e.g., "Aspirin", "Sildenafil")

**What you'll get:**
- List of molecules matching the name
- Multiple candidates if the name is ambiguous
- Selection interface for choosing the correct molecule
- Complete molecular profiles for selected compounds

### **Synonym Search**
Discover molecules using brand names, trade names, or alternative identifiers.

**What you'll enter:**
- **Synonym**: Brand name, trade name, or alternative name (e.g., "Viagra", "Tylenol")

**What you'll get:**
- Molecules associated with the synonym
- Multiple matches if the synonym applies to several compounds
- Cross-reference information showing name relationships
- Access to the primary molecular data

## Usage
---
1. **Select Search Method**: Choose from the dropdown which type of search you want to perform
2. **Enter Search Information**:
   - For ChEMBL ID: Type the complete identifier (validated in real-time)
   - For Multiple IDs: Add each ChEMBL ID to your search list
   - For Names/Synonyms: Enter the common or brand name
3. **Execute Search**: The tool automatically searches as you type or when you submit
4. **Navigate Results**:
   - For single results: View the complete molecular profile
   - For multiple candidates: Select from the list of matches
   - For batch results: Browse through multiple molecules
5. **Explore Data**: Use the expandable data viewer to examine detailed information

## What You'll See
---
When you use the ChEMBL Getter, the interface displays several key sections:

### Search Interface
- **Search Type Dropdown**: Selection between different search methods
- **Input Fields**: Dynamic forms that change based on your search type
- **Validation Feedback**: Real-time checking of ChEMBL ID formats
- **Multiple ID Manager**: Add, remove, and organize ChEMBL IDs for batch searches

### Results Display
The results show comprehensive molecular information:

**Basic Information:**
- ChEMBL ID and preferred name
- Molecular type and classification
- Development phase and approval status

**Chemical Structure:**
- Canonical SMILES representation
- Standard InChI and InChI key
- Molecular formula and weight
- 2D structure visualization (molfile format)

**Molecular Properties:**
- Drug-likeness metrics (Lipinski's Rule of Five)
- Physical properties (LogP, polar surface area)
- Structural features (aromatic rings, hydrogen bond donors/acceptors)

**Cross-References:**
- Links to PubChem, PDB, UniProt
- Literature citations and patents
- Regulatory database references

**Synonyms and Names:**
- Trade names and brand names
- International nonproprietary names (INN)
- Research codes and alternative identifiers

## Example Searches
---

### Find Aspirin by ChEMBL ID
```
Search Type: ChEMBL ID
ChEMBL ID: CHEMBL25
```
*Returns complete information about aspirin including structure, properties, and cross-references*

### Search for Aspirin by Name
```
Search Type: Preferred Name
Preferred Name: Aspirin
```
*Finds aspirin and related compounds with similar names*

### Find Viagra by Brand Name
```
Search Type: Synonym
Synonym: Viagra
```
*Discovers sildenafil and related compounds using the brand name*

### Batch Analysis of Multiple Compounds
```
Search Type: Multiple ChEMBL IDs
ChEMBL IDs: CHEMBL25, CHEMBL192, CHEMBL1737
```
*Retrieves information for aspirin, sildenafil, and vardenafil simultaneously*

## Data Fields Description
---
The ChEMBL Getter displays rich molecular information:

### Structure Information
- **Canonical SMILES**: Simplified molecular representation
- **Standard InChI**: International chemical identifier
- **Molfile**: 2D structure coordinates and connectivity
- **Molecular Formula**: Chemical composition

### Drug Development
- **Max Phase**: Highest development phase reached
- **First Approval**: Year of first regulatory approval
- **Indication Class**: Therapeutic areas and uses
- **Withdrawn Information**: Withdrawal status and reasons

### Molecular Properties
- **Molecular Weight**: Exact and average molecular weights
- **LogP/LogD**: Lipophilicity measurements
- **Polar Surface Area**: Membrane permeability predictor
- **Rule of Five**: Drug-likeness assessment

### Cross-References
- **Database Links**: PubChem, ChEBI, DrugBank connections
- **Literature**: Associated publications and patents
- **Regulatory**: FDA, EMA, and other agency references

## Advanced Features
---

### Validation and Error Handling
- **Format Checking**: Real-time ChEMBL ID validation
- **Duplicate Detection**: Prevents adding the same molecule twice
- **Error Messages**: Clear feedback for invalid inputs or failed searches

### Multiple Result Handling
- **Candidate Selection**: Choose from multiple matches
- **Result Navigation**: Browse through batch search results
- **Comparison Mode**: Side-by-side molecular comparison

### Data Export and Integration
- **Structured Data**: Complete molecular profiles in JSON format
- **Cross-Database Links**: Direct access to external resources
- **Batch Processing**: Efficient handling of multiple molecules

## Use Cases
---
This tool is particularly valuable for:
- **Drug Discovery**: Finding existing drugs and analyzing their properties
- **Chemical Research**: Exploring molecular structures and relationships
- **Pharmaceutical Analysis**: Investigating drug development histories
- **Brand Name Resolution**: Connecting trade names to scientific identifiers
- **Batch Analysis**: Comparing multiple compounds simultaneously
- **Database Mining**: Exploring chemical space and molecular diversity
- **Regulatory Research**: Investigating approval status and indications
- **Academic Studies**: Educational exploration of medicinal chemistry
- **Patent Research**: Finding molecules by various naming conventions
- **Competitive Intelligence**: Analyzing pharmaceutical landscapes
