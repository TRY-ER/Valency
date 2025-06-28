# PubChem Getter Documentation
---

## Overview
---

The **PubChem Getter** is a comprehensive chemical compound search and retrieval tool that provides access to the PubChem database, one of the world's largest collections of chemical information. This versatile tool allows users to search for and retrieve detailed compound information using multiple search methods including compound IDs, chemical names, SMILES strings, InChI keys, cross-references, and molecular mass queries.

## Features
---
- **Six Search Methods**: Multiple ways to find compounds in the PubChem database
- **PubChem CID Search**: Direct lookup using unique compound identifiers
- **Compound Name Search**: Find compounds by their chemical names
- **SMILES String Search**: Search using molecular structure notation
- **InChI Key Search**: Retrieve compounds using InChI identifiers
- **Cross-Reference Search**: Find compounds via external database references
- **Molecular Mass Search**: Discover compounds by mass ranges
- **Real-time Validation**: Automatic format checking for CIDs and mass values
- **Interactive Results**: Browse multiple candidates with structure visualization
- **Comprehensive Data**: Detailed compound information including properties and structures

### Tool Sections
---
1. **Search Configuration**
   - Search type selector dropdown
   - Additional parameter controls (cross-reference types, mass types)
   - Input validation and format guidance
   - Dynamic interface adaptation

2. **Input Interface**
   - Compound ID input with validation
   - Chemical name search fields
   - SMILES and InChI key input areas
   - Cross-reference type and value selection
   - Mass type and value specification

3. **Results Display**
   - Single compound detailed information
   - Multiple candidate selection interface
   - 2D molecular structure visualization
   - Comprehensive property data viewer

## Data Source
---
The tool accesses the PubChem database, providing comprehensive information about:

- **Chemical Structures**: 2D/3D coordinates, bond information, and atomic details
- **Molecular Properties**: Mass, formula, complexity, and physicochemical parameters
- **Identifiers**: SMILES, InChI, InChI keys, and various naming conventions
- **Cross-References**: Links to patents, registries, and other chemical databases
- **Calculated Properties**: LogP, polar surface area, hydrogen bond donors/acceptors
- **Structural Features**: Fingerprints, substructure keys, and molecular descriptors

## Search Types
---

### **PubChem CID Search**
Direct lookup of compounds using their unique PubChem Compound IDs.

**What you'll enter:**
- **CID**: Positive integer identifier (e.g., "2244", "5988")

**What you'll get:**
- Complete compound record with full structural data
- Detailed molecular properties and descriptors
- 2D structure visualization and analysis
- IUPAC names, SMILES, and InChI representations
- Comprehensive physicochemical properties

### **Compound Name Search**
Find compounds using their chemical or common names.

**What you'll enter:**
- **Chemical Name**: Common or systematic name (e.g., "Aspirin", "Acetylsalicylic acid")

**What you'll get:**
- List of matching compounds with multiple candidates
- Interactive selection interface for choosing the correct compound
- Structure visualization for selected candidates
- Complete compound profiles for chosen molecules

### **SMILES String Search**
Search for compounds using their SMILES (Simplified Molecular Input Line Entry System) notation.

**What you'll enter:**
- **SMILES**: Molecular structure string (e.g., "CC(=O)OC1=CC=CC=C1C(=O)O")

**What you'll get:**
- Compounds matching the exact SMILES structure
- Structural and property information
- Alternative representations and identifiers
- Related compound suggestions

### **InChI Key Search**
Retrieve compounds using their International Chemical Identifier keys.

**What you'll enter:**
- **InChI Key**: Standard chemical identifier (e.g., "BSYNRYMUTXBXSQ-UHFFFAOYSA-N")

**What you'll get:**
- Exact compound matches for the InChI key
- Complete molecular data and properties
- Structure visualization and analysis
- Cross-references to other databases

### **Cross-Reference Search**
Find compounds through external database references and identifiers.

**What you'll enter:**
- **Reference Type**: Database type (Patent ID, Registry ID, PubMed ID, etc.)
- **Reference Value**: Specific identifier in the chosen database

**What you'll get:**
- Compounds linked to the external reference
- Cross-database connectivity information
- Complete compound profiles
- Related patent and literature information

### **Molecular Mass Search**
Discover compounds by their molecular mass or weight.

**What you'll enter:**
- **Mass Type**: Exact mass, monoisotopic mass, or molecular weight
- **Mass Value**: Target mass in Daltons (e.g., "180.16")

**What you'll get:**
- Compounds matching the specified mass
- Mass-based similarity rankings
- Molecular formula and structural information
- Isomer and analog identification

## Usage
---
1. **Select Search Method**: Choose from the six available search types
2. **Configure Parameters**: Set additional options for cross-reference or mass searches
3. **Enter Search Data**: Input the compound identifier or search criteria
4. **Validate Input**: The tool automatically checks format and provides guidance
5. **Execute Search**: Submit your query to retrieve compound information
6. **Navigate Results**:
   - For single results: View complete compound data immediately
   - For multiple candidates: Select from the list to view detailed information
   - Use arrow keys for quick navigation through candidates
7. **Analyze Data**: Examine molecular structures, properties, and relationships

## What You'll See
---
When you use the PubChem Getter, the interface displays several key sections:

### Search Interface
- **Search Type Dropdown**: Selection between six different search methods
- **Parameter Controls**: Additional inputs for cross-reference types and mass types
- **Input Validation**: Real-time format checking with clear error messages
- **Dynamic Forms**: Interface adapts based on selected search type

### Results Display
The search results show comprehensive compound information:

**Single Compound View:**
- Complete PubChem compound record
- 2D molecular structure visualization
- Interactive property information panel
- Molecular descriptors and calculated values

**Multiple Candidates View:**
- List of matching compounds with navigation
- Candidate selection interface
- Structure preview for selected compounds
- Detailed information for chosen molecules

**Compound Information:**
- IUPAC names in multiple formats
- Molecular formula and exact mass
- SMILES and InChI representations
- Physicochemical properties (LogP, TPSA, etc.)
- Hydrogen bond donors and acceptors
- Molecular complexity and fingerprints
- Cross-references to external databases

## Example Searches
---

### Find Aspirin by CID
```
Search Type: PubChem CID
CID: 2244
```
*Returns complete information about aspirin including structure and properties*

### Search for Caffeine by Name
```
Search Type: Compound Name
Name: Caffeine
```
*Finds caffeine and related compounds with similar names*

### Find Ethanol by SMILES
```
Search Type: SMILES String
SMILES: CCO
```
*Discovers ethanol using its molecular structure notation*

### Search by InChI Key
```
Search Type: InChI Key
InChI Key: BSYNRYMUTXBXSQ-UHFFFAOYSA-N
```
*Retrieves aspirin using its international chemical identifier*

### Find Compounds by Patent
```
Search Type: Cross-Reference
Reference Type: Patent ID
Reference Value: US1234567
```
*Discovers compounds referenced in specific patents*

### Search by Molecular Weight
```
Search Type: Molecular Mass
Mass Type: Molecular Weight
Mass Value: 180.16
```
*Finds compounds with approximately 180.16 Da molecular weight*

## Data Fields Description
---
The PubChem Getter displays rich compound information:

### Structural Information
- **Atoms**: Atomic coordinates, elements, and connectivity
- **Bonds**: Bond orders, types, and molecular connectivity
- **Coordinates**: 2D and 3D structural coordinates
- **Conformers**: Alternative molecular conformations
- **Charge**: Overall molecular charge state

### Molecular Identifiers
- **CID**: Unique PubChem Compound identifier
- **SMILES**: Canonical, absolute, and isomeric representations
- **InChI**: Standard international chemical identifier
- **InChI Key**: Condensed chemical identifier hash
- **IUPAC Names**: Multiple systematic naming conventions

### Calculated Properties
- **Molecular Weight**: Average and monoisotopic masses
- **Molecular Formula**: Chemical composition
- **LogP**: Lipophilicity measurement (XLogP3)
- **TPSA**: Topological polar surface area
- **Complexity**: Molecular structural complexity score
- **Heavy Atoms**: Count of non-hydrogen atoms

### Chemical Features
- **Hydrogen Bond Donors**: Count of HB donor sites
- **Hydrogen Bond Acceptors**: Count of HB acceptor sites
- **Rotatable Bonds**: Molecular flexibility indicator
- **Chiral Centers**: Stereochemical information
- **Aromatic Systems**: Aromatic ring analysis

### Database Cross-References
- **Fingerprints**: Substructure keys and molecular fingerprints
- **External IDs**: Links to other chemical databases
- **Patent References**: Associated patent information
- **Literature**: Connected publications and research

## Advanced Features
---

### Multiple Search Types
- **Flexible Input**: Six different ways to search for compounds
- **Format Validation**: Real-time checking of CIDs, SMILES, and InChI keys
- **Parameter Adaptation**: Dynamic interface based on search type
- **Cross-Database Integration**: Links to external chemical resources

### Results Navigation
- **Candidate Selection**: Interactive browsing of multiple search results
- **Keyboard Navigation**: Arrow key support for quick candidate switching
- **Structure Visualization**: 2D molecular structure display
- **Property Analysis**: Comprehensive molecular descriptor examination

### Data Validation
- **Format Checking**: Automatic validation of compound identifiers
- **Error Prevention**: Clear guidance for proper input formatting
- **Search Optimization**: Intelligent query processing and result ranking
- **Quality Assurance**: Verified compound data from PubChem

### Visualization Tools
- **2D Structure Display**: Interactive molecular structure viewer
- **Property Panels**: Organized display of molecular characteristics
- **Information Hierarchy**: Structured presentation of complex data
- **Export Capabilities**: Access to raw compound data

## Use Cases
---
This tool is particularly valuable for:
- **Chemical Research**: Finding detailed compound information and properties
- **Drug Discovery**: Identifying potential drug candidates and analogs
- **Academic Studies**: Educational exploration of chemical structures
- **Patent Research**: Finding compounds referenced in intellectual property
- **Regulatory Submission**: Accessing official compound data for documentation
- **Chemical Safety**: Obtaining molecular property data for risk assessment
- **Synthesis Planning**: Finding structural information for chemical synthesis
- **Database Mining**: Systematic exploration of chemical space
- **Quality Control**: Verifying compound identity and purity
- **Computational Chemistry**: Obtaining structural data for modeling
- **Medicinal Chemistry**: Analyzing molecular properties for drug design
- **Chemical Biology**: Understanding molecular interactions and properties
