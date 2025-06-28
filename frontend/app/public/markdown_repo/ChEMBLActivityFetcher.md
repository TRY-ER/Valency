# ChEMBL Activity Fetcher Documentation
---

## Overview
---

The **ChEMBL Activity Fetcher** is a specialized bioactivity data retrieval tool that provides access to comprehensive biological activity information from the ChEMBL database. This tool enables researchers to explore bioactivity data from two perspectives: finding all compounds active against a specific biological target, or discovering all biological activities for a particular molecule. It's essential for drug discovery, target analysis, and understanding structure-activity relationships.

## Features
---
- **Dual Search Modes**: Target-based and molecule-based activity searches
- **Target Activity Search**: Find all compounds and their activities against a biological target
- **Molecule Activity Search**: Discover all recorded biological activities for a specific molecule
- **Activity Type Filtering**: Filter by specific bioactivity measurements (IC50, Ki, EC50, etc.)
- **pChEMBL Value Options**: Control inclusion of standardized potency measurements
- **Real-time Validation**: Automatic ChEMBL ID format checking for molecule searches
- **Interactive Results**: Browse and analyze multiple activity records
- **Comprehensive Data**: Detailed bioactivity profiles with target and compound information

### Tool Sections
---
1. **Search Configuration**
   - Activity search type selector (Target or Molecule)
   - Activity type filtering options
   - pChEMBL value inclusion controls
   - Input validation and guidance

2. **Input Interface**
   - Target name input for biological target searches
   - ChEMBL ID input with validation for molecule searches
   - Activity standard type selection dropdown
   - Clear error messaging and format guidance

3. **Results Display**
   - Target information and associated activity molecules
   - Molecule activity profiles with target details
   - Interactive activity record navigation
   - Comprehensive bioactivity data visualization

## Data Source
---
The tool accesses the ChEMBL bioactivity database, providing:

- **Target Information**: Biological target details, ChEMBL target IDs, and target classifications
- **Activity Records**: Comprehensive bioactivity measurements and experimental data
- **Standardized Values**: pChEMBL values for normalized potency comparisons
- **Assay Details**: Experimental conditions, protocols, and data quality indicators
- **Cross-References**: Links to targets, molecules, and related database entries
- **Activity Types**: IC50, Ki, EC50, Kd, and many other bioactivity measurements

## Search Types
---

### **Target Activities Search**
Find all compounds and their biological activities against a specific target.

**What you'll enter:**
- **Target Name**: Biological target name (e.g., "hERG", "Dopamine D2 receptor")
- **Standard Type**: Optional activity type filter (IC50, Ki, EC50, etc.)

**What you'll get:**
- Target information and details
- List of active compounds against the target
- Activity measurements and experimental data
- Structure visualization for active molecules
- Comparative activity analysis across compounds

**Use cases:**
- Drug discovery target analysis
- Finding known inhibitors or modulators
- Competitive intelligence research
- SAR analysis across compound series

### **Molecule Activities Search**
Discover all recorded biological activities for a specific ChEMBL molecule.

**What you'll enter:**
- **ChEMBL ID**: Molecule identifier (e.g., "CHEMBL25", "CHEMBL192")
- **pChEMBL Filter**: Include only activities with standardized potency values

**What you'll get:**
- Complete bioactivity profile for the molecule
- List of all tested biological targets
- Activity measurements across different assays
- Target selectivity and off-target effects
- Comprehensive pharmacological profile

**Use cases:**
- Compound profiling and characterization
- Safety assessment and off-target analysis
- Mechanism of action studies
- Lead optimization guidance

## Usage
---
1. **Select Search Mode**: Choose between Target Activities or Molecule Activities
2. **Configure Search Parameters**:
   - For Target Search: Enter target name and optionally select activity type
   - For Molecule Search: Enter ChEMBL ID and choose pChEMBL filtering
3. **Validate Input**: The tool automatically validates ChEMBL IDs and provides format guidance
4. **Execute Search**: Submit your query to retrieve bioactivity data
5. **Browse Results**:
   - For Target Search: Review active compounds and select molecules for detailed analysis
   - For Molecule Search: Explore activity records and navigate through different targets
   - Use arrow keys for quick navigation through results
6. **Analyze Data**: Examine detailed bioactivity information, experimental conditions, and relationships

## What You'll See
---
When you use the ChEMBL Activity Fetcher, the interface displays several key sections:

### Search Interface
- **Search Type Dropdown**: Selection between Target and Molecule activity searches
- **Input Fields**: Dynamic forms that change based on search type selection
- **Activity Type Selector**: Dropdown for filtering specific bioactivity measurements
- **Validation Indicators**: Real-time feedback for ChEMBL ID format checking
- **pChEMBL Options**: Checkbox for including standardized potency values

### Target Activities Results
The target search results show comprehensive target and activity information:

**Target Information:**
- Target ChEMBL ID and preferred name
- Target classification and organism
- Target type and functional description
- Cross-references to other databases

**Active Compounds:**
- List of molecules active against the target
- ChEMBL IDs and preferred names
- Activity measurements and standard types
- Structure visualization and molecular properties

**Activity Details:**
- Bioactivity values (IC50, Ki, EC50, etc.)
- Experimental conditions and assay information
- Data quality indicators and confidence scores
- Publication references and data sources

### Molecule Activities Results
The molecule search results show comprehensive bioactivity profiles:

**Activity Records:**
- List of all biological activities for the molecule
- Target information for each activity
- Activity measurements and standard types
- Experimental details and assay conditions

**Target Selectivity:**
- Comparison of activities across different targets
- Selectivity profiles and off-target effects
- Therapeutic area and target class analysis
- Safety and toxicity considerations

**Pharmacological Profile:**
- Comprehensive activity summary
- Mechanism of action insights
- Drug development implications
- Literature and patent references

## Example Searches
---

### Find hERG Channel Inhibitors
```
Search Type: Target Activities
Target Name: hERG
Standard Type: IC50
```
*Returns compounds that inhibit the hERG potassium channel with IC50 values*

### Discover Aspirin's Biological Activities
```
Search Type: Molecule Activities
ChEMBL ID: CHEMBL25
pChEMBL Value Exists: Yes
```
*Shows all recorded biological activities for aspirin across different targets*

### Explore Dopamine Receptor Modulators
```
Search Type: Target Activities
Target Name: Dopamine D2 receptor
Standard Type: All
```
*Finds all compounds tested against dopamine D2 receptor with any activity type*

### Profile Sildenafil's Selectivity
```
Search Type: Molecule Activities
ChEMBL ID: CHEMBL192
pChEMBL Value Exists: Yes
```
*Reveals sildenafil's activity across multiple targets for selectivity analysis*

## Data Fields Description
---
The ChEMBL Activity Fetcher displays rich bioactivity and target information:

### Activity Measurements
- **Standard Type**: Type of bioactivity measurement (IC50, Ki, EC50, Kd, etc.)
- **Standard Value**: Numerical activity value in standard units
- **Standard Units**: Units of measurement (nM, μM, mg/L, etc.)
- **pChEMBL Value**: Standardized potency measurement (-log10 of molar IC50/EC50/Ki)
- **Activity Comment**: Additional experimental notes and observations

### Target Information
- **Target ChEMBL ID**: Unique ChEMBL identifier for the biological target
- **Preferred Name**: Common or systematic name of the target
- **Target Type**: Classification (protein, enzyme, receptor, etc.)
- **Organism**: Species in which the target was tested
- **Target Classification**: Therapeutic area and protein family information

### Experimental Details
- **Assay ChEMBL ID**: Identifier for the specific bioassay
- **Assay Description**: Experimental protocol and conditions
- **Assay Type**: Type of assay (binding, functional, ADMET, etc.)
- **Confidence Score**: Data quality and reliability indicator
- **Document ChEMBL ID**: Reference to source publication or patent

### Molecular Information
- **Molecule ChEMBL ID**: Unique identifier for the tested compound
- **Preferred Name**: Common name of the molecule
- **Molecular Structure**: SMILES representation when available
- **Molecular Properties**: Basic physicochemical properties
- **Synonyms**: Alternative names and identifiers

## Advanced Features
---

### Activity Type Filtering
- **Comprehensive Options**: Filter by specific bioactivity types (IC50, Ki, EC50, etc.)
- **All Types Available**: Option to retrieve all activity measurements
- **Standard Descriptions**: Clear explanations of each activity type
- **Selective Analysis**: Focus on relevant bioactivity endpoints

### pChEMBL Value Control
- **Standardized Potency**: Include only activities with pChEMBL values
- **Quality Filtering**: Focus on high-quality, standardized measurements
- **Comparative Analysis**: Enable meaningful potency comparisons
- **Research Standards**: Align with pharmaceutical industry practices

### Results Navigation
- **Keyboard Controls**: Arrow key navigation through activity records
- **Interactive Selection**: Click-based browsing of molecules and activities
- **Data Organization**: Structured presentation of complex bioactivity data
- **Comprehensive Views**: Detailed information panels for selected records

### Validation and Error Handling
- **ChEMBL ID Validation**: Real-time format checking for molecule searches
- **Target Name Guidance**: Suggestions for proper target name formatting
- **Error Prevention**: Clear feedback prevents invalid searches
- **Data Availability**: Informative messages about search results

## Use Cases
---
This tool is particularly valuable for:
- **Drug Discovery**: Identifying lead compounds and analyzing target activity
- **Target Validation**: Understanding compound activity against specific targets
- **Compound Profiling**: Comprehensive bioactivity characterization
- **Safety Assessment**: Off-target activity and selectivity analysis
- **Competitive Intelligence**: Analyzing competitor compound activities
- **SAR Analysis**: Structure-activity relationship studies
- **Mechanism Studies**: Understanding compound mechanisms of action
- **Lead Optimization**: Guiding medicinal chemistry efforts
- **Academic Research**: Educational exploration of bioactivity data
- **Regulatory Submission**: Supporting drug development documentation
- **Biomarker Discovery**: Identifying activity patterns and target relationships
- **Therapeutic Area Research**: Exploring disease-relevant target activities
