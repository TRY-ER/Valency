# UniProt Summary Viewer Documentation
---

## Overview
---

The **UniProt Summary Viewer** is a specialized protein analysis tool that provides comprehensive summary information about proteins and their structural data. This tool retrieves detailed metadata about protein structures, experimental methods, and associated structural files from the UniProt and AlphaFold databases, presenting the information in an easy-to-navigate format.

## Features
---
- **Comprehensive Summary Data**: Displays detailed protein entry information including sequence details and structural coverage
- **Structure Information**: Shows available protein structures with metadata about experimental methods and quality scores
- **Interactive Data Viewer**: Expandable JSON data viewer for exploring complete API responses
- **File Downloads**: Direct access to structural files in various formats
- **Real-time Search**: Instant retrieval of summary data upon entering UniProt accessions

### Tool Sections
---
1. **Search Input**
   - UniProt accession key input field
   - "Fetch Summary" button for manual data retrieval
   - Loading indicator during data fetching
   - Error messages for invalid inputs

2. **Complete Summary Data Viewer**
   - Expandable/collapsible JSON data display
   - Full API response visualization
   - Structured data presentation
   - Search and navigation within large datasets

3. **Associated Files Section**
   - Download links for structural files
   - Direct access to CIF format files
   - File format indicators and descriptions

## Data Source
---
The tool retrieves protein summary information from the UniProt database through the AlphaFold API, providing:

- **Entry Information**: UniProt accession details, sequence length, and segment information
- **Structural Data**: Available protein structures with quality metrics
- **Experimental Details**: Methods used for structure determination
- **Coverage Information**: How much of the protein sequence is covered by available structures
- **Confidence Scores**: Quality assessments for predicted structures
- **File Access**: Direct links to downloadable structural data files

## Usage
---
1. **Enter a UniProt Accession**: Type or paste a UniProt accession key into the search field (e.g., "Q5VSL9")
2. **Click "Fetch Summary"**: Press the button to retrieve comprehensive summary data
3. **Explore Summary Data**: 
   - Click to expand the "Complete UniProt Summary Data" section
   - Browse through the structured JSON response
   - Navigate different sections of the summary information
4. **Download Files**: Click on available file download icons to access structural data
5. **Try Different Proteins**: Use the search field to explore different protein summaries

## What You'll See
---
When you use the UniProt Summary Viewer, the interface displays several key sections:

### Search Interface
- **Input Field**: Text box for entering UniProt accession numbers
- **Fetch Button**: Button to initiate the summary data retrieval
- **Loading Indicator**: Spinning animation with status text during data fetching
- **Error Messages**: Clear feedback if the accession is invalid or if there are connection issues

### Summary Data Display
The main data viewer shows comprehensive information including:

**UniProt Entry Details:**
- **Accession (AC)**: The unique UniProt identifier
- **Entry ID**: Alternative identifier for the protein
- **Checksum**: Sequence verification code
- **Sequence Length**: Total number of amino acids
- **Segment Information**: Start and end positions of analyzed segments

**Structure Information:**
For each available structure, you'll see:
- **Model Identifier**: Unique ID for the structural model
- **Model Category**: Whether experimentally determined or computationally predicted
- **Provider**: Source database or organization
- **Experimental Method**: Technique used (X-ray crystallography, NMR, Cryo-EM, etc.)
- **Resolution**: Quality measure for experimental structures
- **Confidence Scores**: Quality metrics for predicted structures
- **Coverage**: Percentage of protein sequence covered by the structure
- **Creation Date**: When the structure was deposited
- **Oligomeric State**: Whether the protein functions as a monomer, dimer, etc.

### File Downloads Section
Available file formats for download:
- **CIF Files**: Crystallographic Information Files containing detailed structural data
- **Model URLs**: Direct links to 3D structure files
- **Ensemble Data**: Multiple conformations for flexible structures

## Data Fields Description
---
The summary viewer displays rich metadata organized into categories:

### Entry Information
- **AC (Accession)**: Primary UniProt identifier
- **ID**: Entry name in UniProt
- **Uniprot Checksum**: Sequence verification code
- **Sequence Length**: Total amino acids in the protein
- **Segment Start/End**: Specific region analyzed

### Structure Details
- **Model Category**: EXPERIMENTALLY DETERMINED or PREDICTED
- **Model Format**: File format (typically PDB)
- **Model Type**: ATOMIC level detail
- **Sequence Identity**: Similarity between structure and UniProt sequence
- **Coverage**: Fraction of protein covered by structure
- **Confidence Type**: Scoring method (e.g., pLDDT for AlphaFold)
- **Resolution**: Experimental resolution in Angstroms
- **Experimental Method**: X-ray, NMR, Electron Crystallography, etc.

## Example
---
### Sample UniProt Accessions to Try
Explore these example proteins to see comprehensive summary data:

```
Q5VSL9
```
*Human protein with multiple structural entries*

```
P04637
```
*Tumor suppressor p53 - extensively studied protein*

```
P0DTC2
```
*SARS-CoV-2 spike protein - recent structural studies*

### What Happens When You Search
1. Type "Q5VSL9" into the search field
2. Click "Fetch Summary" button
3. Watch the loading indicator while data is retrieved
4. The tool displays:
   - Complete summary data in an expandable viewer
   - Detailed structure information for all available models
   - Download options for structural files
   - Quality metrics and experimental details

## Use Cases
---
This tool is particularly valuable for:
- **Structure Analysis**: Comparing different structural models of the same protein
- **Quality Assessment**: Evaluating confidence scores and experimental methods
- **Research Planning**: Identifying the best structural data for specific research needs
- **Database Mining**: Exploring comprehensive metadata for protein families
- **Method Comparison**: Understanding different experimental approaches used
- **Coverage Analysis**: Determining how well a protein sequence is structurally characterized
- **File Management**: Accessing and organizing structural data files
- **Educational Research**: Learning about protein structure determination methods
