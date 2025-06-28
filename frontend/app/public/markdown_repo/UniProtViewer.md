# UniProt Viewer Documentation
---

## Overview
---

The **UniProt Viewer** is a comprehensive protein visualization tool that displays detailed protein information from UniProt databases integrated with AlphaFold structure predictions. This interactive component allows users to explore protein data, view 3D molecular structures, and access associated files for further analysis.

## Features
---
- **Protein Information Display**: Shows comprehensive protein metadata including gene information, organism details, and UniProt identifiers
- **3D Structure Visualization**: Interactive 3D protein structure viewer powered by AlphaFold predictions
- **File Access**: Direct download links for various protein data formats (PDB, CIF, BCIF, PAE data)

### Tool Sections
---
1. **Input Section**
   - UniProt accession key input field
   - Real-time validation and feedback
   - Optional hide/show functionality

2. **Protein Information Panel**
   - Gene name and description
   - Organism scientific name
   - UniProt ID and entry details
   - Sequence version and model creation dates
   - Taxonomic information

3. **3D Structure Viewer**
   - Interactive molecular structure visualization
   - Powered by AlphaFold predictions
   - Supports PDB format structures
   - Zoom, rotate, and pan functionality

4. **Associated Files Section**
   - Download links for multiple file formats
   - PDB files for structure data
   - CIF/BCIF files for crystallographic data
   - PAE (Predicted Aligned Error) data and images
   - Direct file access with download icons

## Data Source
---
The tool retrieves protein information from the AlphaFold Protein Structure Database, which provides AI-predicted protein structures. The system automatically fetches:

- **Basic Protein Information**: Gene name, description, and organism details
- **Structural Data**: 3D protein models and coordinates
- **Sequence Information**: Complete amino acid sequences
- **Quality Metrics**: Prediction confidence scores and error estimates
- **File Downloads**: Multiple file formats for further analysis

## Usage
---
1. **Enter a UniProt Accession**: Type or paste a UniProt accession key into the search field (e.g., "Q5VSL9")
2. **View Protein Information**: The tool automatically displays detailed protein metadata including:
   - Gene name and protein description
   - Source organism information
   - UniProt identifiers and version details
   - Model creation dates
3. **Explore 3D Structure**: Interact with the 3D protein model using mouse controls:
   - **Left click + drag**: Rotate the structure
   - **Scroll wheel**: Zoom in and out
   - **Right click + drag**: Pan the view
4. **Download Files**: Click on any file type icon to download associated data files
5. **Try Different Proteins**: Use the input field to explore different protein structures

## What You'll See
---
When you use the UniProt Viewer, the interface displays several information panels:

### Protein Information Panel
Displays key details about the selected protein:
- **Gene**: The gene name that codes for this protein
- **Description**: What the protein does and its biological function
- **Organism**: The scientific name of the species this protein comes from
- **UniProt ID**: The unique identifier in the UniProt database
- **Entry ID**: Additional database reference
- **Sequence Version Date**: When the protein sequence was last updated
- **Model Created Date**: When the AlphaFold prediction was generated
- **Latest Version**: Version number of the current model
- **Tax ID**: Taxonomic identifier for the organism

### 3D Structure Display
An interactive molecular viewer showing:
- The predicted 3D structure of the protein
- Color-coded regions indicating prediction confidence
- Ability to rotate, zoom, and examine the structure from all angles
- Real-time rendering of the molecular model

### File Downloads Section
Icons and links for downloading various file formats:
- **PDB File**: Standard 3D structure format for molecular visualization software
- **CIF File**: Crystallographic format with detailed structural information
- **BCIF File**: Binary compressed version of CIF files
- **PAE Image**: Visual representation of prediction confidence
- **PAE Data**: Raw confidence scores in JSON format

## File Formats Supported
---
1. **PDB Files (.pdb)**
   - 3D structure coordinates
   - Standard protein structure format
   - Compatible with most molecular viewers

2. **CIF Files (.cif)**
   - Crystallographic Information File
   - More detailed structural information
   - Includes experimental data

3. **BCIF Files (.bcif)**
   - Binary CIF format
   - Compressed version of CIF
   - Faster loading for large structures

4. **PAE Data (.json, .png)**
   - Predicted Aligned Error information
   - Confidence scores for structure prediction
   - Visual representation of prediction reliability

## Example
---
### Sample UniProt Accessions to Try
Explore these example proteins to see how the tool works:

```
Q5VSL9
```
*Human protein - good starting example*

```
P04637
```
*Tumor suppressor p53 - important cancer research protein*

```
P0DTC2
```
*SARS-CoV-2 spike protein - virus research example*

### What Happens When You Enter an Accession
1. Type "Q5VSL9" into the search field
2. The tool will automatically fetch the data and display:
   - Protein information in the details panel
   - Interactive 3D structure on the right
   - Download options at the bottom
3. You can then rotate the 3D model, read about the protein, and download files for further analysis

## Use Cases
---
This tool is particularly useful for:
- **Protein Research**: Detailed analysis of protein structures and properties
- **Drug Discovery**: Examining target proteins and binding sites
- **Structural Biology**: Comparing predicted vs. experimental structures
- **Educational Purposes**: Teaching protein structure and function
- **Bioinformatics Workflows**: Integration with protein analysis pipelines
- **Quality Assessment**: Evaluating AlphaFold prediction confidence through PAE data
