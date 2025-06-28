# Protein Explorer Documentation
---

## Overview
---

The **Protein Explorer** tool allows users to input a Protein Data Bank ID (**PDB ID**) for a given protein and provides an interactive interface for visualizing its 3D structure. In addition, the component displays key protein details in an information panel.

## Features
---
- **3D Protein Visualization**: A rotatable, scalable 3D view of the protein structure.
- **Information Panel**: Displays additional protein details, such as authors, journal references, molecular weight, organism data, and more.

### Protein Properties Displayed
---
1. **PDB ID**: The unique identifier for the protein structure.  
2. **Title, Authors, Journal, Year**: Citation details for the protein.  
3. **Experiment Method**: Describes how the protein structure was determined (e.g., X-ray crystallography, NMR).  
4. **Molecular Weight (kDa)**: The total molecular weight of the protein.  
5. **Deposited Model Count**: How many structural models are stored under this PDB ID.  
6. **Polymer Entity Count**: Number of distinct polymer chains in the protein.  
7. **Polymer Monomer Count**: Total monomers that make up each polymer chain.  
8. **Resolution**: Indicates the level of detail in the determined structure.  
9. **Release Date**: When the structure was deposited or released.  

## Usage
---
The `Protein Explorer` component can be used by providing a PDB ID. Once the PDB ID is supplied, the tool will display a 3D model of the protein and populate the information panel with relevant properties.