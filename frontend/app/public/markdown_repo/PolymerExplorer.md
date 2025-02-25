# Polymer Explorer Documentation
---

## Overview
---

The **Polymer Explorer** tool accepts a PSMILES (Polymer SMILES) string and provides an interactive interface for visualizing the polymer. It also displays key details about the polymer’s monomer in an information panel.

## Features
---
- **PSMILES Input**: The tool accepts a PSMILES string as input to represent the monomer structure.  
- **Interactive Visualization**: A side panel displays the polymer monomer, optionally highlighting potential bonding sites for polymer chain formation.  
- **Information Panel**: Shows relevant monomer properties such as molecular formula, molecular weight, and ring count.

### Polymer Monomer Properties Displayed
---
1. **Molecular Formula**: The chemical formula of the monomer.  
2. **Monomer Molecular Weight**: The total molecular weight of the monomer.  
3. **Number of Rings**: The count of ring structures present in the monomer.  
4. **Open Bond Indexes**: Positions where the next monomer can attach to form the polymer chain.  

## Usage
---
The `Polymer Explorer` component can be used by providing a PSMILES string. Once the PSMILES is submitted, the tool renders the monomer structure in a side panel and populates the information panel with relevant properties.

## Example
---
You can pass a PSMILES string like so:

```
[*]CC[*]
```

Here, `"[*]CC[*]"` is a sample PSMILES string for poly-ethylene.