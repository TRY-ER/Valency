# BRICS PSMILES Generator Documentation
---

## Overview
---

The **BRICS Polymer Generator** tool enables users to generate novel polymer structures using the Breaking Retrosynthetically Interesting Chemical Substructures (BRICS) algorithm. The tool takes PSMILES (Polymer SMILES) strings as input and generates hypothetical polymers by recombining monomer fragments.

## Features
---
- **Batch Input**: Accepts multiple PSMILES strings via text file upload
- **Progress Tracking**: Visual progress bar showing completion status 
- **File Preview**: Preview panel to verify uploaded PSMILES strings
- **Results Panel**: Interactive output display showing generated polymers
- **Export Options**: Download generated structures or reset the process

### Tool Sections
---
1. **Progress Section**
   - Progress bar showing completed steps
   - Current step indication
   - Overall process status

2. **Input Section**
   - File upload interface
   - PSMILES preview panel
   - File validation feedback

3. **Output Section**
   - Generated polymers display
   - Download results button
   - Reset operation button
   - Step-by-step process details

## Usage
---
1. Prepare a text file containing PSMILES strings (one per line)
2. Upload the file through the interface
3. Verify the PSMILES strings in the preview panel
4. Monitor generation progress
5. Download or copy the generated polymers

## Example
---
Input file format example:

```
[*]CC[*]
[*]CC(C)[*]
[*]CCC(=O)[*]
```

The tool will process these PSMILES strings and generate new polymers by:
1. Breaking them into fragments
2. Recombining fragments in novel ways
3. Generating valid polymer structures

### Use Cases
---
This tool is particularly useful when:
- Exploring polymer space around known monomers
- Generating libraries of potential new materials
- Discovering new polymer scaffolds
- Expanding existing polymer databases