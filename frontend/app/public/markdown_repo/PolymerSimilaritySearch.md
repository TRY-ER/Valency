# Polymer Similarity Search Documentation
---

## Overview
---

The **Polymer Similarity Search** tool enables users to find similar polymers by providing a PSMILES (Polymer SMILES) string query. The tool returns a collection of structurally similar polymers ranked by their similarity scores.

## Features
---
- **PSMILES Input**: Takes a PSMILES string as the query polymer monomer
- **Candidate Count**: Allows users to specify the number of similar polymers to retrieve
- **Results Visualization**: Displays 2D structures of similar polymer monomers with their PSMILES strings
- **Similarity Metrics**: Shows similarity distance scores for each retrieved polymer
- **Export Options**: Results can be downloaded as a text file or copied to clipboard

### Results Display
---
For each similar polymer found, the following information is shown:
1. **2D Structure**: Visual representation of the polymer monomer
2. **PSMILES String**: The monomer structure in PSMILES format
3. **Similarity Score**: Numerical value indicating structural similarity to the query
4. **Rank**: Position in the similarity-ordered results

## Usage
---
The tool can be used by:
1. Entering a valid PSMILES string
2. Specifying the desired number of similar polymers to retrieve
3. Initiating the search
4. Viewing and optionally exporting the results

## Example
---
You can search for similar polymers like so:

```
Query PSMILES: [*]CC[*]
Number of candidates: 10
```

This will return polymers similar to polyethylene, ranked by their structural similarity.

### Export Options
---
- **Download**: Save results as a text file
- **Copy**: Copy all results to clipboard for further use