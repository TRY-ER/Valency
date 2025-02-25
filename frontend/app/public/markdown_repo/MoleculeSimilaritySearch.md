# Molecule Similarity Search Documentation
---

## Overview
---

The **Similarity Search** tool enables users to find similar molecules by providing a SMILES string query. The tool returns a collection of structurally similar molecules ranked by their similarity scores.

## Features
---
- **SMILES Input**: Takes a SMILES string as the query molecule.
- **Candidate Count**: Allows users to specify the number of similar molecules to retrieve.
- **Results Visualization**: Displays 2D structures of similar molecules with their SMILES strings.
- **Similarity Metrics**: Shows similarity distance scores for each retrieved molecule.
- **Export Options**: Results can be downloaded as a text file or copied to clipboard.

### Results Display
---
For each similar molecule found, the following information is shown:
1. **2D Structure**: Visual representation of the molecule
2. **SMILES String**: The molecular structure in SMILES format
3. **Similarity Score**: Numerical value indicating structural similarity to the query
4. **Rank**: Position in the similarity-ordered results

## Usage
---
The tool can be used by:
1. Entering a valid SMILES string
2. Specifying the desired number of similar molecules to retrieve
3. Initiating the search
4. Viewing and optionally exporting the results

## Example
---
You can search for similar molecules like so:

```
Query SMILES: CCO
Number of candidates: 10
```

This will return molecules similar to ethanol, ranked by their structural similarity.

### Export Options
---
- **Download**: Save results as a text file
- **Copy**: Copy all results to clipboard for further use