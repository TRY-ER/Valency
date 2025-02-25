# Protein Similarity Search Documentation
---

## Overview
---

The **Protein Similarity Search** tool enables users to find similar proteins by providing a PDB ID (Protein Data Bank ID). The tool returns a collection of structurally similar proteins ranked by their similarity scores.

## Features
---
- **PDB ID Input**: Takes a PDB ID as the query protein
- **Candidate Count**: Allows users to specify the number of similar proteins to retrieve
- **Results Visualization**: Displays 3D structures of similar proteins
- **Similarity Metrics**: Shows similarity scores for each retrieved protein
- **Export Options**: Results can be downloaded as a text file or copied to clipboard

### Results Display
---
For each similar protein found, the following information is shown:
1. **3D Structure**: Interactive visualization of the protein structure
2. **PDB ID**: The unique identifier for the protein
3. **Similarity Score**: Numerical value indicating structural similarity to the query
4. **Basic Properties**: Key information such as:
   - Title
   - Authors
   - Release Date
   - Resolution
   - Experimental Method

## Usage
---
The tool can be used by:
1. Entering a valid PDB ID
2. Specifying the desired number of similar proteins to retrieve
3. Initiating the search
4. Viewing and optionally exporting the results

## Example
---
You can search for similar proteins like so:

```
Query PDB ID: 1YCR
Number of candidates: 10
```

This will return proteins with similar structural features to the query protein, ranked by their similarity scores.

### Export Options
---
- **Download**: Save results as a text file
- **Copy**: Copy all results to clipboard for further use