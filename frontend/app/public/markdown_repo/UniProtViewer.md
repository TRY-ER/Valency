# UniProt Viewer Documentation
---

## Overview
---

The **UniProt Viewer** tool allows users to fetch and visualize protein structure data from the AlphaFold Protein Structure Database. By providing a UniProt Accession Key, the tool retrieves detailed information, 3D models, and sequence data for the specified protein. It presents this information through an interactive interface comprising several specialized sections.

## Features
---
- **UniProt Accession Input**: A dedicated input section where users can enter a UniProt Accession Key. The input includes validation to ensure the key is in the correct format.
- **Information Panel**: Displays comprehensive details about the protein, including:
    - Gene Name
    - Protein Description
    - Organism (Scientific Name)
    - UniProt ID
    - Entry ID
    - Sequence Version Date
    - Model Created Date
    - Latest Version
    - Tax ID
- **3D Viewer**: An interactive 3D visualization of the protein structure, loaded from the PDB URL obtained via the AlphaFold API.
- **Sequence Viewer**: Displays the protein's amino acid sequence using the RCSB Saguaro feature viewer, allowing for an interactive exploration of the sequence.
- **File Downloader**: Provides direct download links for various files associated with the protein model, such as:
    - PDB File (.pdb)
    - CIF File (.cif)
    - BCIF File (.bcif)
    - PAE Image (.png) - Predicted Aligned Error Image
    - PAE Data (.json) - Predicted Aligned Error Data

## API Endpoint
---
The UniProt Viewer utilizes the AlphaFold Protein Structure Database API to fetch protein data. More information about the AlphaFold service can be found at [https://alphafold.ebi.ac.uk/](https://alphafold.ebi.ac.uk/).

**Endpoint**: `GET /prediction/{qualifier}`

**Description**: Retrieves all AlphaFold models and associated metadata for a given UniProt accession.

**Parameters**:
- `qualifier` (path, string, required): The UniProt accession key for the desired protein. Example: `Q5VSL9`.
- `sequence_checksum` (query, string, optional): CRC64 checksum of the UniProt sequence.

The API returns a JSON array containing objects with detailed information for each model.

## Usage
---
To use the `UniProt Viewer`:
1. Navigate to the UniProt Viewer section within the application.
2. Enter a valid UniProt Accession Key into the input field (e.g., `Q5VSL9`).
3. The tool will automatically fetch the data from the AlphaFold API.
4. Once the data is retrieved:
    - The **Information Panel** will populate with protein details.
    - The **3D Viewer** will render the protein structure.
    - The **Sequence Viewer** will display the protein's amino acid sequence.
    - The **File Downloader** section will show links to download associated files.

## Example
---
You can try the UniProt Viewer with the following UniProt Accession Key:

```
Q5VSL9
```

Upon entering this key, the tool will fetch and display the relevant information and visualizations for the corresponding protein.
