# AlphaFold Annotations Viewer Documentation
---

## Overview
---

The **AlphaFold Annotations Viewer** tool allows users to fetch and visualize annotations for proteins from the AlphaFold database. By providing a UniProt Accession Key and optionally specifying residue ranges, the tool retrieves detailed annotation information and presents it through an interactive interface with visualization tracks.

## Features
---
- **UniProt Accession Input**: A dedicated input section where users can enter a UniProt Accession Key along with optional start and end positions for residue ranges.
- **Annotations Summary**: Displays a summary count of different annotation types found for the protein.
- **Interactive Annotations Viewer**: Uses RCSB Saguaro to display annotations in colored tracks, allowing users to:
  - Scroll and zoom to explore different regions
  - View different annotation types in separate color-coded tracks
  - See detailed annotation information
- **Complete Data Viewer**: Provides access to the full API response data in a collapsible, structured format.

### Annotation Types Displayed
---
The viewer can display various types of annotations depending on what's available for the protein, such as:
- Functional domains
- Structural features
- Binding sites
- Secondary structure elements
- Post-translational modifications
- And other annotation types provided by AlphaFold

## API Endpoint
---
The AlphaFold Annotations Viewer utilizes the AlphaFold Protein Structure Database API to fetch annotation data. More information about the AlphaFold service can be found at [https://alphafold.ebi.ac.uk/](https://alphafold.ebi.ac.uk/).

**Endpoint**: `POST /mcp/alphafold/get_alphafold_annotations`

**Description**: Retrieves AlphaFold annotations for a UniProt residue range.

**Parameters**:
- `uniprot_accession` (string, required): The UniProt accession ID for the protein. Example: `Q5VSL9`.
- `start` (number, optional): The start position of the residue range.
- `end` (number, optional): The end position of the residue range.

The API returns a JSON object containing detailed annotation information for the specified protein and range.

## Usage
---
To use the `AlphaFold Annotations Viewer`:
1. Navigate to the AlphaFold Annotations Viewer section within the application.
2. Enter a valid UniProt Accession Key into the input field (e.g., `Q5VSL9`).
3. Optionally, specify start and end positions to focus on a specific residue range.
4. Click "Fetch Annotations" to retrieve the data from the AlphaFold API.
5. Once the data is retrieved:
    - The **Annotations Summary** will show counts of different annotation types.
    - The **Annotations Viewer** will display interactive tracks with color-coded annotations.
    - The **Complete Data Viewer** will provide access to the full API response.

## Example
---
You can try the AlphaFold Annotations Viewer with the following UniProt Accession Key:

```
Q5VSL9
```

You can also specify a residue range:
- Start Position: `1`
- End Position: `100`

Upon entering this information and fetching annotations, the tool will display the relevant annotation information and visualizations for the corresponding protein region.

## Interactive Features
---
- **Zoom and Pan**: Use the mouse to scroll and zoom within the annotations viewer to explore different regions of the protein.
- **Color-Coded Tracks**: Different annotation types are displayed in separate tracks with distinct colors for easy identification.
- **Detailed Information**: Hover over annotations to see detailed descriptions and positions.
- **Data Export**: Access the complete annotation data through the expandable data viewer for further analysis.
