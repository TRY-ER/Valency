# AlphaFold Annotations Viewer Documentation
---

## Overview
---

The **AlphaFold Annotations Viewer** is a specialized protein analysis tool that visualizes functional annotations mapped onto protein sequences. This tool allows users to explore specific types of annotations such as mutagenesis data, binding sites, and other functional features, providing both interactive visual representations and detailed annotation data for comprehensive protein analysis.

## Features
---
- **Dual Input System**: Accepts both UniProt accession keys and annotation type specifications
- **Interactive Annotations Viewer**: Visual sequence track display with color-coded annotation regions
- **Annotation Type Filtering**: Search for specific annotation types (e.g., MUTAGEN, BINDING, DOMAIN)
- **Summary Statistics**: Quick overview of annotation counts by type
- **Visual Track Display**: Interactive sequence viewer with zoom and navigation capabilities
- **Complete Data Access**: Expandable JSON viewer for detailed annotation information
- **Real-time Search**: Instant annotation retrieval and visualization

### Tool Sections
---
1. **Input Interface**
   - UniProt accession key input field
   - Annotation type specification field
   - "Fetch Annotations" button with loading feedback
   - Error messages and validation

2. **Annotations Summary Panel**
   - Count display for different annotation types
   - Grid layout showing annotation distribution
   - Quick statistics overview

3. **Interactive Annotations Viewer**
   - Sequence-based track visualization
   - Color-coded annotation regions
   - Zoom and pan functionality
   - Detailed annotation labels

4. **Complete Data Viewer**
   - Expandable JSON data display
   - Full API response exploration
   - Structured annotation details

## Data Source
---
The tool retrieves annotation information from the AlphaFold database, providing:

- **Sequence Information**: Complete protein sequence with positional mapping
- **Annotation Details**: Functional annotations with precise start/end positions
- **Evidence Types**: Source and evidence classification for each annotation
- **Annotation Values**: Specific values and measurements for quantitative annotations
- **Source Attribution**: Original databases and publications for each annotation
- **Regional Data**: Continuous regions with associated functional information

## Usage
---
1. **Enter Protein Information**: 
   - Type a UniProt accession key (e.g., "Q5VSL9")
   - Specify annotation type (e.g., "MUTAGEN")
2. **Click "Fetch Annotations"**: Press the button to retrieve annotation data
3. **View Summary Statistics**: Check the overview panel for annotation counts
4. **Explore Interactive Viewer**:
   - **Mouse scroll**: Zoom in and out on the sequence track
   - **Click and drag**: Pan along the sequence
   - **Hover**: View detailed annotation information
5. **Examine Detailed Data**: Expand the complete data viewer for comprehensive information
6. **Try Different Types**: Change annotation type to explore different functional features

## What You'll See
---
When you use the AlphaFold Annotations Viewer, the interface displays several specialized sections:

### Input Controls
- **UniProt Accession Field**: Text input for protein identifier
- **Annotation Type Field**: Specification of annotation category to retrieve
- **Fetch Button**: Initiates the annotation data retrieval
- **Loading Indicator**: Visual feedback during data fetching
- **Error Display**: Clear messages for invalid inputs or connection issues

### Annotations Summary
A grid display showing:
- **Annotation Counts**: Number of each annotation type found
- **Type Categories**: Different functional annotation classes
- **Quick Statistics**: Overview of annotation distribution
- **Visual Organization**: Grid layout for easy comparison

### Detailed Information Display
For each annotation, you can view:
- **Position Information**: Exact start and end residue numbers
- **Annotation Description**: Detailed functional information
- **Evidence Type**: Source and confidence classification
- **Source Attribution**: Original database or publication
- **Annotation Values**: Specific measurements or classifications
- **Regional Details**: Continuous functional regions

## Data Fields Description
---
The annotation viewer displays comprehensive information:

### Protein Information
- **Accession**: UniProt identifier for the protein
- **ID**: UniProt entry name
- **Sequence**: Complete amino acid sequence

### Annotation Details
- **Type**: Category of annotation (MUTAGEN, BINDING, etc.)
- **Description**: Detailed functional description
- **Source Name**: Originating database or study
- **Source URL**: Link to original data source
- **Evidence**: Classification of evidence type (COMPUTATIONAL/PREDICTED, EXPERIMENTAL)

### Positional Information
- **Residues**: Specific amino acid positions affected
- **Regions**: Continuous sequence regions with start/end coordinates
- **Annotation Values**: Associated measurements or classifications
- **Units**: Measurement units for quantitative annotations

## Example
---
### Sample Searches to Try
Explore these example protein and annotation combinations:

**Mutagenesis Annotations:**
```
UniProt: Q5VSL9
Annotation Type: MUTAGEN
```
*Shows experimental mutagenesis results and predicted effects*

### What Happens When You Search
1. Type "Q5VSL9" in the UniProt field
2. Enter "MUTAGEN" in the annotation type field
3. Click "Fetch Annotations"
4. The tool displays:
   - Summary statistics showing annotation counts
   - Interactive sequence viewer with color-coded annotation regions
   - Detailed annotation information in expandable viewer
   - Positional mapping of all annotations

## Visual Features
---
### Color-Coded Display
- **Position Accuracy**: Precise mapping to sequence positions
- **Type Distinction**: Easy visual separation of annotation categories

## Use Cases
---
This tool is particularly valuable for:
- **Mutation Analysis**: Studying effects of amino acid substitutions
- **Functional Mapping**: Identifying active sites and binding regions
- **Domain Analysis**: Understanding protein architecture and organization
- **Comparative Studies**: Comparing annotations across protein families
- **Drug Target Analysis**: Identifying potential therapeutic intervention sites
- **Research Planning**: Identifying regions of interest for further study
- **Education**: Teaching protein structure-function relationships
- **Quality Assessment**: Evaluating prediction confidence and experimental evidence
