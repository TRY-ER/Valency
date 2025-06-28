# RCSB PDB Explorer Documentation
---

## Overview
---

The **RCSB PDB Explorer** is a comprehensive protein structure database search tool that provides multiple search methods to discover and explore protein structures from the Research Collaboratory for Structural Bioinformatics Protein Data Bank. This powerful interface allows users to search through thousands of protein structures using various criteria including text searches, specific attributes, sequence similarity, and structural motifs.

## Features
---
- **Multiple Search Types**: Six different search methods to find protein structures
- **Text Search**: Natural language search through protein descriptions and metadata
- **PDB ID Search**: Direct lookup of specific protein structures by their PDB identifiers
- **Attribute Search**: Precise filtering using specific database attributes and values
- **Combined Search**: Sophisticated queries combining text and attribute filters
- **Sequence Similarity**: Find proteins with similar amino acid sequences
- **Motif Search**: Discover proteins containing specific sequence patterns
- **Flexible Results**: Multiple output formats and verbosity levels
- **Interactive Interface**: Dynamic form fields that adapt based on search type

### Tool Sections
---
1. **Search Type Selector**
   - Dropdown menu to choose search method
   - Dynamic interface that changes based on selection
   - Clear descriptions for each search type

2. **Dynamic Input Forms**
   - Text fields for queries and identifiers
   - Dropdown selectors for operators and types
   - Multi-filter support for complex searches
   - Sequence input areas for biological sequences

3. **Results Display**
   - Expandable JSON data viewer
   - Structured presentation of search results
   - Multiple verbosity levels (compact, minimal, verbose)
   - Direct access to PDB identifiers and metadata

## Data Source
---
The tool searches the RCSB Protein Data Bank, providing access to:

- **Structure Information**: Detailed metadata about protein structures
- **Experimental Data**: Information about how structures were determined
- **Sequence Data**: Amino acid, DNA, and RNA sequences
- **Attribute Details**: Specific properties like resolution, experimental methods
- **Search Scores**: Relevance rankings for search results
- **Cross-References**: Links to related databases and publications

## Search Types
---

### **Text Search**
Search through protein descriptions, names, and associated text content.

**What you'll enter:**
- **Text Query**: Natural language search terms (e.g., "hemoglobin", "DNA binding protein")

**What you'll get:**
- List of PDB identifiers matching your search terms
- Relevance scores showing how well each result matches
- Option for detailed metadata about each structure

### **PDB ID Search**
Direct lookup of specific protein structures using their unique identifiers.

**What you'll enter:**
- **PDB ID**: Four-character PDB identifier (e.g., "6M0J", "1TIM")

**What you'll get:**
- Detailed information about the specific protein structure
- Complete metadata and experimental details
- Cross-references to related structures

### **Attribute Search**
Precise filtering using specific database attributes and comparison operators.

**What you'll enter:**
- **Attribute Path**: Database field to search (e.g., "exptl.method", "rcsb_entry_info.resolution_combined")
- **Operator**: How to compare values (exact match, greater than, less than, etc.)
- **Value**: What to compare against (e.g., "X-RAY DIFFRACTION", "2.0")

**What you'll get:**
- Structures matching your specific criteria
- Precise filtering based on experimental or structural properties
- Ability to find structures with specific characteristics

### **Combined Text & Attribute Search**
Sophisticated queries that combine text searches with multiple attribute filters.

**What you'll enter:**
- **Main Text Query**: Primary search terms
- **Attribute Filters**: Multiple filters with different criteria
- **Add/Remove Filters**: Dynamic filter management

**What you'll get:**
- Results that match both text criteria and specific attributes
- Complex queries for detailed structure discovery
- Highly targeted search results

### **Sequence Identity Search**
Find proteins with similar amino acid, DNA, or RNA sequences.

**What you'll enter:**
- **Sequence**: Complete biological sequence (protein, DNA, or RNA)
- **Sequence Type**: Protein, DNA, or RNA
- **Identity Cutoff**: Minimum similarity percentage (0.0 to 1.0)
- **E-value Cutoff**: Statistical significance threshold

**What you'll get:**
- Structures with sequences similar to your input
- Similarity scores and statistical measures
- Related proteins from different organisms or variants

### **Sequence Motif Search**
Discover proteins containing specific sequence patterns or functional motifs.

**What you'll enter:**
- **Motif Pattern**: Sequence pattern (PROSITE format, regex, or simple patterns)
- **Pattern Type**: PROSITE, regex, or simple pattern formats
- **Sequence Type**: Protein, DNA, or RNA

**What you'll get:**
- Proteins containing your specified motif
- Functional families sharing common sequence patterns
- Binding sites, active sites, or structural motifs

## Usage
---
1. **Select Search Type**: Choose from the dropdown menu which type of search you want to perform
2. **Fill in Parameters**: Complete the form fields that appear based on your search type selection
3. **Configure Options**: Set any additional parameters like result limits or verbosity levels
4. **Execute Search**: Click "Search PDB" to run your query
5. **Explore Results**: Review the results in the expandable data viewer
6. **Refine Search**: Modify parameters and search again to refine your results

## What You'll See
---
When you use the RCSB PDB Explorer, the interface displays several key sections:

### Search Interface
- **Search Type Dropdown**: Selection menu for different search methods
- **Dynamic Form Fields**: Input areas that change based on your selected search type
- **Parameter Controls**: Options for result formatting and limits
- **Search Button**: Execute button with loading feedback

### Results Display
The results viewer shows different information based on your search:

**Compact Results:**
- Simple list of PDB identifiers
- Quick overview of matching structures
- Minimal data for fast browsing

**Minimal Results:**
- PDB identifiers with relevance scores
- Ranking information showing match quality
- Balance between detail and performance

**Verbose Results:**
- Complete metadata for each result
- Detailed scoring information
- Service-specific match details
- Comprehensive structural information

### Search Examples
**Text Search Results:**
```
['2PGH', '3PEL', '3GOU', '6IHX', '1G08', '1G09', '1G0A', '2QSP', '5C6E', '3CIU']
```

**Scored Results:**
```
[
  {'identifier': '2PGH', 'score': 1.0},
  {'identifier': '3PEL', 'score': 0.9995626684955836},
  {'identifier': '3GOU', 'score': 0.9989625551654878}
]
```

## Advanced Features
---

### Operator Types
Different comparison methods for attribute searches:
- **Exact Match**: Find structures with exactly matching values
- **Greater Than**: Numerical comparisons for resolution, size, etc.
- **Less Than**: Find structures below certain thresholds
- **Contains**: Text matching within fields
- **In List**: Match any value from a provided list

### Sequence Patterns
Support for multiple pattern formats:
- **PROSITE**: Standard protein motif notation
- **Regex**: Regular expression patterns
- **Simple**: Basic wildcard patterns

### Result Types
Different identifier formats returned:
- **Entry**: Complete PDB entry identifiers
- **Polymer Entity**: Specific molecular chain identifiers
- **Assembly**: Biological assembly identifiers

## Example Searches
---

### Find Hemoglobin Structures
```
Search Type: Text Search
Text Query: hemoglobin
```
*Returns structures related to hemoglobin proteins*

### High-Resolution X-ray Structures
```
Search Type: Attribute Search
Attribute Path: exptl.method
Operator: Exact Match
Value: X-RAY DIFFRACTION
```
*Finds structures determined by X-ray crystallography*

### DNA-Binding Proteins with Good Resolution
```
Search Type: Combined Search
Text Query: DNA binding
Attribute Filter 1: rcsb_entry_info.resolution_combined, less_than, 2.0
```
*Combines text search with resolution filtering*

### Similar Sequences
```
Search Type: Sequence Identity Search
Sequence: MVLSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSFPTTKTYFPHF
Sequence Type: Protein
Identity Cutoff: 0.8
```
*Finds proteins with at least 80% sequence similarity*

## Use Cases
---
This tool is particularly valuable for:
- **Structure Research**: Finding specific protein structures for analysis
- **Comparative Studies**: Discovering related structures across different organisms
- **Method Analysis**: Filtering structures by experimental techniques
- **Sequence Analysis**: Finding structural representatives of sequence families
- **Motif Studies**: Exploring functional sequence patterns in structures
- **Quality Assessment**: Filtering by resolution and experimental quality
- **Drug Discovery**: Finding target proteins and binding sites
- **Educational Research**: Exploring structural biology concepts
- **Database Mining**: Large-scale analysis of protein structure data
- **Family Studies**: Discovering evolutionary relationships through structure
