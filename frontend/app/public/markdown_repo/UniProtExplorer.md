# UniProt Explorer Documentation
---

## Overview
---

The **UniProt Explorer** is a comprehensive protein database search and retrieval tool that provides direct access to the UniProtKB (Universal Protein Knowledgebase). This powerful interface allows users to search through millions of protein sequences and annotations using sophisticated query syntax, or retrieve specific protein entries by their identifiers. It serves as a gateway to one of the world's most comprehensive protein information resources.

## Features
---
- **Flexible Search Methods**: Two distinct search approaches for different research needs
- **Advanced Query Syntax**: Sophisticated search capabilities with field-specific queries
- **Database Search**: Comprehensive searches across the entire UniProtKB database
- **Direct Entry Retrieval**: Quick access to specific protein entries by ID
- **Multiple Output Formats**: Various data formats for different analysis needs
- **Isoform Support**: Option to include protein isoforms in search results
- **Pagination Control**: Manage large result sets with customizable page sizes
- **Real-time Results**: Interactive interface with immediate data display

### Tool Sections
---
1. **Search Method Selector**
   - Dropdown to choose between database search or direct entry retrieval
   - Clear descriptions for each search type
   - Dynamic interface adaptation

2. **Search Configuration**
   - Query string input with advanced syntax support
   - UniProt ID input for direct lookups
   - Result size controls and pagination
   - Isoform inclusion options

3. **Results Display**
   - Structured presentation of protein data
   - Expandable JSON data viewer
   - Multiple format support for different analysis needs
   - Direct access to protein annotations and sequences

## Data Source
---
The tool accesses UniProtKB, the world's leading high-quality protein database, providing:

- **Protein Sequences**: Complete amino acid sequences for millions of proteins
- **Functional Annotations**: Detailed information about protein function and structure
- **Cross-References**: Links to hundreds of biological databases
- **Literature Citations**: Connection to scientific publications
- **Evolutionary Information**: Taxonomic and phylogenetic data
- **Structural Data**: Links to 3D structures and structural annotations
- **Disease Associations**: Medical and clinical relevance information

## Search Types
---

### **Search UniProtKB**
Comprehensive database search using powerful query syntax to find proteins based on various criteria.

**What you'll enter:**
- **Query String**: Advanced search terms with field-specific syntax
- **Number of Results**: How many entries to retrieve (1-500)
- **Include Isoforms**: Whether to include protein variants

**Query Examples:**
- Simple terms: `insulin`, `hemoglobin`
- Organism-specific: `human kinase`, `organism_name:"Escherichia coli"`
- Gene names: `gene:BRCA1`, `gene_exact:TP53`
- Protein families: `family:immunoglobulin`
- Disease-related: `disease:cancer`, `keyword:tumor`
- Structural features: `length:[100 TO 500]`, `mass:[50000 TO *]`

**What you'll get:**
- List of matching protein entries
- Complete protein annotations and metadata
- Sequence information and cross-references
- Functional classifications and descriptions

### **Get UniProtKB Entry**
Direct retrieval of specific protein information using UniProt identifiers.

**What you'll enter:**
- **UniProt ID**: Specific protein identifier (accession number or entry name)

**ID Examples:**
- Accession numbers: `P12345`, `Q9Y6R4`
- Entry names: `SPIKE_SARS2`, `INS_HUMAN`
- Swiss-Prot IDs: `P53_HUMAN`, `BRCA1_HUMAN`

**What you'll get:**
- Complete protein entry information
- Detailed functional annotations
- Sequence data and structural information
- Cross-references to related databases

## Advanced Query Syntax
---

### **Basic Search Operations**
- **Simple terms**: `insulin` (searches all text fields)
- **Multiple terms**: `human antigen` or `human AND antigen`
- **Exact phrases**: `"serum albumin"` (use quotes for exact matches)
- **Exclusion**: `human NOT mouse` or `human -mouse`
- **OR searches**: `human OR mouse`
- **Grouping**: `(human OR mouse) AND kinase`
- **Wildcards**: `anti*` (starts with), `*kinase*` (contains)

### **Field-Specific Searches**
The tool supports sophisticated field-specific queries:

**Identification Fields:**
- **Accession**: `accession:P62988`
- **Gene names**: `gene:HPSE` (partial match), `gene_exact:HPSE` (exact match)
- **Protein names**: `protein_name:CD233`
- **Secondary accession**: `sec_acc:P02023`

**Organism and Taxonomy:**
- **Organism name**: `organism_name:"Homo sapiens"`
- **Organism ID**: `organism_id:9606`
- **Taxonomy**: `taxonomy_name:mammal`, `taxonomy_id:40674`
- **Strain**: `strain:wistar`
- **Virus host**: `virus_host_id:10090`

**Functional Classification:**
- **Protein families**: `family:serpin`
- **Keywords**: `keyword:toxin`, `keyword:KW-0800`
- **EC numbers**: `ec:3.2.1.23`
- **Gene Ontology**: `go:0015629`
- **Protein existence**: `existence:3`

**Physical Properties:**
- **Sequence length**: `length:[500 TO 700]`, `length:[50 TO *]`
- **Molecular mass**: `mass:[500000 TO *]`
- **Fragment status**: `fragment:true`

**Database Cross-References:**
- **PDB structures**: `database:pdb`, `xref:pdb-1aut`
- **Pfam domains**: `database:pfam`
- **ChEBI compounds**: `chebi:18420`
- **InChIKey**: `inchikey:WQZGKKKJIJFFOK-GASJEMHNSA-N`

**Dates and Status:**
- **Creation date**: `date_created:[2012-10-01 TO *]`
- **Modification date**: `date_modified:[2019-01-01 TO 2019-03-01]`
- **Review status**: `reviewed:true` (Swiss-Prot entries)
- **Active status**: `active:false` (obsolete entries)

**Literature and Evidence:**
- **Authors**: `lit_author:ashburner`
- **Publication scope**: `scope:mutagenesis`
- **Mass spectrometry**: `cc_mass_spectrometry:maldi`

### **Range Queries**
Use square brackets for numerical ranges:
- `length:[100 TO 500]` (proteins 100-500 amino acids)
- `mass:[50000 TO *]` (proteins over 50 kDa)
- `date_created:[2020-01-01 TO *]` (entries from 2020 onwards)

## Usage
---
1. **Select Search Method**: Choose between "Search UniProtKB" or "Get UniProtKB Entry"
2. **Configure Search Parameters**:
   - For database search: Enter query string and set result preferences
   - For direct retrieval: Enter specific UniProt ID
3. **Set Additional Options**: Configure result size and isoform inclusion
4. **Execute Search**: Click "Search" to retrieve protein information
5. **Explore Results**: Use the data viewer to examine detailed protein information
6. **Refine Search**: Modify parameters to narrow or broaden your results

## What You'll See
---
When you use the UniProt Explorer, the interface displays:

### Search Interface
- **Search Type Dropdown**: Selection between database search and direct entry retrieval
- **Query Input Fields**: Dynamic forms based on selected search method
- **Parameter Controls**: Options for result size and isoform inclusion
- **Search Controls**: Execute and clear buttons with loading feedback

### Results Display
The results show comprehensive protein information:

**Database Search Results:**
- List of matching protein entries
- Entry identifiers and primary names
- Organism information and gene names
- Functional descriptions and classifications
- Cross-references to related databases

**Direct Entry Results:**
- Complete protein entry information
- Detailed sequence and structural data
- Functional annotations and features
- Literature references and evidence
- Cross-references to external databases

## Example Searches
---

### Find Human Kinases
```
Search Type: Search UniProtKB
Query String: human kinase
Number of Results: 50
Include Isoforms: false
```
*Returns human protein kinases with functional annotations*

### Search for BRCA1 Gene Products
```
Search Type: Search UniProtKB
Query String: gene:BRCA1
Number of Results: 30
Include Isoforms: true
```
*Finds all proteins from the BRCA1 gene including variants*

### Get Specific Protein Entry
```
Search Type: Get UniProtKB Entry
UniProt ID: P04637
```
*Retrieves complete information for the p53 tumor suppressor protein*

### Find Large Membrane Proteins
```
Search Type: Search UniProtKB
Query String: length:[500 TO *] AND keyword:membrane
Number of Results: 100
Include Isoforms: false
```
*Searches for large membrane proteins over 500 amino acids*

### Explore COVID-19 Proteins
```
Search Type: Search UniProtKB
Query String: organism_name:"SARS-CoV-2"
Number of Results: 25
Include Isoforms: false
```
*Finds all proteins from the SARS-CoV-2 virus*

## Use Cases
---
This tool is particularly valuable for:
- **Protein Research**: Comprehensive protein database exploration and analysis
- **Gene Function Studies**: Understanding protein products of specific genes
- **Comparative Biology**: Studying protein evolution across species
- **Drug Discovery**: Identifying and characterizing therapeutic targets
- **Structural Biology**: Finding proteins for structural studies
- **Systems Biology**: Mapping protein networks and pathways
- **Disease Research**: Exploring disease-associated proteins
- **Educational Research**: Learning about protein structure and function
- **Literature Mining**: Finding proteins mentioned in scientific publications
- **Database Integration**: Collecting protein data for bioinformatics pipelines
