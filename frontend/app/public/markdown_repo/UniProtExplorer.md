# UniProt Explorer

The UniProt Explorer is a comprehensive tool for searching and retrieving protein information from the UniProt database. It provides an intuitive interface for accessing the UniProtKB (Universal Protein Knowledgebase) using various search methods and identifiers.

## Features

### Search Types

#### 1. Search UniProtKB
- **Query-based search**: Use natural language queries to search the UniProtKB database
- **Advanced filtering**: Filter results by various criteria
- **Flexible result formats**: Choose from JSON, TSV, FASTA, XML, TXT, List, GFF, OBO, RDF, and XLSX formats
- **Custom fields**: Specify which data fields to retrieve (e.g., 'id,xref_pdb,gene_names')
- **Pagination support**: Control the number of results and navigate through pages
- **Isoform inclusion**: Option to include protein isoforms in search results

#### 2. Get UniProtKB Entry
- **Direct retrieval**: Fetch specific protein entries using UniProt IDs
- **Multiple ID formats**: Support for accession numbers (e.g., P12345) and entry names (e.g., SPIKE_SARS2)
- **Various output formats**: JSON, FASTA, TXT, XML, RDF, GFF formats available

## Usage

### Search UniProtKB

1. **Select Search Type**: Choose "Search UniProtKB" from the dropdown
2. **Enter Query**: Input your search query in the query string field
   - Examples:
     - `gene:BRCA1` - Search for BRCA1 gene
     - `organism:human` - Search for human proteins
     - `protein kinase` - General text search
     - `reviewed:true AND organism:9606` - Reviewed human proteins
3. **Configure Options**:
   - **Result Format**: Choose your preferred output format
   - **Fields**: Specify which data columns to retrieve (optional)
   - **Number of Results**: Set how many results to return (max 500)
   - **Include Isoforms**: Check to include protein isoforms
4. **Execute Search**: Click the "Search" button to retrieve results

### Get UniProtKB Entry

1. **Select Search Type**: Choose "Get UniProtKB Entry" from the dropdown
2. **Enter UniProt ID**: Input the specific UniProt identifier
   - Examples:
     - `P12345` (accession number)
     - `SPIKE_SARS2` (entry name)
3. **Choose Format**: Select the desired output format
4. **Retrieve Entry**: Click "Search" to fetch the protein entry

## Query Examples

### Common Search Patterns

- **Gene-based search**: `gene:TP53`
- **Organism-specific**: `organism:"Homo sapiens"`
- **Protein family**: `family:"protein kinase"`
- **Disease-related**: `annotation:(type:disease diabetes)`
- **Subcellular location**: `locations:(location:"membrane")`
- **Combined search**: `gene:BRCA1 AND organism:human AND reviewed:true`

### Advanced Queries

- **Date range**: `created:[2020-01-01 TO 2023-12-31]`
- **Length range**: `length:[100 TO 500]`
- **Mass range**: `mass:[10000 TO 50000]`
- **Taxonomic range**: `taxonomy:mammalia`

## Result Fields

When using the "Fields" option in Search UniProtKB, you can specify which data to retrieve:

### Common Fields
- `id` - UniProt ID
- `gene_names` - Gene names
- `protein_name` - Protein name
- `organism_name` - Organism name
- `length` - Sequence length
- `mass` - Molecular mass
- `xref_pdb` - PDB cross-references
- `cc_function` - Function annotation
- `cc_subcellular_location` - Subcellular location
- `ft_domain` - Domain features

### Example Field Specification
```
id,gene_names,protein_name,organism_name,xref_pdb,length,mass
```

## Output Formats

### JSON Format
- Structured data suitable for programmatic access
- Complete information with nested objects
- Default format for web applications

### FASTA Format
- Protein sequences in FASTA format
- Includes sequence headers with metadata
- Ideal for sequence analysis tools

### TSV/XLSX Format
- Tabular data for spreadsheet applications
- Customizable columns via fields parameter
- Easy to import into analysis tools

### XML Format
- Structured markup with complete annotations
- UniProtKB XML schema compliant
- Suitable for data exchange

## Tips for Effective Searching

1. **Use specific terms**: More specific queries yield better results
2. **Combine filters**: Use Boolean operators (AND, OR, NOT) for complex queries
3. **Check reviewed entries**: Add `reviewed:true` for manually curated entries
4. **Limit organism**: Specify organism to reduce irrelevant results
5. **Use field restrictions**: Specify fields to get only needed data and improve performance

## Error Handling

The tool provides clear error messages for:
- Invalid UniProt IDs
- Malformed query strings
- Network connectivity issues
- API rate limiting
- Invalid parameter combinations

## Integration

The UniProt Explorer integrates with the broader drug discovery toolkit, allowing seamless transition between:
- Protein structure analysis (PDB Explorer)
- Sequence similarity searches
- Functional annotation analysis
- Cross-database references

This tool serves as a crucial component for protein research, drug target identification, and biochemical pathway analysis within the integrated platform.
