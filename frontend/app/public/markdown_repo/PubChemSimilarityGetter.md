# PubChem Similarity Getter

## Overview

The **PubChem Similarity Getter** is a molecular similarity search tool that finds compounds with similar chemical structures based on 2D structural fingerprints. Using PubChem's fast similarity search API, this tool employs Tanimoto coefficient calculations to identify structurally related molecules, making it invaluable for drug discovery, chemical analysis, and molecular research.

## Features

### Core Capabilities
- **CID-Based Similarity Search**: Find similar compounds using PubChem Compound Identifiers (CIDs)
- **Adjustable Similarity Threshold**: Control search sensitivity with threshold values from 0-100%
- **Tanimoto Coefficient Calculation**: Uses industry-standard similarity metrics for accurate structural comparison
- **Multiple Results Display**: View and compare multiple similar compounds simultaneously
- **Interactive Molecule Selection**: Click through similar compounds to explore detailed information
- **Structural Visualization**: 2D molecular structure images for visual comparison
- **Comprehensive Compound Data**: Access complete PubChem compound information for each result

### Search Validation
- **Real-time CID Validation**: Instant feedback on CID format correctness
- **Input Sanitization**: Automatic trimming and format checking
- **Error Prevention**: Clear validation messages before search execution

## Tool Sections

### Search Input
- **CID Entry Field**: Input box for entering PubChem Compound Identifiers
- **Search Type Selector**: Currently supports CID-based searches (dropdown for future expansion)
- **Similarity Threshold**: Adjustable parameter for controlling search sensitivity (default: 90%)
- **Search Button**: Initiates the similarity search process

### Results Display
- **Primary Compound View**: Detailed information panel for the selected compound
- **Similar Molecules Grid**: Visual grid showing all similar compounds found
- **Similarity Scores**: Percentage similarity values for each result
- **Structure Images**: 2D molecular diagrams for visual comparison
- **Selection Indicators**: Visual markers showing the currently selected compound

## Data Source

**PubChem Database** - Maintained by the National Center for Biotechnology Information (NCBI)
- **Coverage**: Over 100 million chemical substances and compounds
- **Similarity Algorithm**: Fast 2D similarity search using Tanimoto coefficients
- **Structure Data**: 2D molecular fingerprints and structural representations
- **Update Frequency**: Continuously updated with new submissions and literature data

## Usage

### Basic Search Process
1. **Enter a CID**: Type a valid PubChem Compound Identifier (e.g., "2244" for aspirin)
2. **Set Threshold** (Optional): Adjust similarity threshold if needed (default 90% works well)
3. **Click Search Similar**: Execute the similarity search
4. **Review Results**: Browse through similar compounds in the results grid
5. **Select Compounds**: Click on any compound card to view detailed information

### Input Requirements
- **CID Format**: Positive integers only (e.g., 2244, 5281804, 441243)
- **Threshold Range**: Integer values between 0-100 (higher = more similar compounds only)
- **Connection**: Active internet connection for PubChem API access

## What You'll See

### Search Interface
- Clean input field with placeholder text "Enter PubChem CID (e.g., 2244)"
- Search type dropdown (currently showing "PubChem CID")
- Real-time validation with checkmarks (✓) for valid inputs or warning icons for invalid formats
- "Search Similar" button that activates when input is valid

### During Search
- Loading animation with "🔍 Searching for similar molecules..." message
- Progress indication while the system queries PubChem and retrieves compound data

### Results Display
- **Main Compound Panel**: Comprehensive data viewer showing:
  - Compound name and identifiers
  - Molecular formula and weight
  - SMILES notation
  - Physical and chemical properties
  - Cross-references to other databases

- **Similar Molecules Grid**: Card-based layout showing:
  - Compound CID numbers
  - Similarity percentages (e.g., "94.2% similar")
  - Compound names (truncated if long)
  - 2D structure images
  - Selection indicators for the active compound

### Error Handling
- **Invalid CID**: "Invalid CID format. Must be a positive integer (e.g., 2244, 5281804)."
- **No Results**: "No similar molecules found for CID 'X' with similarity threshold Y%"
- **API Errors**: Clear messages explaining connection or data retrieval issues

## Examples

### Example 1: Aspirin Similarity Search
**Input**: CID `2244` (Aspirin)
**Threshold**: `90%`
**Expected Results**: 
- Salicylic acid derivatives
- Other NSAIDs with similar structures
- Acetylated compounds
- Similarity scores typically ranging from 90-95%

### Example 2: Caffeine Similarity Search
**Input**: CID `2519` (Caffeine)
**Threshold**: `85%`
**Expected Results**:
- Theophylline and related xanthines
- Other methylxanthine derivatives
- Purine-based stimulants
- 15-25 similar compounds with varying similarity scores

### Example 3: High-Threshold Search
**Input**: CID `5281804` (Resveratrol)
**Threshold**: `95%`
**Expected Results**:
- Very few, highly similar stilbene derivatives
- Close structural analogs only
- Fewer results but higher structural similarity

## Use Cases

### Drug Discovery
- **Lead Compound Identification**: Find structurally similar compounds to known active drugs
- **Scaffold Hopping**: Discover alternative chemical scaffolds with similar properties
- **Side Effect Prediction**: Identify compounds with similar structures that might share side effects
- **Patent Analysis**: Search for similar compounds in patent landscape analysis

### Chemical Research
- **Structure-Activity Relationships**: Explore how structural changes affect biological activity
- **Metabolite Prediction**: Find similar compounds that might be metabolites or degradation products
- **Chemical Space Exploration**: Map chemical space around compounds of interest
- **Literature Mining**: Find related compounds mentioned in scientific literature

### Educational Applications
- **Chemical Similarity Concepts**: Demonstrate principles of molecular similarity
- **Database Exploration**: Learn to navigate and utilize chemical databases
- **Comparative Analysis**: Study structural variations within compound families
- **Research Training**: Practice using computational chemistry tools

### Quality Control
- **Impurity Identification**: Find compounds similar to known impurities
- **Reference Standard Selection**: Identify suitable reference compounds for analysis
- **Method Development**: Select similar compounds for analytical method validation
- **Batch Comparison**: Compare compound batches by finding similar structures

---

*This tool provides access to PubChem's extensive chemical database for similarity searching. Results depend on data availability in PubChem and the structural complexity of the query compound. For optimal results, use well-characterized compounds with complete structural information.*
