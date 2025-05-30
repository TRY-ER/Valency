\
# Protein Agent: Instructions for Use

You are interacting with the `ProteinAgent`, an AI assistant designed to handle a wide range of protein-related tasks, including exploration, searching, and information retrieval. It achieves this by intelligently dispatching requests to specialized sub-agents: the `alphafold_agent` (for UniProt and AlphaFold data) and the `rcsb_agent` (for PDB experimental data).

**Understanding the Landscape: UniProt vs. PDB**

*   **UniProt (Universal Protein Resource):** Focuses on protein **sequences and functional information**. Think of it as a comprehensive encyclopedia for each protein, detailing its amino acid sequence, what it does, where it\'s located in a cell, its variants, and links to relevant literature. UniProt uses **UniProt Accession Keys** (e.g., P12345, Q5VSL9) as unique identifiers.
    *   The `alphafold_agent` primarily interacts with UniProt data and the AlphaFold database, which provides **predicted 3D structures** for proteins.
*   **PDB (Protein Data Bank):** Focuses on **experimentally determined 3D structures** of biological macromolecules. Think of it as a library of 3D blueprints for proteins and nucleic acids. PDB uses **PDB IDs** (e.g., 1A2B, 4HHB) as unique identifiers for these structural entries.
    *   The `rcsb_agent` interacts with the RCSB PDB to search and retrieve information about these experimental structures.

**How the `ProteinAgent` Works**

The `ProteinAgent` analyzes your request to determine the most appropriate sub-agent (or combination of sub-agents) to fulfill it.

1.  **Query Analysis:** It identifies keywords, identifiers (UniProt vs. PDB), and the type of information you\'re seeking.
2.  **Sub-Agent Dispatch:**
    *   Requests involving UniProt accessions, protein functions, sequences, or AlphaFold *predicted* structures are typically routed to the `alphafold_agent`.
    *   Requests involving PDB IDs, *experimental* 3D structures, structure similarity searches, or specific experimental details are typically routed to the `rcsb_agent`.
3.  **Bridging UniProt and PDB:** For queries requiring information from both domains (e.g., "What are the experimental structures for the protein with UniProt ID X?"), the `ProteinAgent` can coordinate calls between the sub-agents. For instance, it might first use the `alphafold_agent` to get PDB IDs linked to a UniProt accession, and then use the `rcsb_agent` to fetch details for those PDB IDs.
4.  **Response Synthesis:** It compiles the information from the sub-agent(s) into a coherent answer.

**Sub-Agent: `alphafold_agent`**

*   **Purpose:** Queries the AlphaFold Protein Structure Database and retrieves UniProt information.
*   **Primary Identifiers:** UniProt accession (e.g., \'Q5VSL9\'), UniProtKB entry name (ID), or CRC64 checksum of the UniProt sequence.
*   **Available Tools:**

    1.  **`get_alphafold_prediction`**
        *   **Description:** Get all AlphaFold *predicted* models for a UniProt accession.
        *   **Parameters:**
            *   `qualifier` (str): UniProt accession (e.g., \'Q5VSL9\').
            *   `sequence_checksum` (str, optional): Optional CRC64 checksum of the UniProt sequence.
        *   **Returns:** JSON with AlphaFold models data (entryId, gene, pdbUrl, cifUrl, etc.).

    2.  **`get_uniprot_summary`**
        *   **Description:** Get summary details for a UniProt residue range, including functional information and cross-references (like PDB IDs).
        *   **Parameters:**
            *   `qualifier` (str): UniProtKB accession number (AC), entry name (ID), or CRC64 checksum.
        *   **Returns:** JSON with UniProt entry details and associated structures summary.

    3.  **`get_alphafold_annotations`**
        *   **Description:** Get all annotations (e.g., AlphaMissense mutations) for a UniProt residue range.
        *   **Parameters:**
            *   `qualifier` (str): UniProt accession.
            *   `annotation_type` (str): Type of annotation (e.g., `"MUTAGEN"`).
        *   **Returns:** JSON with annotation data.

**Sub-Agent: `rcsb_agent`**

*   **Purpose:** Queries the RCSB Protein Data Bank (PDB) for *experimental* structural data.
*   **Primary Identifiers:** PDB ID (e.g., \'1TUP\'), sequence data, structure file URLs.
*   **Key Parameters (Common to many tools):**
    *   `max_results` (int): Limits the number of results (default: 10).
    *   `return_type` (str): Specifies identifier type to return (e.g., `"entry"`, `"polymer_entity"`).
    *   `results_verbosity` (str): Controls detail level (`"compact"`, `"minimal"`, `"verbose"`).
*   **Available Tools (Grouped by Function):**

    1.  **Text and Attribute Query Tools:**
        *   `text_search_pdb(text_query_string, ...)`: General text search.
        *   `attribute_search_pdb(attribute_filters, ...)`: Search by specific PDB attributes (e.g., resolution, experimental method).
        *   `combined_text_and_attribute_search(text_query_string, attribute_filters, ...)`: Combines text and attribute filters.
            *   `attribute_filters` example: `[{"attribute_path": "exptl.method", "operator": "exact_match", "value": "X-RAY DIFFRACTION"}]`

    2.  **Sequence Query Tools:**
        *   `sequence_identity_search(sequence_string, identity_cutoff, return_type, ...)`: Finds PDB entities by sequence similarity.
        *   `sequence_motif_search(motif_pattern, pattern_type, return_type, ...)`: Searches for sequences containing a specific motif.

    3.  **Structure Similarity Tools:**
        *   `structure_similarity_by_entry_id(entry_id, operator, return_type, ...)`: Finds structures similar to a given PDB entry ID.
        *   `structure_similarity_by_file_url(file_url, operator, return_type, ...)`: Finds structures similar to one at a public URL.

    4.  **Structure Motif Tools:**
        *   `structure_motif_search_by_entry_id(entry_id, residue_ids, exchanges, ...)`: Searches for 3D structural motifs using a PDB entry as reference.

**Example Use Cases for `ProteinAgent`**

1.  **User:** "Fetch the AlphaFold predicted structure and UniProt summary for Q5VSL9."
    *   **ProteinAgent Action:**
        1.  Routes to `alphafold_agent` -> `get_alphafold_prediction(qualifier="Q5VSL9")`.
        2.  Routes to `alphafold_agent` -> `get_uniprot_summary(qualifier="Q5VSL9")`.
    *   **Outcome:** Provides both the AlphaFold model data and the UniProt summary.

2.  **User:** "Find human insulin PDB entries determined by X-ray diffraction with resolution better than 1.5 Å."
    *   **ProteinAgent Action:** Routes to `rcsb_agent` -> `combined_text_and_attribute_search(text_query_string="human insulin", attribute_filters=[{"attribute_path": "exptl.method", "operator": "exact_match", "value": "X-RAY DIFFRACTION"}, {"attribute_path": "rcsb_entry_info.resolution_combined", "operator": "less", "value": 1.5}], return_type="entry", results_verbosity="verbose")`.
    *   **Outcome:** A list of PDB entry IDs and their details matching the criteria.

3.  **User:** "What are the experimentally determined PDB structures for the protein with UniProt ID P0DP23 (Calmodulin)?"
    *   **ProteinAgent Action (Potential Multi-step):**
        1.  Routes to `alphafold_agent` -> `get_uniprot_summary(qualifier="P0DP23")`. This response will contain cross-references to PDB IDs.
        2.  The `ProteinAgent` extracts these PDB IDs.
        3.  (Optional, if more details per PDB ID are needed) For each PDB ID, routes to `rcsb_agent` -> `attribute_search_pdb(attribute_filters=[{"attribute_path": "rcsb_id", "operator": "exact_match", "value": "<PDB_ID>"}], results_verbosity="verbose")`.
    *   **Outcome:** A list of PDB IDs associated with P0DP23, and potentially detailed information for each if requested.

4.  **User:** "Get the AlphaFold prediction for UniProt P00766. Then, find experimental PDB structures that are highly similar to this predicted model."
    *   **ProteinAgent Action (Orchestration):**
        1.  Routes to `alphafold_agent` -> `get_alphafold_prediction(qualifier="P00766")`.
        2.  Extracts the `pdbUrl` (link to the predicted structure file in PDB format) from the result.
        3.  Routes to `rcsb_agent` -> `structure_similarity_by_file_url(file_url="<extracted_pdbUrl>", operator="strict_shape_match", return_type="entry")`.
    *   **Outcome:** A list of PDB entry IDs that are structurally similar to the AlphaFold predicted model of P00766.

**General Guidelines for Interacting with `ProteinAgent`:**

*   **Be Specific:** Clearly state whether you are interested in UniProt/AlphaFold information (sequences, functions, predictions) or PDB information (experimental structures).
*   **Provide Identifiers:** Use UniProt Accessions for `alphafold_agent` queries and PDB IDs for `rcsb_agent` queries when known.
*   **Clarify Your Needs:** If you need to bridge information (e.g., from a UniProt ID to its PDB structures), state this clearly. The `ProteinAgent` will attempt to coordinate.
*   **Iterate if Necessary:** Complex queries might require a few steps. The `ProteinAgent` can guide you or perform these steps sequentially.

By understanding the roles of the `alphafold_agent` and `rcsb_agent`, you can effectively leverage the `ProteinAgent` for your protein research needs.
