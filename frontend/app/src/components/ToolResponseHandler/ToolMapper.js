
import React from 'react';

// Import Explorer Components
import MolEComponent from '../../pages/Explorer/MolExplorer/MolExplorer';
import PolyEComponent from '../../pages/Explorer/PolyExplorer/PolyExplorer';
import BasicProtViewer from '../../pages/Explorer/ProtExplorer/BasicProtViewer';
import MolStarViewer from '../../pages/Explorer/ProtExplorer/MolStarViewer';

// Import Alphafold Components
import UniProtViewer from '../../pages/Explorer/ProtExplorer/UniProtViewer';
import UniProtSummaryViewer from '../../pages/Explorer/ProtExplorer/UniProtSummaryViewer';
import AlphafoldAnnotationsViewer from '../../pages/Explorer/ProtExplorer/AlphafoldAnnotationsViewer';

// Import ChemBL Components
import PSMILESListComponent from '../../pages/Generator/BRICS/PSMILESListComponent';
import ChemBLGetter from '../../pages/Explorer/ChemBLExplorer/ChemBLGetter';
import ChemBLSimilarityGetter from '../../pages/Explorer/ChemBLExplorer/ChemBLSimilarityGetter';
import ChemBLActivityFetcher from '../../pages/Explorer/ChemBLExplorer/ChemBLActivityFetcher';
import TargetByGeneName from '../../pages/Explorer/ChemBLExplorer/UtilityComponents/TargetByGeneName';
import SmilesToCTAB from '../../pages/Explorer/ChemBLExplorer/UtilityComponents/SmilesToCTAB';
import ComputeMolecularDescriptor from '../../pages/Explorer/ChemBLExplorer/UtilityComponents/ComputeMolecularDescriptor';
import ComputeStructuralAlerts from '../../pages/Explorer/ChemBLExplorer/UtilityComponents/ComputeStructuralAlerts';
import StandardizeMolecules from '../../pages/Explorer/ChemBLExplorer/UtilityComponents/StandardizeMolecules';
import ParentMolecule from '../../pages/Explorer/ChemBLExplorer/UtilityComponents/ParentMolecule';

// Import PubChem Components
import PubChemGetter from '../../pages/Explorer/PubChemExplorer/PubChemGetter';
import CompoundProperties from '../../pages/Explorer/PubChemExplorer/UtilityComponents/CompoundProperties';
import CompoundSynonyms from '../../pages/Explorer/PubChemExplorer/UtilityComponents/CompoundSynonyms';
import IdentitySearch from '../../pages/Explorer/PubChemExplorer/UtilityComponents/IdentitySearch';
import PubChemSimilarityGetter from '../../pages/Explorer/PubChemExplorer/PubChemSimilarityGetter';
import SubstructureSearch from '../../pages/Explorer/PubChemExplorer/UtilityComponents/SubstructureSearch';

// Import RCSB Components
import RCSBPDBExplorer from '../../pages/Explorer/RCSBPDBExplorer/RCSBPDBExplorer';
import RCSBStructureSimilaritySearch from '../../pages/Explorer/RCSBPDBExplorer/RCSBStructureSimilaritySearch'; // Import Generator Components

// Import Uniprot Components
import UniProtExplorer from '../../pages/Explorer/UniProtExplorer/UniProtExplorer';

import BRICSComponent from '../../pages/Generator/BRICS/PSMILESContent';
import LSTMComponent from '../../pages/Generator/LSTM/LSTMContent';
import ADMETComponent from '../../pages/Generator/ADMET/ADMET';


// Import other UI components that might be useful for tools
import DataViewer from '../UI/DataViewer';

/**
 * ToolMapper - Maps tool names to their corresponding React components
 * This object serves as a registry for all available tools in the application
 */
const ToolMapper = {
    "get_admet_prediction": ADMETComponent,
    "get_alphafold_prediction": UniProtViewer,   
    "get_uniprot_summary": UniProtSummaryViewer,
    "get_alphafold_annotations": AlphafoldAnnotationsViewer,
    "get_brics_candidates": PSMILESListComponent,
    "find_molecule_by_pref_name": ChemBLGetter,
    "get_molecule_by_chembl_id": ChemBLGetter,
    "find_molecule_by_synonym": ChemBLGetter,
    "get_molecules_by_chembl_ids": ChemBLGetter,

    "find_similar_molecules_by_smiles": ChemBLSimilarityGetter,
    "find_similar_molecules_by_chembl_id": ChemBLSimilarityGetter,
    
    "get_activities_for_target": ChemBLActivityFetcher,
    "get_activities_for_molecule": ChemBLActivityFetcher,
    "find_target_by_gene_name": TargetByGeneName, 
    "smiles_to_ctab": SmilesToCTAB,
    "compute_molecular_descriptors": ComputeMolecularDescriptor,
    "compute_structural_alerts": ComputeStructuralAlerts,
    "standardize_molecule_from_smiles": StandardizeMolecules,
    "get_parent_molecule_from_smiles": ParentMolecule,

    "get_compound_by_cid": PubChemGetter,
    "get_cids_by_name": PubChemGetter,
    "get_cids_by_smiles": PubChemGetter,
    "get_cids_by_inchikey": PubChemGetter,
    "get_compound_properties": CompoundProperties,
    "get_compound_synonyms_by_cid": CompoundSynonyms,
    "fast_identity_search_by_cid": IdentitySearch,
    "fast_similarity_2d_search_by_cid": PubChemSimilarityGetter,
    "fast_substructure_search_by_smiles": SubstructureSearch, 


    "text_search_pdb": RCSBPDBExplorer,
    "get_protein_details_by_id_pypdb": RCSBPDBExplorer,
    "structure_similarity_by_entry_id":  RCSBStructureSimilaritySearch,

    "search_uniprotkb": UniProtExplorer,
    "get_uniprotkb_entry": UniProtExplorer,

    // Utility Components
    'DataViewer': DataViewer,
    'data_viewer': DataViewer,
    
    // Fallback component for unknown tools
    'UnknownTool': ({ toolName, toolData }) => (
        <div className="unknown-tool-container" style={{ 
            padding: '20px', 
            border: '2px dashed #ccc', 
            borderRadius: '8px',
            textAlign: 'center',
            backgroundColor: '#f9f9f9'
        }}>
            <h3 style={{ color: '#666', marginBottom: '10px' }}>
                Tool Not Found: {toolName}
            </h3>
            <p style={{ color: '#888', marginBottom: '15px' }}>
                The requested tool "{toolName}" is not available in the ToolMapper.
            </p>
            {toolData && (
                <details style={{ marginTop: '15px', textAlign: 'left' }}>
                    <summary style={{ cursor: 'pointer', fontWeight: 'bold' }}>
                        View Raw Tool Data
                    </summary>
                    <pre style={{
                        fontFamily: 'monospace',
                        padding: '10px',
                        backgroundColor: '#f0f0f0',
                        borderRadius: '4px',
                        overflow: 'auto',
                        fontSize: '12px',
                        marginTop: '10px'
                    }}>
                        {JSON.stringify(toolData, null, 2)}
                    </pre>
                </details>
            )}
        </div>
    )
};

/**
 * Helper function to get a component by tool name
 * @param {string} toolName - The name of the tool
 * @returns {React.Component|null} - The corresponding React component
 */
export const getToolComponent = (toolName, toolData) => {
    if (!toolName) return null;
    
    // Try exact match first
    if (ToolMapper[toolName]) {
        const Component = ToolMapper[toolName];
        return Component; 
    }
    
    // Try case-insensitive match
    const normalizedToolName = toolName.toLowerCase();
    const matchingKey = Object.keys(ToolMapper).find(key => 
        key.toLowerCase() === normalizedToolName
    );
    
    if (matchingKey) {
        const Component = ToolMapper[matchingKey];
        return Component;

    }
    
    // Try partial match
    const partialMatchKey = Object.keys(ToolMapper).find(key => 
        key.toLowerCase().includes(normalizedToolName) || 
        normalizedToolName.includes(key.toLowerCase())
    );
    
    if (partialMatchKey) {
        return ToolMapper[partialMatchKey];
    }
    
    // Return unknown tool component as fallback
    return ToolMapper.UnknownTool;
};

/**
 * Helper function to get all available tool names
 * @returns {Array<string>} - Array of all available tool names
 */
export const getAvailableTools = () => {
    return Object.keys(ToolMapper).filter(key => key !== 'UnknownTool');
};

/**
 * Helper function to check if a tool exists
 * @param {string} toolName - The name of the tool to check
 * @returns {boolean} - Whether the tool exists
 */
export const isToolAvailable = (toolName) => {
    return getToolComponent(toolName) !== ToolMapper.UnknownTool;
};

export default ToolMapper;