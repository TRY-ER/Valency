// Example usage of ToolResponseHandler with ToolMapper
// This file demonstrates how to use the ToolResponseHandler component

import React, { useState } from 'react';
import ToolResponseHandler from './ToolResponseHandler';
import { getAvailableTools, isToolAvailable } from './ToolMapper';

const ToolResponseHandlerExample = () => {
    const [selectedTool, setSelectedTool] = useState('');
    const [toolData, setToolData] = useState(null);
    
    // Get all available tools for the dropdown
    const availableTools = getAvailableTools();
    
    // Example tool data for different tool types
    const exampleToolData = {
        'MoleculeExplorer': {
            smiles: 'CCO',
            molecularWeight: 46.07,
            formula: 'C2H6O'
        },
        'PolymerExplorer': {
            psmiles: '[*]CC[*]',
            monomerWeight: 28.05,
            bondingIndices: [0, 2]
        },
        'ProteinExplorer': {
            pdbId: '1MO8',
            resolution: '2.0',
            organism: 'Homo sapiens'
        },
        'ADMET': {
            smiles: 'CC(C)CC1=CC=C(C=C1)C(C)C(=O)O',
            predictions: {
                absorption: 0.85,
                distribution: 0.72,
                metabolism: 0.68,
                excretion: 0.79,
                toxicity: 0.23
            }
        }
    };
    
    const handleToolChange = (toolName) => {
        setSelectedTool(toolName);
        setToolData(exampleToolData[toolName] || { message: `Example data for ${toolName}` });
    };
    
    return (
        <div style={{ padding: '20px', maxWidth: '1200px', margin: '0 auto' }}>
            <h2>ToolResponseHandler Example</h2>
            
            <div style={{ marginBottom: '20px' }}>
                <label htmlFor="tool-select" style={{ display: 'block', marginBottom: '10px', fontWeight: 'bold' }}>
                    Select a Tool:
                </label>
                <select
                    id="tool-select"
                    value={selectedTool}
                    onChange={(e) => handleToolChange(e.target.value)}
                    style={{
                        padding: '8px 12px',
                        borderRadius: '4px',
                        border: '1px solid #ccc',
                        fontSize: '14px',
                        minWidth: '200px'
                    }}
                >
                    <option value="">-- Select a Tool --</option>
                    {availableTools.map(tool => (
                        <option key={tool} value={tool}>{tool}</option>
                    ))}
                </select>
                
                {selectedTool && (
                    <div style={{ marginTop: '10px', fontSize: '14px' }}>
                        <span style={{ color: isToolAvailable(selectedTool) ? 'green' : 'red' }}>
                            {isToolAvailable(selectedTool) ? '✓ Tool Available' : '✗ Tool Not Available'}
                        </span>
                    </div>
                )}
            </div>
            
            <div style={{ marginBottom: '20px' }}>
                <h3>Available Tools ({availableTools.length}):</h3>
                <div style={{ display: 'flex', flexWrap: 'wrap', gap: '8px' }}>
                    {availableTools.map(tool => (
                        <span
                            key={tool}
                            onClick={() => handleToolChange(tool)}
                            style={{
                                padding: '4px 8px',
                                backgroundColor: selectedTool === tool ? '#007bff' : '#f8f9fa',
                                color: selectedTool === tool ? 'white' : '#333',
                                borderRadius: '12px',
                                fontSize: '12px',
                                cursor: 'pointer',
                                border: '1px solid #dee2e6'
                            }}
                        >
                            {tool}
                        </span>
                    ))}
                </div>
            </div>
            
            <div style={{ border: '1px solid #dee2e6', borderRadius: '8px', padding: '20px' }}>
                <h3>Tool Interface:</h3>
                <ToolResponseHandler 
                    ToolName={selectedTool} 
                    ToolData={toolData}
                />
            </div>
        </div>
    );
};

export default ToolResponseHandlerExample;
