import React, { useState, useEffect } from "react";
import { motion } from "framer-motion";
import { fadeInUpVariantStatic } from "../../../../components/animations/framerAnim";
import GlassyContainer from "../../../../components/glassy_container/gc";

import DataViewer from "../../../../components/UI/DataViewer/DataViewer";
import { computeMolecularDescriptors } from "../../../../services/api/mcpToolsService";
import { FaCheckDouble } from "react-icons/fa6";
import { IoWarningOutline } from "react-icons/io5";
import { call_endpoint_async } from "../../../../endpoints/caller";
import { endpoints } from "../../../../endpoints/endpoints";

const ComputeMolecularDescriptor = ({ toolData = null, initialSmiles = "" }) => {
    const [smiles, setSmiles] = useState(initialSmiles);
    const [isValidSmiles, setIsValidSmiles] = useState(false);
    const [isLoading, setIsLoading] = useState(false);
    const [results, setResults] = useState(toolData);
    const [error, setError] = useState(null);

    // Handle initial results from props
    useEffect(() => {
        if (toolData) {
            setResults(toolData);
        }
    }, [toolData]);

    // Handle initial SMILES from props
    useEffect(() => {
        if (initialSmiles) {
            setSmiles(initialSmiles);
        }
    }, [initialSmiles]);

    // SMILES validation using the same endpoint as MolExplorer
    useEffect(() => {
        if (smiles.length > 0) {
            // Check if the SMILES is valid
            const payload = {
                type: "MOL",
                value: smiles
            }
            call_endpoint_async(endpoints.validate, payload).then((response) => {
                if (response.data.status === "success") {
                    setIsValidSmiles(response.data.valid);
                }
            }).catch((error) => {
                console.log('SMILES validation error>>', error);
                setIsValidSmiles(false);
            })
        }
        else {
            setIsValidSmiles(false);
        }
    }, [smiles]);

    const handleFormSubmit = async () => {
        if (!smiles.trim()) {
            setError('Please enter a SMILES string.');
            return;
        }
        
        if (!isValidSmiles) {
            setError('Please enter a valid SMILES string.');
            return;
        }
        
        setIsLoading(true);
        setError(null);
        setResults(null);
        
        try {
            const response = await computeMolecularDescriptors({ smiles: smiles.trim() });
            if (response && response.status === "success") {
                var processedResults = response.result["0"];
                processedResults = JSON.parse(processedResults);
                setResults(processedResults.data);
            } else {
                setError(response.message || 'Failed to compute molecular descriptors');
            }
        } catch (err) {
            console.error('Error computing molecular descriptors:', err);
            setError(err.message || 'Failed to compute molecular descriptors');
        } finally {
            setIsLoading(false);
        }
    };

    const handleReset = () => {
        setResults(null);
        setSmiles("");
        setError(null);
        setIsValidSmiles(false);
    };

    return (
        <motion.div
            variants={fadeInUpVariantStatic}
            initial="initial"
            animate="animate"
            className="utility-component-container"
        >
            <div className="utility-input-section">
                <motion.div
                    initial={{ opacity: 0, y: 40 }}
                    animate={{ opacity: 1, y: 0 }}
                    transition={{ duration: 0.5, ease: "easeOut" }}
                >
                    <GlassyContainer>
                    <h3 style={{ marginBottom: '15px', fontWeight: '700' }}>Molecular Descriptor Calculator</h3>
                    <p style={{ marginBottom: '15px', color: 'var(--color-text-secondary)' }}>
                        Calculate molecular descriptors and properties from SMILES structures.
                    </p>
                    
                    {/* Add inline styles for input placeholder */}
                    <style>
                        {`
                            .smiles-input::placeholder {
                                color: var(--color-text-secondary);
                                font-weight: 400;
                            }
                            .smiles-input:focus {
                                outline: none;
                                background-color: var(--glassy-color);
                            }
                        `}
                    </style>
                    
                    {/* Custom SMILES Input with Validation */}
                    <div style={{ display: 'flex', flexDirection: 'column', gap: '10px' }}>
                        <div style={{ display: 'flex', alignItems: 'center', gap: '10px', width: '100%' }}>
                            <input
                                type="text"
                                className="smiles-input"
                                value={smiles}
                                onChange={(e) => setSmiles(e.target.value)}
                                onKeyDown={(e) => {
                                    if (e.key === 'Enter' && isValidSmiles && !isLoading) {
                                        handleFormSubmit();
                                    }
                                }}
                                placeholder="Enter SMILES string (e.g., CC(C)(C)OC(=O)NCc1ccccc1)"
                                disabled={isLoading}
                                style={{
                                    flexGrow: 1,
                                    padding: '10px',
                                    height: 'auto',
                                    borderRadius: '15px',
                                    border: 'none',
                                    outline: 'none',
                                    backgroundColor: 'var(--glassy-color)',
                                    color: 'var(--color-text-primary)',
                                    fontSize: '16px',
                                    fontWeight: '600',
                                    minHeight: '30px',
                                    transition: 'background-color 0.2s ease'
                                }}
                                onFocus={(e) => {
                                    e.target.style.backgroundColor = 'var(--glassy-color)';
                                }}
                                onBlur={(e) => {
                                    e.target.style.backgroundColor = 'var(--glassy-color)';
                                }}
                            />
                            {/* Validation Icon */}
                            {isValidSmiles ? (
                                <FaCheckDouble style={{ fontSize: '2rem', color: 'var(--color-success)' }} />
                            ) : (
                                <IoWarningOutline style={{ 
                                    fontSize: '2rem', 
                                    color: smiles.length > 0 ? 'var(--color-alert)' : 'var(--glassy-color)' 
                                }} />
                            )}
                            <button
                                onClick={handleFormSubmit}
                                disabled={!isValidSmiles || isLoading}
                                style={{
                                    padding: '10px 16px',
                                    borderRadius: '15px',
                                    border: 'none',
                                    backgroundColor: (!isValidSmiles || isLoading) ? 'var(--color-disabled, #6c757d)' : 'var(--color-success, #28a745)',
                                    color: '#fff',
                                    cursor: (!isValidSmiles || isLoading) ? 'not-allowed' : 'pointer',
                                    fontSize: '1em',
                                    fontWeight: '500',
                                    minHeight: '50px',
                                    whiteSpace: 'nowrap'
                                }}
                            >
                                {isLoading ? 'Computing...' : 'Compute Descriptors'}
                            </button>
                        </div>
                        {/* Validation message */}
                        {smiles.length > 0 && !isValidSmiles && (
                            <p style={{ 
                                color: 'var(--color-alert, #ff6b6b)', 
                                fontSize: '0.9em', 
                                margin: '5px 0 0 0' 
                            }}>
                                Please enter a valid SMILES string
                            </p>
                        )}
                    </div>
                    {(results || toolData) && (
                        <button 
                            onClick={handleReset}
                            style={{
                                marginTop: '10px',
                                padding: '8px 16px',
                                backgroundColor: 'var(--color-secondary)',
                                color: 'var(--color-text-primary)',
                                border: 'none',
                                borderRadius: '4px',
                                cursor: 'pointer',
                                transition: 'all 0.2s ease'
                            }}
                        >
                            {toolData ? "Reset to Initial" : "Reset"}
                        </button>
                    )}
                </GlassyContainer>
                </motion.div>
            </div>

            {error && (
                <motion.div 
                    className="utility-results-section" 
                    style={{ marginTop: '20px' }}
                    initial={{ opacity: 0, y: 50 }}
                    animate={{ opacity: 1, y: 0 }}
                    exit={{ opacity: 0, y: 50 }}
                    transition={{ duration: 0.6, ease: "easeOut" }}
                >
                    <GlassyContainer>
                        <div style={{
                            padding: '16px',
                            backgroundColor: 'rgba(255, 71, 87, 0.1)',
                            borderRadius: '8px',
                            border: '1px solid rgba(255, 71, 87, 0.3)'
                        }}>
                            <h4 style={{ marginBottom: '8px', color: 'var(--color-error, #ff6b6b)', margin: '0 0 8px 0' }}>
                                ✗ Error
                            </h4>
                            <p style={{ color: 'var(--color-error, #ff6b6b)', margin: 0 }}>
                                {error}
                            </p>
                        </div>
                    </GlassyContainer>
                </motion.div>
            )}

            {results && (
                <motion.div 
                    className="utility-results-section" 
                    style={{ marginTop: '20px' }}
                    initial={{ opacity: 0, y: 60 }}
                    animate={{ opacity: 1, y: 0 }}
                    exit={{ opacity: 0, y: 60 }}
                    transition={{ duration: 0.7, ease: "easeOut", delay: 0.2 }}
                >
                    <GlassyContainer>
                        <DataViewer 
                            data={results} 
                            title={`Molecular Descriptors for SMILES: ${smiles || initialSmiles}`}
                        />
                    </GlassyContainer>
                </motion.div>
            )}
        </motion.div>
    );
};

export default ComputeMolecularDescriptor;
