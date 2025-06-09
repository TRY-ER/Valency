import React, { useState, useEffect } from "react";
import { motion } from "framer-motion";
import { fadeInUpVariantStatic } from "../../../../components/animations/framerAnim";
import GlassyContainer from "../../../../components/glassy_container/gc";
import DataViewer from "../../../../components/UI/DataViewer/DataViewer";
import { getParentMoleculeFromSmiles } from "../../../../services/api/mcpToolsService";
import { FaCheckDouble } from "react-icons/fa6";
import { IoWarningOutline } from "react-icons/io5";
import { call_endpoint_async } from "../../../../endpoints/caller";
import { endpoints } from "../../../../endpoints/endpoints";

const ParentMolecule = ({ toolData = null, initialSmiles = "" }) => {
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
            const response = await getParentMoleculeFromSmiles({ smiles: smiles.trim() });
            if (response && response.status === "success") {
                var processedResults = response.result["0"];
                processedResults = JSON.parse(processedResults);
                setResults(processedResults);
            } else {
                setError(response.message || 'Failed to extract parent molecule');
            }
        } catch (err) {
            console.error('Error extracting parent molecule:', err);
            setError(err.message || 'Failed to extract parent molecule');
        } finally {
            setIsLoading(false);
        }
    };

    const handleReset = () => {
        setResults(toolData);
        setSmiles(initialSmiles);
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
                    <h3 style={{ marginBottom: '15px', fontWeight: '700' }}>Parent Molecule Extractor</h3>
                    <p style={{ marginBottom: '15px', color: 'var(--color-text-secondary)' }}>
                        Extract parent molecule structures from SMILES representations by removing salts, solvents, and counterions.
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
                                placeholder="Enter SMILES string (including salts/solvents)"
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
                                    outline: 'none',
                                    backgroundColor: isValidSmiles && !isLoading ? 'var(--color-accent)' : 'var(--c-light-border)',
                                    color: isValidSmiles && !isLoading ? 'white' : 'var(--color-text-secondary)',
                                    cursor: isValidSmiles && !isLoading ? 'pointer' : 'not-allowed',
                                    fontSize: '16px',
                                    fontWeight: '700',
                                    transition: 'all 0.2s ease'
                                }}
                            >
                                {isLoading ? 'Extracting...' : 'Extract Parent'}
                            </button>
                        </div>
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
                            title={`Parent Molecule Extraction for SMILES: ${smiles || initialSmiles}`}
                        />
                    </GlassyContainer>
                </motion.div>
            )}
        </motion.div>
    );
};

export default ParentMolecule;
