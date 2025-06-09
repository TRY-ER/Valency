import React, { useState, useEffect } from "react";
import { motion } from "framer-motion";
import { fadeInUpVariantStatic } from "../../../../components/animations/framerAnim";
import GlassyContainer from "../../../../components/glassy_container/gc";

import { standardizeMoleculeFromSmiles } from "../../../../services/api/mcpToolsService";
import { FaCheckDouble } from "react-icons/fa6";
import { IoWarningOutline } from "react-icons/io5";
import { call_endpoint_async } from "../../../../endpoints/caller";
import { endpoints } from "../../../../endpoints/endpoints";

// CTABViewer component (inline, similar to SmilesToCTAB)
const CTABViewer = ({ data, title }) => {
    const [activeTab, setActiveTab] = useState('raw');

    if (!data) return null;

    const ctabData = data.standard_molblock || data.ctab || (typeof data === 'string' ? data : JSON.stringify(data, null, 2));

    const copyToClipboard = (text) => {
        const textToCopy = typeof text === 'string' ? text : JSON.stringify(text, null, 2);
        navigator.clipboard.writeText(textToCopy).then(() => {
            // You could add a toast notification here if needed
            console.log('CTAB copied to clipboard');
        }).catch(err => {
            console.error('Failed to copy CTAB: ', err);
        });
    };

    return (
        <div style={{ width: '100%', background: 'var(--color-bg-primary, #ffffff)', borderRadius: '8px', overflow: 'hidden' }}>
            <div style={{ padding: '16px 20px 0 20px', borderBottom: '1px solid var(--c-light-border, #e1e5e9)' }}>
                <h4 style={{ margin: '0 0 16px 0', color: 'var(--color-text-primary, #333333)', fontWeight: '600', fontSize: '1.1rem' }}>
                    {title}
                </h4>
                <div style={{ display: 'flex', gap: '8px', marginBottom: '-1px' }}>
                    <button 
                        onClick={() => setActiveTab('raw')}
                        style={{
                            padding: '8px 16px',
                            background: activeTab === 'raw' ? 'var(--color-bg-primary, #ffffff)' : 'transparent',
                            border: '1px solid var(--c-light-border, #e1e5e9)',
                            borderBottom: 'none',
                            borderRadius: '6px 6px 0 0',
                            cursor: 'pointer',
                            fontSize: '0.9rem',
                            color: activeTab === 'raw' ? 'var(--color-accent, #0066cc)' : 'var(--color-text-secondary, #666666)'
                        }}
                    >
                        Standardized CTAB
                    </button>
                    {(data.smiles || data.molecular_formula) && (
                        <button 
                            onClick={() => setActiveTab('info')}
                            style={{
                                padding: '8px 16px',
                                background: activeTab === 'info' ? 'var(--color-bg-primary, #ffffff)' : 'transparent',
                                border: '1px solid var(--c-light-border, #e1e5e9)',
                                borderBottom: 'none',
                                borderRadius: '6px 6px 0 0',
                                cursor: 'pointer',
                                fontSize: '0.9rem',
                                color: activeTab === 'info' ? 'var(--color-accent, #0066cc)' : 'var(--color-text-secondary, #666666)'
                            }}
                        >
                            Molecule Info
                        </button>
                    )}
                </div>
            </div>

            <div style={{ padding: '20px', minHeight: '200px' }}>
                {activeTab === 'raw' && (
                    <motion.div
                        initial={{ opacity: 0, y: 20 }}
                        animate={{ opacity: 1, y: 0 }}
                        transition={{ duration: 0.4, ease: "easeOut" }}
                    >
                        <div style={{ display: 'flex', justifyContent: 'space-between', alignItems: 'center', marginBottom: '16px', paddingBottom: '12px', borderBottom: '1px solid var(--c-light-border, #e1e5e9)' }}>
                            <h5 style={{ margin: '0', color: 'var(--color-text-primary, #333333)', fontWeight: '600' }}>
                                Standardized CTAB Content
                            </h5>
                            <button 
                                style={{
                                    padding: '8px 16px',
                                    background: 'var(--color-accent, #0066cc)',
                                    color: 'white',
                                    border: 'none',
                                    borderRadius: '6px',
                                    cursor: 'pointer',
                                    fontSize: '0.9rem',
                                    fontWeight: '500'
                                }}
                                onClick={() => copyToClipboard(ctabData)}
                            >
                                Copy CTAB
                            </button>
                        </div>
                        <pre style={{ 
                            background: 'var(--color-bg-secondary, #f8f9fa)', 
                            border: '1px solid var(--c-light-border, #e1e5e9)', 
                            borderRadius: '6px', 
                            padding: '16px', 
                            fontFamily: 'Monaco, Menlo, monospace', 
                            fontSize: '0.85rem', 
                            lineHeight: '1.4', 
                            color: 'var(--color-text-primary, #333333)', 
                            overflow: 'auto', 
                            whiteSpace: 'pre', 
                            maxHeight: '400px' 
                        }}>
                            {ctabData}
                        </pre>
                    </motion.div>
                )}

                {activeTab === 'info' && (data.smiles || data.molecular_formula) && (
                    <motion.div
                        initial={{ opacity: 0, y: 20 }}
                        animate={{ opacity: 1, y: 0 }}
                        transition={{ duration: 0.4, ease: "easeOut" }}
                    >
                        <div style={{ background: 'var(--color-bg-secondary, #f8f9fa)', borderRadius: '8px', padding: '16px', border: '1px solid var(--c-light-border, #e1e5e9)' }}>
                            <h5 style={{ margin: '0 0 12px 0', color: 'var(--color-text-primary, #333333)', fontWeight: '600', fontSize: '1rem' }}>
                                Standardized Molecule Information
                            </h5>
                            <div style={{ background: 'var(--color-bg-primary, #ffffff)', borderRadius: '6px', padding: '12px' }}>
                                <div style={{ display: 'grid', gridTemplateColumns: '1fr', gap: '12px' }}>
                                    {data.molecular_formula && (
                                        <div style={{ display: 'flex', justifyContent: 'space-between', alignItems: 'center', padding: '12px', background: 'var(--color-bg-secondary, #f8f9fa)', borderRadius: '6px', border: '1px solid var(--c-light-border, #e1e5e9)' }}>
                                            <span style={{ fontWeight: '500', color: 'var(--color-text-secondary, #666666)' }}>
                                                Molecular Formula:
                                            </span>
                                            <span style={{ 
                                                fontWeight: '600', 
                                                color: 'var(--color-text-primary, #333333)', 
                                                fontFamily: 'Monaco, Menlo, monospace', 
                                                background: 'var(--color-bg-primary, #ffffff)', 
                                                padding: '4px 8px', 
                                                borderRadius: '4px', 
                                                border: '1px solid var(--c-light-border, #e1e5e9)' 
                                            }}>
                                                {data.molecular_formula}
                                            </span>
                                        </div>
                                    )}
                                    {data.smiles && (
                                        <div style={{ display: 'flex', justifyContent: 'space-between', alignItems: 'center', padding: '12px', background: 'var(--color-bg-secondary, #f8f9fa)', borderRadius: '6px', border: '1px solid var(--c-light-border, #e1e5e9)' }}>
                                            <span style={{ fontWeight: '500', color: 'var(--color-text-secondary, #666666)' }}>
                                                Standardized SMILES:
                                            </span>
                                            <span style={{ 
                                                fontWeight: '600', 
                                                color: 'var(--color-text-primary, #333333)', 
                                                fontFamily: 'Monaco, Menlo, monospace', 
                                                background: 'var(--color-bg-primary, #ffffff)', 
                                                padding: '4px 8px', 
                                                borderRadius: '4px', 
                                                border: '1px solid var(--c-light-border, #e1e5e9)' 
                                            }}>
                                                {data.smiles}
                                            </span>
                                        </div>
                                    )}
                                </div>
                            </div>
                        </div>
                    </motion.div>
                )}
            </div>
        </div>
    );
};

const StandardizeMolecules = ({ toolData = null, initialSmiles = "" }) => {
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
            const response = await standardizeMoleculeFromSmiles({ smiles: smiles.trim() });
            if (response && response.status === "success") {
                var processedResults = response.result["0"];
                processedResults = JSON.parse(processedResults);
                
                // Extract the first element from the list which contains the standardized data
                if (processedResults.data && Array.isArray(processedResults.data) && processedResults.data.length > 0) {
                    setResults(processedResults.data[0]);
                } else {
                    setError('No standardization data received from the server');
                }
            } else {
                setError(response.message || 'Failed to standardize molecule');
            }
        } catch (err) {
            console.error('Error standardizing molecule:', err);
            setError(err.message || 'Failed to standardize molecule');
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
            <motion.div 
                className="utility-input-section"
                initial={{ opacity: 0, y: 20 }}
                animate={{ opacity: 1, y: 0 }}
                transition={{ duration: 0.5, ease: "easeOut" }}
            >
                <GlassyContainer>
                    <h3 style={{ marginBottom: '15px', fontWeight: '700' }}>Molecule Standardizer</h3>
                    <p style={{ marginBottom: '15px', color: 'var(--color-text-secondary)' }}>
                        Standardize and normalize molecular structures from SMILES strings using ChEMBL protocols.
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
                                {isLoading ? 'Standardizing...' : 'Standardize Molecule'}
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

            {error && (
                <motion.div 
                    className="utility-results-section" 
                    style={{ marginTop: '20px' }}
                    initial={{ opacity: 0, y: 30 }}
                    animate={{ opacity: 1, y: 0 }}
                    transition={{ duration: 0.5, ease: "easeOut" }}
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
                    initial={{ opacity: 0, y: 30 }}
                    animate={{ opacity: 1, y: 0 }}
                    transition={{ duration: 0.6, ease: "easeOut" }}
                >
                    <GlassyContainer>
                        <CTABViewer 
                            data={results[0]?.standard_molblock || results.standard_molblock} 
                            title={`Standardized Molecule for SMILES: ${smiles || initialSmiles}`}
                        />
                    </GlassyContainer>
                </motion.div>
            )}
        </motion.div>
    );
};

export default StandardizeMolecules;
