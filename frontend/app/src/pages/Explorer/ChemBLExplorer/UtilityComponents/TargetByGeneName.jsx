import React, { useState, useEffect } from "react";
import { motion } from "framer-motion";
import { fadeInUpVariantStatic } from "../../../../components/animations/framerAnim";
import GlassyContainer from "../../../../components/glassy_container/gc";
import SimpleInputBox from "../../../../components/UI/SimpleInputBox/SimpleInputBox";
import { findTargetByGeneName } from "../../../../services/api/mcpToolsService";
import DataViewer from "../../../../components/UI/DataViewer/DataViewer";

const TargetByGeneName = ({ toolData = null, initialGeneName = "" }) => {
    const [geneName, setGeneName] = useState(initialGeneName);
    const [isLoading, setIsLoading] = useState(false);
    const [results, setResults] = useState(toolData);
    const [error, setError] = useState(null);

    // Handle initial results from props
    useEffect(() => {
        if (toolData) {
            setResults(toolData);
        }
    }, [toolData]);

    // Handle initial gene name from props
    useEffect(() => {
        if (initialGeneName) {
            setGeneName(initialGeneName);
        }
    }, [initialGeneName]);

    const handleFormSubmit = async () => {
        if (!geneName.trim()) return;
        
        setIsLoading(true);
        setError(null);
        setResults(null);
        
        try {
            const response = await findTargetByGeneName({ gene_name: geneName.trim() });
            if (response && response.status === "success"){
                var components = response.result["0"]
                components = JSON.parse(components);
                setResults(components.data);
            }            
            
        } catch (err) {
            console.error('Error fetching target data:', err);
            setError(err.message || 'Failed to fetch target data');
        } finally {
            setIsLoading(false);
        }
    };

    const handleReset = () => {
        setResults(null);
        setGeneName("");
        setError(null);
    };

    return (
        <motion.div
            variants={fadeInUpVariantStatic}
            initial="initial"
            animate="animate"
            className="utility-component-container"
        >
            <div className="utility-input-section">
                <GlassyContainer>
                    <h3 style={{ marginBottom: '15px', fontWeight: '700' }}>Target by Gene Name</h3>
                    <p style={{ marginBottom: '15px', color: 'var(--color-text-secondary)' }}>
                        Enter a gene name to find associated biological targets in ChEMBL database.
                    </p>
                    <SimpleInputBox
                        value={geneName}
                        onChange={setGeneName}
                        onSubmit={handleFormSubmit}
                        placeholder="Enter gene name (e.g., EGFR, TP53)"
                        loading={isLoading}
                        disabled={isLoading}
                    />
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
            </div>

            {error && (
                <div className="utility-results-section" style={{ marginTop: '20px' }}>
                    <GlassyContainer>
                        <h4 style={{ marginBottom: '15px', color: 'var(--color-error, #ff6b6b)' }}>
                            Error
                        </h4>
                        <p style={{ color: 'var(--color-error, #ff6b6b)', margin: 0 }}>
                            {error}
                        </p>
                    </GlassyContainer>
                </div>
            )}

            {results && (
                <div className="utility-results-section" style={{ marginTop: '20px' }}>
                    <GlassyContainer>
                        <DataViewer 
                            data={results} 
                            title={`ChEMBL Targets for Gene: ${geneName}`}
                            initiallyExpanded={false}
                            maxDepth={8}
                        />
                    </GlassyContainer>
                </div>
            )}
        </motion.div>
    );
};

export default TargetByGeneName;
