import React, { useState, useEffect, useCallback } from "react";
// import "./ProtExplorer.css"; // REMOVED
import "./UniProtViewer.css"; // Specific styles for UniProtViewer (now contains shared styles too)
import UniProtInputBox from "../../../components/UI/UniProtInputBox/UniProtInputBox";
import InfoBox from "../../../components/UI/InfoBox/InfoBox";
import ThreeDViewer from "../../../components/UI/ThreeDViewer/ThreeDViewer";
import { motion } from "framer-motion";
import { fadeInUpVariantStatic } from "../../../components/animations/framerAnim";
import GlassyContainer from "../../../components/glassy_container/gc";

const UniProtViewer = () => {
    const [accessionKey, setAccessionKey] = useState("");
    const [apiData, setApiData] = useState(null);
    const [isLoading, setIsLoading] = useState(false);
    const [error, setError] = useState(null);
    const [pdbUrl, setPdbUrl] = useState("");

    const handleValidationComplete = useCallback((isValid, data, errorMessage) => {
        setIsLoading(false);
        if (isValid && data) {
            setApiData(data);
            setError(null);
            if (data.pdbUrl) {
                setPdbUrl(data.pdbUrl);
            } else {
                setError("PDB URL not found in API response.");
                setPdbUrl("");
            }
        } else {
            setApiData(null);
            setPdbUrl("");
            setError(errorMessage || "Validation failed.");
        }
    }, []);

    useEffect(() => {
        if (!accessionKey || accessionKey.trim() === "") {
            setApiData(null);
            setPdbUrl("");
            setError(null);
        }
    }, [accessionKey]);

    const renderFileInfo = (label, url) => (
        url ? <p key={label}><a href={url} target="_blank" rel="noopener noreferrer">{label}</a></p> : null
    );

    return (
        <>
            <motion.div
                initial="hidden"
                animate="visible"
                variants={fadeInUpVariantStatic}
                className="explore-container"
            >
                <div className="explorer-row-1">
                    <GlassyContainer>
                        <UniProtInputBox
                            value={accessionKey}
                            onChange={setAccessionKey}
                            onValidationComplete={handleValidationComplete}
                            header="Enter UniProt Accession Key"
                            placeholder="Try UniProt Accession: Q5VSL9"
                        />
                    </GlassyContainer>
                </div>

                {error && <div className="error-message" style={{color: 'red', marginTop: '10px'}}>{error}</div>}

                {apiData && (
                    <div className="explorer-row-2" style={{marginTop: '20px'}}>
                        <GlassyContainer style={{marginRight: '10px'}}
                         className="exp-2-row-1"> {/* Changed flex to 2 for 40% width */}
                            <h3 style={{ marginBottom: '15px' }}>Protein Information</h3>
                            {
                                (() => {
                                    const displayInfo = {
                                        "Gene": apiData.gene,
                                        "Description": apiData.uniprotDescription,
                                        "Organism": apiData.organismScientificName,
                                        "UniProt ID": apiData.uniprotId,
                                        "Entry ID": apiData.entryId,
                                        "Sequence Version Date": apiData.sequenceVersionDate,
                                        "Model Created Date": apiData.modelCreatedDate,
                                        "Latest Version": apiData.latestVersion,
                                        "Tax ID": apiData.taxId
                                    };
                                    return Object.entries(displayInfo).map(([key, value]) => (
                                        value ? (
                                            <div key={key} className="info-row">
                                                <span className="info-label">{key}</span>
                                                <span className="info-value">{String(value)}</span>
                                            </div>
                                        ) : null
                                    ));
                                })()
                            }
                        </GlassyContainer>

                        {pdbUrl && 
                            <div className="exp-2-row-2"> {/* Changed flex to 3 for 60% width */}
                                <ThreeDViewer 
                                    activeMol={pdbUrl}
                                    isValidMol={!!pdbUrl} 
                                    setIsValidMol={() => {}} 
                                    visType={"URL_PDB"} 
                                />
                            </div>
                        }
                    </div>
                )}
                
                {apiData && apiData.uniprotSequence && (
                     <div className="explorer-row-3" style={{marginTop: '20px'}}>
                        <GlassyContainer>
                            <h3 style={{ marginBottom: '15px' }}>Protein Sequence</h3>
                            <pre style={{whiteSpace: 'pre-wrap', wordBreak: 'break-all', maxHeight: '200px', overflowY: 'auto'}}>
                                {apiData.uniprotSequence}
                            </pre>
                        </GlassyContainer>
                    </div>
                )}

                {apiData && (
                    <div className="explorer-row-4" style={{marginTop: '20px'}}>
                        <GlassyContainer>
                            <h3 style={{ marginBottom: '15px' }}>Associated Files</h3>
                            {renderFileInfo("PDB File (.pdb)", apiData.pdbUrl)}
                            {renderFileInfo("CIF File (.cif)", apiData.cifUrl)}
                            {renderFileInfo("BCIF File (.bcif)", apiData.bcifUrl)}
                            {renderFileInfo("PAE Image (.png)", apiData.paeImageUrl)}
                            {renderFileInfo("PAE Data (.json)", apiData.paeDocUrl)}
                        </GlassyContainer>
                    </div>
                )}
            </motion.div>
        </>
    );
};

export default UniProtViewer;
