import React, { useState, useEffect, useCallback, useRef } from "react";
import PropTypes from "prop-types";
import "./UniProtViewer.css"; // Specific styles for UniProtViewer (now contains shared styles too)
import UniProtInputBox from "../../../components/UI/UniProtInputBox/UniProtInputBox";
import InfoBox from "../../../components/UI/InfoBox/InfoBox";
import ThreeDViewer from "../../../components/UI/ThreeDViewer/ThreeDViewer";
import { motion } from "framer-motion";
import { fadeInUpVariantStatic } from "../../../components/animations/framerAnim";
import GlassyContainer from "../../../components/glassy_container/gc";
import { RcsbFv } from "@rcsb/rcsb-saguaro";
import { FaDownload } from 'react-icons/fa'; // Import the download icon

const UniProtViewer = ({ 
    toolData = null, 
    initialAccessionKey = "", 
    hideInputBox = false 
}) => {
    const [accessionKey, setAccessionKey] = useState(initialAccessionKey);
    const [apiData, setApiData] = useState(toolData);
    const [isLoading, setIsLoading] = useState(false);
    const [error, setError] = useState(null);
    const [pdbUrl, setPdbUrl] = useState(toolData?.pdbUrl || "");
    const rcsbFvInstance = useRef(null);
    const sequenceViewerRef = useRef(null);

    useEffect(() => {
        console.log('tool data recieved >>', toolData)
    }, [toolData])

    useEffect(() => {
        console.log("apiData updated >>", apiData);
    }, [apiData])

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

    // Effect to handle initial data passed as props
    useEffect(() => {
        if (toolData) {
            var tempData = null;
            if (typeof toolData === "object"){
                tempData = toolData["0"];
            }
            else{
                tempData = toolData;
            }
            setApiData(tempData);
            setError(null);
            if (tempData.pdbUrl) {
                setPdbUrl(tempData.pdbUrl);
            } else {
                setError("PDB URL not found in provided data.");
                setPdbUrl("");
            }
        }
    }, [toolData]);

    useEffect(() => {
        if (!accessionKey || accessionKey.trim() === "") {
            // Only clear data if no initial data was provided
            if (!toolData) {
                setApiData(null);
                setPdbUrl("");
                setError(null);
            }
            if (rcsbFvInstance.current) {
                rcsbFvInstance.current.unmount();
                rcsbFvInstance.current = null;
            }
        }
    }, [accessionKey, toolData]);

    useEffect(() => {
        if (apiData && apiData.uniprotSequence && apiData.uniprotSequence.length > 0 && sequenceViewerRef.current) {
            if (rcsbFvInstance.current) {
                rcsbFvInstance.current.unmount();
            }

            const sequence = apiData.uniprotSequence;
            console.log("sequence >>", sequence)
            const boardConfig = {
                range: {
                    min: 1,
                    max: sequence.length
                },
                rowTitleWidth: 170,
                includeAxis: true
            };
            const sequenceTrack = {
                trackId: "sequenceTrack",
                trackHeight: 40,
                trackColor: "#FFFFFF", // Changed to whitish
                displayType: "sequence",
                displayColor: "#333333", // Color for the sequence letters
                rowTitle: "SEQUENCE",
                titleFlagColor: "#FFFFFF", // Set title flag color to white
                trackData: [{
                    begin: 1,
                    label: sequence
                }]
            };

            rcsbFvInstance.current = new RcsbFv({
                boardConfigData: boardConfig,
                rowConfigData: [sequenceTrack],
                elementId: sequenceViewerRef.current.id
            });
        } else {
            if (rcsbFvInstance.current) {
                rcsbFvInstance.current.unmount();
                rcsbFvInstance.current = null;
            }
        }

        return () => {
            if (rcsbFvInstance.current) {
                rcsbFvInstance.current.unmount();
                rcsbFvInstance.current = null;
            }
        };
    }, [apiData]);

    const renderFileInfo = (label, url) => {
        if (!url) return null;

        // Extract file type for display, e.g., "PDB File" from "PDB File (.pdb)"
        // This will correctly handle labels like "PAE Image (.png)" -> "PAE Image"
        const fileType = label.split(" (")[0];

        return (
            <div key={label} className="file-download-item" style={{ textAlign: 'center', marginRight: '20px' }}>
                <a href={url} target="_blank" rel="noopener noreferrer" download className="file-download-link">
                    {/* Icon size reduced from 48px to 32px */}
                    <div className="file-icon-placeholder" style={{ fontSize: '32px', marginBottom: '5px' }}>
                        <FaDownload />
                    </div>
                    <span className="file-name">{fileType}</span>
                </a>
            </div>
        );
    };

    return (
        <>
            <motion.div
                initial="hidden"
                animate="visible"
                variants={fadeInUpVariantStatic}
                className="explore-container"
            >
                {!hideInputBox && (
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
                )}

                {error && <div className="error-message" style={{ color: 'red', marginTop: '10px' }}>{error}</div>}

                {apiData && (
                    <div className="explorer-row-2" style={{ marginTop: '20px' }}>
                        <GlassyContainer style={{ marginRight: '10px' }}
                            className="exp-2-row-1">
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
                            <div className="exp-2-row-2">
                                <ThreeDViewer
                                    activeMol={pdbUrl}
                                    isValidMol={!!pdbUrl}
                                    setIsValidMol={() => { }}
                                    visType={"URL_PDB"}
                                />
                            </div>
                        }
                    </div>
                )}

                {/* {apiData && apiData.uniprotSequence && (
                    <div className="explorer-row-3" style={{ marginTop: '20px' }}>
                        <GlassyContainer>
                            <h3 style={{ marginBottom: '15px', fontWeight: '700'}}>Sequence Viewer</h3>
                            <p style={{ marginBottom: '15px' }}>Scroll on the track to zoom in and view !</p>
                            <div
                                key={accessionKey} // Force re-mount of this div when accessionKey changes
                                id="rcsb-saguaro-sequence-viewer"
                                ref={sequenceViewerRef}
                                style={{ minHeight: '270px', backgroundColor: 'gray', borderRadius: '10px', padding: '10px', zIndex: 100 }}
                            ></div>
                        </GlassyContainer>
                    </div>
                )} */}

                {apiData && (
                    <div className="explorer-row-4" style={{ marginTop: '20px' }}>
                        <GlassyContainer className="exp-4-col-1">
                            <h3 style={{ marginBottom: '35px', fontWeight: 'bold' }}>Associated Files</h3>
                            <div style={{ display: 'flex', flexDirection: 'row', justifyContent: 'flex-start', flexWrap: 'wrap' }}>
                                {renderFileInfo("PDB File (.pdb)", apiData.pdbUrl)}
                                {renderFileInfo("CIF File (.cif)", apiData.cifUrl)}
                                {renderFileInfo("BCIF File (.bcif)", apiData.bcifUrl)}
                                {renderFileInfo("PAE Image (.png)", apiData.paeImageUrl)}
                                {renderFileInfo("PAE Data (.json)", apiData.paeDocUrl)}
                            </div>
                        </GlassyContainer>
                    </div>
                )}
            </motion.div>
        </>
    );
};

// PropTypes for better type checking and documentation
UniProtViewer.propTypes = {
    toolData: PropTypes.shape({
        gene: PropTypes.string,
        uniprotDescription: PropTypes.string,
        organismScientificName: PropTypes.string,
        uniprotId: PropTypes.string,
        entryId: PropTypes.string,
        sequenceVersionDate: PropTypes.string,
        modelCreatedDate: PropTypes.string,
        latestVersion: PropTypes.string,
        taxId: PropTypes.string,
        uniprotSequence: PropTypes.string,
        pdbUrl: PropTypes.string,
        cifUrl: PropTypes.string,
        bcifUrl: PropTypes.string,
        paeImageUrl: PropTypes.string,
        paeDocUrl: PropTypes.string,
    }),
    initialAccessionKey: PropTypes.string,
    hideInputBox: PropTypes.bool,
};

export default UniProtViewer;
