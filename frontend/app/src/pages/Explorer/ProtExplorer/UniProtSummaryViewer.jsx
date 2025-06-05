import React, { useState, useEffect, useCallback, useRef } from "react";
import PropTypes from "prop-types";
import "./UniProtViewer.css"; // Specific styles for UniProtViewer (now contains shared styles too)
import SimpleInputBox from "../../../components/UI/SimpleInputBox/SimpleInputBox";
import InfoBox from "../../../components/UI/InfoBox/InfoBox";
import ThreeDViewer from "../../../components/UI/ThreeDViewer/ThreeDViewer";
import DataViewer from "../../../components/UI/DataViewer/DataViewer";
import { motion } from "framer-motion";
import { fadeInUpVariantStatic } from "../../../components/animations/framerAnim";
import GlassyContainer from "../../../components/glassy_container/gc";
import { RcsbFv } from "@rcsb/rcsb-saguaro";
import { FaDownload } from 'react-icons/fa'; // Import the download icon
import { getUniprotSummary } from "../../../services/api/mcpToolsService";

const UniProtSummaryViewer = ({ 
    toolData = null, 
    initialAccessionKey = "", 
    hideInputBox = false 
}) => {
    const [accessionKey, setAccessionKey] = useState(initialAccessionKey);
    const [inputValue, setInputValue] = useState(initialAccessionKey);
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

    const handleSubmit = useCallback(async (accession) => {
        if (!accession || accession.trim() === "") {
            setError("Please enter a valid UniProt accession key");
            return;
        }

        try {
            setIsLoading(true);
            setError(null);
            setAccessionKey(accession);
            
            console.log('Calling MCP UniProt summary with accession:', accession);
            const response = await getUniprotSummary({ qualifier: accession.trim() });
            const result = response?.result;
            
            if (!result) {
                console.error('No UniProt summary result returned');
                setError("No UniProt summary data available");
                setApiData(null);
                setPdbUrl("");
                return;
            }

            console.log('UniProt summary result:', result);
            console.log('type of UniProt summary result:', typeof result);
            let actualResults = null;
            
            if (typeof result === "object") {
                // Handle different response formats
                if (result["0"]) {
                    // If result is wrapped in an object with "0" key
                    const result_data_json = JSON.parse(result["0"]);
                    if (result_data_json.error === null) {
                        actualResults = result_data_json.data;
                    }
                    else{
                        actualResults = result_data_json;
                    }
                } else {
                    // Direct result object
                    actualResults = result;
                }
            }

            if (actualResults) {
                setApiData(actualResults);
                setError(null);
                if (actualResults.pdbUrl) {
                    setPdbUrl(actualResults.pdbUrl);
                } else {
                    console.log("PDB URL not found in UniProt summary response");
                    setPdbUrl("");
                }
            } else {
                setError("Failed to process UniProt summary data");
                setApiData(null);
                setPdbUrl("");
            }
        } catch (error) {
            console.error('Error fetching UniProt summary:', error);
            setError(`Failed to fetch UniProt summary: ${error.message}`);
            setApiData(null);
            setPdbUrl("");
        } finally {
            setIsLoading(false);
        }
    }, []);

    // Effect to handle initial data passed as props
    useEffect(() => {
        if (toolData) {
            var tempData = null;
            if (typeof toolData === "object"){
                tempData = toolData;
            }
            setApiData(tempData);
            setError(null);
            // if (tempData.pdbUrl) {
            //     setPdbUrl(tempData.pdbUrl);
            // } else {
            //     setError("PDB URL not found in provided data.");
            //     setPdbUrl("");
            // }
        }
    }, [toolData]);

    // Effect to handle initial accession key
    useEffect(() => {
        if (initialAccessionKey && initialAccessionKey.trim() !== "") {
            setInputValue(initialAccessionKey);
            setAccessionKey(initialAccessionKey);
            // Optionally auto-submit if we have initial data
            if (!toolData) {
                handleSubmit(initialAccessionKey);
            }
        }
    }, [initialAccessionKey, toolData, handleSubmit]);

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
                            <SimpleInputBox
                                value={inputValue}
                                onChange={setInputValue}
                                onSubmit={handleSubmit}
                                header="Enter UniProt Accession Key"
                                placeholder="Try UniProt Accession: Q5VSL9"
                                buttonText="Fetch Summary"
                                isLoading={isLoading}
                                error={error}
                            />
                        </GlassyContainer>
                    </div>
                )}

                {isLoading && (
                    <div className="loading-indicator" style={{ marginTop: '20px', textAlign: 'center' }}>
                        <GlassyContainer>
                            <div style={{ padding: '20px' }}>
                                <div className="loading-spinner" style={{ 
                                    border: '4px solid #f3f3f3', 
                                    borderTop: '4px solid rgb(3, 196, 3)', 
                                    borderRadius: '50%', 
                                    width: '40px', 
                                    height: '40px', 
                                    animation: 'spin 2s linear infinite',
                                    margin: '0 auto'
                                }}></div>
                                <p style={{ marginTop: '10px' }}>Fetching UniProt summary data...</p>
                            </div>
                        </GlassyContainer>
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
                    <div className="explorer-row-5" style={{ marginTop: '20px' }}>
                        <DataViewer 
                            data={apiData} 
                            title="Complete UniProt Summary Data"
                            initiallyExpanded={false}
                        />
                    </div>
                )}
                {apiData && (
                    <div className="explorer-row-4" style={{ marginTop: '20px' }}>
                        <GlassyContainer className="exp-4-col-1">
                            <h3 style={{ marginBottom: '35px', fontWeight: 'bold' }}>Associated Files</h3>
                            <div style={{ display: 'flex', flexDirection: 'row', justifyContent: 'flex-start', flexWrap: 'wrap' }}>
                                {renderFileInfo("CIF File (.cif)", apiData?.structures[0].summary?.model_url)}
                            </div>
                        </GlassyContainer>
                    </div>
                )}
            </motion.div>
        </>
    );
};

// // PropTypes for better type checking and documentation
// UniProtViewer.propTypes = {
//     toolData: PropTypes.shape({
//         gene: PropTypes.string,
//         uniprotDescription: PropTypes.string,
//         organismScientificName: PropTypes.string,
//         uniprotId: PropTypes.string,
//         entryId: PropTypes.string,
//         sequenceVersionDate: PropTypes.string,
//         modelCreatedDate: PropTypes.string,
//         latestVersion: PropTypes.string,
//         taxId: PropTypes.string,
//         uniprotSequence: PropTypes.string,
//         pdbUrl: PropTypes.string,
//         cifUrl: PropTypes.string,
//         bcifUrl: PropTypes.string,
//         paeImageUrl: PropTypes.string,
//         paeDocUrl: PropTypes.string,
//     }),
//     initialAccessionKey: PropTypes.string,
//     hideInputBox: PropTypes.bool,
// };

export default UniProtSummaryViewer;
