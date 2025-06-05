import React, { useState, useEffect, useCallback, useRef } from "react";
import PropTypes from "prop-types";
import "./UniProtViewer.css"; // Reusing the shared styles
import SimpleInputBox from "../../../components/UI/SimpleInputBox/SimpleInputBox";
import InfoBox from "../../../components/UI/InfoBox/InfoBox";
import ThreeDViewer from "../../../components/UI/ThreeDViewer/ThreeDViewer";
import DataViewer from "../../../components/UI/DataViewer/DataViewer";
import { motion } from "framer-motion";
import { fadeInUpVariantStatic } from "../../../components/animations/framerAnim";
import GlassyContainer from "../../../components/glassy_container/gc";
import { RcsbFv } from "@rcsb/rcsb-saguaro";
import { FaDownload } from 'react-icons/fa';
import { getAlphafoldAnnotations } from "../../../services/api/mcpToolsService";

const AlphafoldAnnotationsViewer = ({ 
    toolData = null, 
    initialAccessionKey = "", 
    hideInputBox = false,
    initialAnnotationType = "MUTAGEN"
}) => {
    const [accessionKey, setAccessionKey] = useState(initialAccessionKey);
    const [inputValue, setInputValue] = useState(initialAccessionKey);
    const [annotationType, setAnnotationType] = useState(initialAnnotationType);
    const [annotationTypeInput, setAnnotationTypeInput] = useState(initialAnnotationType);
    const [apiData, setApiData] = useState(toolData);
    const [isLoading, setIsLoading] = useState(false);
    const [error, setError] = useState(null);
    const rcsbFvInstance = useRef(null);
    const annotationsViewerRef = useRef(null);

    useEffect(() => {
        console.log('AlphaFold annotations tool data received >>', toolData);
    }, [toolData]);

    useEffect(() => {
        console.log("AlphaFold annotations apiData updated >>", apiData);
    }, [apiData]);

    const handleSubmit = useCallback(async (accession, annotation_type) => {
        if (!accession || accession.trim() === "") {
            setError("Please enter a valid UniProt accession key");
            return;
        }

        if (!annotation_type || annotation_type.trim() === "") {
            setError("Please enter a valid annotation type");
            return;
        }

        try {
            setIsLoading(true);
            setError(null);
            setAccessionKey(accession);
            setAnnotationType(annotation_type);
            
            const args = { 
                qualifier: accession.trim(),
                annotation_type: annotation_type.trim()
            };
            
            console.log('Calling MCP AlphaFold annotations with args:', args);
            const response = await getAlphafoldAnnotations(args);
            const result = response?.result;
            
            if (!result) {
                console.error('No AlphaFold annotations result returned');
                setError("No AlphaFold annotations data available");
                setApiData(null);
                return;
            }

            console.log('AlphaFold annotations result:', result);
            console.log('type of AlphaFold annotations result:', typeof result);
            let actualResults = null;
            
            if (typeof result === "object") {
                // Handle different response formats
                if (result["0"]) {
                    // If result is wrapped in an object with "0" key
                    const result_data_json = JSON.parse(result["0"]);
                    if (result_data_json.error === null) {
                        actualResults = result_data_json.data;
                    } else {
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
            } else {
                setError("Failed to process AlphaFold annotations data");
                setApiData(null);
            }
        } catch (error) {
            console.error('Error fetching AlphaFold annotations:', error);
            setError(`Failed to fetch AlphaFold annotations: ${error.message}`);
            setApiData(null);
        } finally {
            setIsLoading(false);
        }
    }, []);

    // Effect to handle initial data passed as props
    useEffect(() => {
        if (toolData) {
            var tempData = null;
            if (typeof toolData === "object") {
                tempData = toolData;
            }
            setApiData(tempData);
            setError(null);
        }
    }, [toolData]);

    // Effect to handle initial accession key
    useEffect(() => {
        if (initialAccessionKey && initialAccessionKey.trim() !== "") {
            setInputValue(initialAccessionKey);
            setAccessionKey(initialAccessionKey);
            // Optionally auto-submit if we have initial data
            if (!toolData) {
                handleSubmit(initialAccessionKey, initialAnnotationType);
            }
        }
    }, [initialAccessionKey, toolData, handleSubmit, initialAnnotationType]);

    useEffect(() => {
        if (!accessionKey || accessionKey.trim() === "") {
            // Only clear data if no initial data was provided
            if (!toolData) {
                setApiData(null);
                setError(null);
            }
            if (rcsbFvInstance.current) {
                rcsbFvInstance.current.unmount();
                rcsbFvInstance.current = null;
            }
        }
    }, [accessionKey, toolData]);

    // Effect to create annotations visualization
    useEffect(() => {
        if (apiData && apiData.annotations && Array.isArray(apiData.annotations) && annotationsViewerRef.current) {
            if (rcsbFvInstance.current) {
                rcsbFvInstance.current.unmount();
            }

            const annotations = apiData.annotations;
            console.log("annotations >>", annotations);
            
            // Find the range for the board configuration
            const minPosition = Math.min(...annotations.map(ann => ann.start || 1));
            const maxPosition = Math.max(...annotations.map(ann => ann.end || ann.start || 1));
            
            const boardConfig = {
                range: {
                    min: minPosition,
                    max: maxPosition
                },
                rowTitleWidth: 170,
                includeAxis: true
            };

            // Create tracks for different annotation types
            const tracks = [];
            
            // Group annotations by type for better visualization
            const annotationsByType = annotations.reduce((acc, ann) => {
                const type = ann.type || 'Unknown';
                if (!acc[type]) {
                    acc[type] = [];
                }
                acc[type].push(ann);
                return acc;
            }, {});

            // Color palette for different annotation types
            const colors = [
                "#FF6B6B", "#4ECDC4", "#45B7D1", "#96CEB4", "#FFEAA7", 
                "#DDA0DD", "#98D8C8", "#F7DC6F", "#BB8FCE", "#85C1E9"
            ];

            Object.entries(annotationsByType).forEach(([type, typeAnnotations], index) => {
                const trackData = typeAnnotations.map(ann => ({
                    begin: ann.start || 1,
                    end: ann.end || ann.start || 1,
                    label: ann.description || ann.label || type,
                    color: colors[index % colors.length]
                }));

                tracks.push({
                    trackId: `${type.toLowerCase().replace(/\s+/g, '_')}_track`,
                    trackHeight: 40,
                    trackColor: "#FFFFFF",
                    displayType: "block",
                    displayColor: colors[index % colors.length],
                    rowTitle: type.toUpperCase(),
                    titleFlagColor: "#FFFFFF",
                    trackData: trackData
                });
            });

            if (tracks.length > 0) {
                rcsbFvInstance.current = new RcsbFv({
                    boardConfigData: boardConfig,
                    rowConfigData: tracks,
                    elementId: annotationsViewerRef.current.id
                });
            }
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

    const handleFormSubmit = () => {
        handleSubmit(inputValue, annotationTypeInput);
    };

    const renderAnnotationsSummary = () => {
        if (!apiData || !apiData.annotations || !Array.isArray(apiData.annotations)) {
            return null;
        }

        const annotations = apiData.annotations;
        const annotationsByType = annotations.reduce((acc, ann) => {
            const type = ann.type || 'Unknown';
            if (!acc[type]) {
                acc[type] = 0;
            }
            acc[type]++;
            return acc;
        }, {});

        return (
            <div className="annotations-summary" style={{ marginBottom: '20px' }}>
                <h4 style={{ marginBottom: '10px', fontWeight: '600' }}>Annotations Summary</h4>
                <div style={{ display: 'grid', gridTemplateColumns: 'repeat(auto-fit, minmax(200px, 1fr))', gap: '10px' }}>
                    {Object.entries(annotationsByType).map(([type, count]) => (
                        <div key={type} className="annotation-type-summary" style={{
                            padding: '10px',
                            backgroundColor: 'rgba(255, 255, 255, 0.1)',
                            borderRadius: '8px',
                            textAlign: 'center'
                        }}>
                            <div style={{ fontWeight: 'bold', fontSize: '1.2em' }}>{count}</div>
                            <div style={{ fontSize: '0.9em', opacity: 0.8 }}>{type}</div>
                        </div>
                    ))}
                </div>
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
                            <div style={{ marginBottom: '20px' }}>
                                <h3 style={{ marginBottom: '15px', fontWeight: '700' }}>AlphaFold Annotations</h3>
                                <SimpleInputBox
                                    value={inputValue}
                                    onChange={setInputValue}
                                    onSubmit={handleFormSubmit}
                                    header="Enter UniProt Accession Key"
                                    placeholder="Try UniProt Accession: Q5VSL9"
                                    buttonText="Fetch Annotations"
                                    isLoading={isLoading}
                                    error={error}
                                />
                                
                                {/* Annotation type input */}
                                <div style={{ marginTop: '15px' }}>
                                    <label style={{ display: 'block', marginBottom: '5px', fontSize: '0.9em', opacity: 0.8 }}>
                                        Annotation Type
                                    </label>
                                    <input
                                        type="text"
                                        value={annotationTypeInput}
                                        onChange={(e) => setAnnotationTypeInput(e.target.value)}
                                        placeholder="MUTAGEN"
                                        style={{
                                            width: '100%',
                                            padding: '8px',
                                            borderRadius: '5px',
                                            border: '1px solid rgba(255, 255, 255, 0.3)',
                                            backgroundColor: 'rgba(255, 255, 255, 0.1)',
                                            color: 'white'
                                        }}
                                    />
                                    <div style={{ fontSize: '0.8em', opacity: 0.6, marginTop: '5px' }}>
                                        Common types: MUTAGEN, DOMAIN, BINDING, ACTIVE_SITE, TRANSMEM
                                    </div>
                                </div>
                            </div>
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
                                <p style={{ marginTop: '10px' }}>Fetching AlphaFold annotations data...</p>
                            </div>
                        </GlassyContainer>
                    </div>
                )}

                {apiData && apiData.annotations && (
                    <div className="explorer-row-2" style={{ marginTop: '20px' }}>
                        <GlassyContainer>
                            <h3 style={{ marginBottom: '15px', fontWeight: '700' }}>Annotations Summary</h3>
                            {renderAnnotationsSummary()}
                        </GlassyContainer>
                    </div>
                )}

                {apiData && apiData.annotations && Array.isArray(apiData.annotations) && (
                    <div className="explorer-row-3" style={{ marginTop: '20px' }}>
                        <GlassyContainer>
                            <h3 style={{ marginBottom: '15px', fontWeight: '700' }}>Annotations Viewer</h3>
                            <p style={{ marginBottom: '15px' }}>Scroll on the track to zoom in and view annotations!</p>
                            <div
                                key={accessionKey} // Force re-mount of this div when accessionKey changes
                                id="rcsb-saguaro-annotations-viewer"
                                ref={annotationsViewerRef}
                                style={{ minHeight: '300px', backgroundColor: 'gray', borderRadius: '10px', padding: '10px', zIndex: 100 }}
                            ></div>
                        </GlassyContainer>
                    </div>
                )}

                {apiData && (
                    <div className="explorer-row-5" style={{ marginTop: '20px' }}>
                        <DataViewer 
                            data={apiData} 
                            title="Complete AlphaFold Annotations Data"
                            initiallyExpanded={false}
                        />
                    </div>
                )}
            </motion.div>
        </>
    );
};

AlphafoldAnnotationsViewer.propTypes = {
    toolData: PropTypes.object,
    initialAccessionKey: PropTypes.string,
    hideInputBox: PropTypes.bool,
    initialAnnotationType: PropTypes.string,
};

export default AlphafoldAnnotationsViewer;
