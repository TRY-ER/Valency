import React, { useEffect, useState } from "react";
import MolInputBox from "../../../components/UI/InputBox/InputBox";
import InfoBox from "../../../components/UI/InfoBox/InfoBox";
import TwoDViewer from "../../../components/UI/TwoDViewer/TwoDViewer";
import GlassyContainer from "../../../components/glassy_container/gc";
import { getBricsCandidates } from "../../../services/api/mcpToolsService";
import "./PSMILESListComponent.css";

/**
 * PSMILESListComponent - A React component for BRICS candidate generation
 * 
 * @param {Object} props - Component props
 * @param {Object} [props.toolData=null] - Pre-loaded tool data with same schema as API response.
 *                                        When provided, component will initialize with this data
 *                                        and skip input sections.
 * @returns {JSX.Element} The rendered component
 */
/**
 * PSMILESListComponent - BRICS Candidate Generation Component
 * 
 * @param {Object} props - Component properties
 * @param {Object|null} props.toolData - Optional pre-loaded tool data with the same schema as API response.
 *                                       When provided, initializes the component with results while maintaining 
 *                                       full functionality for user interaction (reset, new queries, etc.)
 * 
 * Expected toolData schema:
 * {
 *   result: {
 *     "0": JSON.stringify({
 *       candidates: ["SMILES1", "SMILES2", ...]
 *     })
 *   }
 * }
 */
const PSMILESListComponent = ({ toolData = null }) => {
    const [currentCandidate, setCurrentCandidate] = useState("");
    const [isValidCandidate, setIsValidCandidate] = useState(false);
    const [candidatesList, setCandidatesList] = useState([]);
    const [candidateType, setCandidateType] = useState("molecule"); // "molecule" or "polymer"
    const [status, setStatus] = useState("init"); // "init", "loading", "success", "error"
    const [results, setResults] = useState(null);
    const [error, setError] = useState("");
    const [showVals, setShowValues] = useState(null);
    const [selectedCandidate, setSelectedCandidate] = useState("");
    const [isSelectedCandidateValid, setIsSelectedCandidateValid] = useState(false);
    const [selectedIndex, setSelectedIndex] = useState(0);
    // const [showRawResults, setShowRawResults] = useState(false);

    // Initialize component with toolData if provided
    useEffect(() => {
        if (toolData) {
            console.log("Initializing with toolData:", toolData.candidates);
            setResults(toolData);
            setStatus("success");
            
            // Extract candidates from toolData using the same logic as API response
            if (toolData?.candidates) {
                var candidates = toolData.candidates || [];
                setShowValues(candidates);
                // Automatically select the first candidate
                if (candidates.length > 0) {
                    setSelectedCandidate(candidates[0]);
                    setIsSelectedCandidateValid(true);
                    setSelectedIndex(0);
                }
            }
        }
    }, [toolData]);

    const handleAddCandidate = () => {
        if (isValidCandidate && currentCandidate.trim() !== "") {
            setCandidatesList(prevList => [...prevList, currentCandidate.trim()]);
            setCurrentCandidate("");
            setIsValidCandidate(false); // Reset for next input
        } else {
            console.warn("Cannot add invalid or empty candidate string.");
            // Optionally, set an error message to display to the user
        }
    };

    useEffect(() => {
        // Only process results if they're not from initial toolData
        if(results?.result && !toolData){
            console.log("o vlaues >>", results.result["0"])
            console.log("o vlaues type >>", typeof results.result["0"])
            const candidate_obj = JSON.parse(results.result["0"])
            var candidates = [];
            if (candidate_obj) {
                candidates = candidate_obj.candidates || []; 
            }
            setShowValues(candidates);
            // Automatically select the first candidate
            if (candidates.length > 0) {
                setSelectedCandidate(candidates[0]);
                setIsSelectedCandidateValid(true);
                setSelectedIndex(0);
            }
        }
    }, [results, toolData])

    // Keyboard navigation for candidate selection
    useEffect(() => {
        const handleKeyDown = (event) => {
            if (!showVals || showVals.length === 0) return;
            
            if (event.key === 'ArrowDown') {
                event.preventDefault();
                const nextIndex = (selectedIndex + 1) % showVals.length;
                setSelectedIndex(nextIndex);
                setSelectedCandidate(showVals[nextIndex]);
                setIsSelectedCandidateValid(true);
            } else if (event.key === 'ArrowUp') {
                event.preventDefault();
                const prevIndex = selectedIndex === 0 ? showVals.length - 1 : selectedIndex - 1;
                setSelectedIndex(prevIndex);
                setSelectedCandidate(showVals[prevIndex]);
                setIsSelectedCandidateValid(true);
            }
        };

        window.addEventListener('keydown', handleKeyDown);
        return () => window.removeEventListener('keydown', handleKeyDown);
    }, [selectedIndex, showVals])

    const handleRemoveCandidate = (indexToRemove) => {
        setCandidatesList(prevList => prevList.filter((_, index) => index !== indexToRemove));
    };

    const handleProcessList = async () => {
        if (candidatesList.length === 0) {
            setError("Please add at least one candidate string to the list.");
            setStatus("error"); // Set status to error to display the message
            return;
        }
        setStatus("loading");
        setError("");
        setResults(null);

        try {
            const response = await getBricsCandidates({
                smiles_list: candidatesList,
                is_polymer: candidateType === "polymer" // true for polymer, false for molecule
            });
            setResults(response);
            setStatus("success");
        } catch (err) {
            console.error("Error processing candidates list with BRICS:", err);
            setError(err.message || "Failed to process the list. Check console for details.");
            setStatus("error");
        }
    };

    const handleReset = () => {
        setCurrentCandidate("");
        setIsValidCandidate(false);
        setCandidatesList([]);
        setCandidateType("molecule");
        setStatus("init");
        setResults(null);
        setError("");
        setShowValues(null);
        setSelectedCandidate("");
        setIsSelectedCandidateValid(false);
        setSelectedIndex(0);
        // setShowRawResults(false);
        // Reset toolData effects - allow user to start fresh even if toolData was provided initially
    };

    useEffect(() => {
        console.log('show vals >>', showVals)
    }, [showVals])

    return (
        <div style={{ 
            padding: '20px', 
            backgroundColor: 'var(--color-bg-primary)',
            color: 'var(--color-text-primary)',
            width: '100%',
            border: '1px solid var(--c-light-border)',
            borderRadius: '12px'
        }}>
            <h2 style={{ 
                marginBottom: '16px', 
                color: 'var(--color-text-primary)',
                fontSize: '1.5rem',
                fontWeight: '600'
            }}>
                BRICS Candidate Generation from Candidates List
                {toolData && (
                    <span style={{
                        fontSize: '0.8rem',
                        fontWeight: '400',
                        color: 'var(--color-accent)',
                        marginLeft: '12px',
                        backgroundColor: 'var(--color-bg-secondary)',
                        padding: '4px 8px',
                        borderRadius: '4px',
                        border: '1px solid var(--color-accent)'
                    }}>
                        Pre-loaded Results
                    </span>
                )}
            </h2>

            {/* Input Section */}
            <div style={{ 
                marginBottom: '16px', 
                display: 'flex', 
                gap: '12px', 
                alignItems: 'flex-end',
                flexWrap: 'wrap'
            }}>
                {/* Dropdown for candidate type - moved to the left */}
                <div className="brics-input-container">
                    <label className="brics-type-label">
                        Input Type
                    </label>
                    <div className="brics-custom-dropdown">
                        <select
                            value={candidateType}
                            onChange={(e) => setCandidateType(e.target.value)}
                            className="brics-dropdown"
                        >
                            <option value="molecule">Molecule</option>
                            <option value="polymer">Polymer</option>
                        </select>
                    </div>
                    <div className="brics-type-description">
                        {candidateType === "molecule" ? "Small molecules/SMILES" : "Polymers/PSMILES"}
                    </div>
                </div>
                
                <div style={{ flex: '1', minWidth: '300px' }}>
                    <MolInputBox
                        activeMol={currentCandidate}
                        setActiveMol={setCurrentCandidate}
                        isValidMol={isValidCandidate}
                        setIsValidMol={setIsValidCandidate}
                        inputType={"MOL"}
                    />
                </div>
                
                <button 
                    onClick={handleAddCandidate} 
                    disabled={!isValidCandidate || currentCandidate.trim() === ""}
                    style={{ 
                        padding: '10px 16px',
                        backgroundColor: (!isValidCandidate || currentCandidate.trim() === "") 
                            ? 'var(--c-deep-light)' 
                            : 'var(--color-accent)',
                        color: (!isValidCandidate || currentCandidate.trim() === "") 
                            ? 'var(--color-text-secondary)' 
                            : 'white',
                        border: 'none',
                        borderRadius: '8px',
                        cursor: (!isValidCandidate || currentCandidate.trim() === "") ? 'not-allowed' : 'pointer',
                        fontWeight: '500',
                        transition: 'all 0.2s ease',
                        minWidth: '120px'
                    }}
                >
                    + Add to List
                </button>
            </div>

            {/* List Section */}
            <div style={{ marginBottom: '16px' }}>
                <h3 style={{ 
                    color: 'var(--color-text-primary)',
                    fontSize: '1.1rem',
                    marginBottom: '12px',
                    fontWeight: '500'
                }}>
                    Candidates List ({candidatesList.length}) - {candidateType === "molecule" ? "Molecules" : "Polymers"}
                </h3>
                
                {candidatesList.length === 0 ? (
                    <div style={{ 
                        padding: '20px',
                        backgroundColor: 'var(--color-bg-secondary)',
                        borderRadius: '8px',
                        border: '1px solid var(--c-light-border)',
                        textAlign: 'center',
                        color: 'var(--color-text-secondary)'
                    }}>
                        No candidates added yet. Add some {candidateType === "molecule" ? "SMILES" : "PSMILES"} to get started!
                    </div>
                ) : (
                    <div style={{ 
                        backgroundColor: 'var(--color-bg-secondary)',
                        borderRadius: '8px',
                        border: '1px solid var(--c-light-border)',
                        padding: '12px',
                        display: 'flex',
                        flexWrap: 'wrap',
                        gap: '8px',
                        overflowY: 'auto'
                    }}>
                        {candidatesList.map((candidate, index) => (
                            <div key={index} style={{ 
                                display: 'inline-flex',
                                alignItems: 'center',
                                backgroundColor: 'var(--color-bg-primary)',
                                border: '1px solid var(--c-light-border)',
                                borderRadius: '20px',
                                padding: '6px 12px',
                                maxWidth: '300px',
                                position: 'relative'
                            }}>
                                <span style={{ 
                                    marginRight: '6px',
                                    color: 'var(--color-accent)',
                                    fontWeight: '500',
                                    fontSize: '0.8rem',
                                    minWidth: '18px'
                                }}>
                                    {index + 1}.
                                </span>
                                <span style={{ 
                                    color: 'var(--color-text-primary)',
                                    fontFamily: 'monospace',
                                    fontSize: '0.8rem',
                                    marginRight: '8px',
                                    overflow: 'hidden',
                                    textOverflow: 'ellipsis',
                                    whiteSpace: 'nowrap',
                                    maxWidth: '200px'
                                }}>
                                    {candidate}
                                </span>
                                <button 
                                    onClick={() => handleRemoveCandidate(index)}
                                    style={{ 
                                        width: '18px',
                                        height: '18px',
                                        borderRadius: '50%',
                                        backgroundColor: 'var(--color-alert)',
                                        color: 'white',
                                        border: 'none',
                                        cursor: 'pointer',
                                        fontSize: '12px',
                                        fontWeight: 'bold',
                                        display: 'flex',
                                        alignItems: 'center',
                                        justifyContent: 'center',
                                        transition: 'all 0.2s ease',
                                        flexShrink: 0
                                    }}
                                    onMouseOver={(e) => {
                                        e.target.style.backgroundColor = 'var(--color-alert)';
                                        e.target.style.opacity = '0.8';
                                        e.target.style.transform = 'scale(1.1)';
                                    }}
                                    onMouseOut={(e) => {
                                        e.target.style.opacity = '1';
                                        e.target.style.transform = 'scale(1)';
                                    }}
                                    title={`Remove: ${candidate}`}
                                >
                                    ×
                                </button>
                            </div>
                        ))}
                    </div>
                )}
            </div>

            {/* Action Section */}
            <div style={{ 
                display: 'flex', 
                gap: '12px', 
                marginBottom: '16px',
                flexWrap: 'wrap'
            }}>
                <button 
                    onClick={handleProcessList} 
                    disabled={candidatesList.length === 0 || status === "loading"}
                    style={{ 
                        padding: '12px 20px',
                        backgroundColor: (candidatesList.length === 0 || status === "loading") 
                            ? 'var(--c-deep-light)' 
                            : 'var(--color-accent)',
                        color: (candidatesList.length === 0 || status === "loading") 
                            ? 'var(--color-text-secondary)' 
                            : 'white',
                        border: 'none',
                        borderRadius: '8px',
                        cursor: (candidatesList.length === 0 || status === "loading") ? 'not-allowed' : 'pointer',
                        fontWeight: '500',
                        fontSize: '1rem',
                        transition: 'all 0.2s ease',
                        width: 'auto'
                    }}
                >
                    {status === "loading" ? "Processing..." : "Process List with BRICS"}
                </button>
                <button 
                    onClick={handleReset} 
                    disabled={status === "loading"}
                    style={{ 
                        padding: '12px 20px',
                        backgroundColor: status === "loading" 
                            ? 'var(--c-deep-light)' 
                            : 'var(--color-bg-secondary)',
                        color: 'var(--color-text-primary)',
                        border: '1px solid var(--c-light-border)',
                        borderRadius: '8px',
                        cursor: status === "loading" ? 'not-allowed' : 'pointer',
                        fontWeight: '500',
                        fontSize: '1rem',
                        transition: 'all 0.2s ease',
                        width: 'auto'
                    }}
                >
                    Reset
                </button>
            </div>

            {/* Error Message */}
            {status === "error" && error && (
                <div style={{ 
                    color: 'var(--color-alert)',
                    backgroundColor: 'var(--color-bg-secondary)',
                    padding: '12px',
                    border: '1px solid var(--color-alert)',
                    borderRadius: '8px',
                    marginBottom: '16px'
                }}>
                    <p style={{ margin: '0', fontWeight: '500' }}>
                        <strong>Error:</strong> {error}
                    </p>
                </div>
            )}

            {/* Results Section */}
            {status === "success" && results && (
                <div style={{ 
                    backgroundColor: 'var(--color-bg-secondary)',
                    borderRadius: '8px',
                    border: '1px solid var(--c-light-border)',
                    padding: '16px',
                    marginBottom: '20px'
                }}>
                    <h3 style={{ 
                        color: 'var(--color-text-primary)',
                        fontSize: '1.1rem',
                        marginBottom: '12px',
                        fontWeight: '500'
                    }}>
                        BRICS Algorithm Results {showVals && showVals.length > 0 ? `(${showVals.length} candidates generated)` : ''}:
                    </h3>
                    
                    {/* No candidates message */}
                    {(!showVals || showVals.length === 0) && (
                        <div style={{ 
                            padding: '20px',
                            backgroundColor: 'var(--color-bg-primary)',
                            borderRadius: '8px',
                            border: '1px solid var(--c-light-border)',
                            textAlign: 'center',
                            color: 'var(--color-text-secondary)'
                        }}>
                            <p style={{ margin: '0', fontSize: '1rem' }}>
                                No candidates generated. The BRICS algorithm did not find any valid candidates for the provided input.
                            </p>
                        </div>
                    )}
                    
                    {/* Existing results display - only show when we have candidates */}
                    {showVals && showVals.length > 0 && (
                        <>
                            {/* Results Display Layout - Similar to MolExplorer */}
                            <div className="brics-results-container">
                        {/* Left Panel - Candidate List */}
                        <div className="brics-candidates-panel">
                            <h4 style={{
                                color: 'var(--color-text-primary)',
                                fontSize: '1rem',
                                marginBottom: '10px',
                                fontWeight: '500'
                            }}>
                                Candidate SMILES
                            </h4>
                            <p style={{
                                color: 'var(--color-text-secondary)',
                                fontSize: '0.8rem',
                                marginBottom: '12px',
                                fontStyle: 'italic'
                            }}>
                                💡 Use ↑↓ arrow keys to navigate
                            </p>
                            <div>
                                {showVals.map((candidate, index) => (
                                    <div 
                                        key={index}
                                        onClick={() => {
                                            setSelectedCandidate(candidate);
                                            setIsSelectedCandidateValid(true);
                                            setSelectedIndex(index);
                                        }}
                                        className={`brics-candidate-item ${selectedCandidate === candidate ? 'selected' : ''}`}
                                    >
                                        <span className="brics-candidate-number">
                                            {index + 1}.
                                        </span>
                                        <span className="brics-candidate-smiles">
                                            {candidate}
                                        </span>
                                    </div>
                                ))}
                            </div>
                        </div>

                        {/* Right Panel - Info and 2D Viewer */}
                        <div className="brics-details-panel">
                            {selectedCandidate ? (
                                <>
                                    {/* Selected candidate header */}
                                    <div style={{
                                        backgroundColor: 'var(--color-bg-primary)',
                                        borderRadius: '6px',
                                        border: '1px solid var(--c-light-border)',
                                        padding: '10px 12px',
                                        marginBottom: '8px'
                                    }}>
                                        <h4 style={{
                                            color: 'var(--color-text-primary)',
                                            fontSize: '0.9rem',
                                            margin: '0',
                                            fontWeight: '500'
                                        }}>
                                            Selected ({selectedIndex + 1}/{showVals.length}): <span style={{ 
                                                fontFamily: 'monospace', 
                                                color: 'var(--color-accent)',
                                                fontWeight: '600' 
                                            }}>
                                                {selectedCandidate}
                                            </span>
                                        </h4>
                                    </div>
                                    
                                    <div className="brics-details-row">
                                        {/* InfoBox */}
                                        <div style={{ flex: '4' }}>
                                            <InfoBox 
                                                activeMol={selectedCandidate} 
                                                isValidMol={isSelectedCandidateValid} 
                                                infoType={"MOL"} 
                                            />
                                        </div>
                                        
                                        {/* TwoDViewer */}
                                        <div style={{ flex: '6' }}>
                                            <TwoDViewer 
                                                activeMol={selectedCandidate} 
                                                isValidMol={isSelectedCandidateValid} 
                                                visType={"MOL"} 
                                            />
                                        </div>
                                    </div>
                                </>
                            ) : (
                                <div className="brics-placeholder">
                                    Select a candidate SMILES from the left panel to view details
                                </div>
                            )}
                        </div>
                    </div>
                        </>
                    )}
                </div>
            )}

            {/* Raw JSON Results (Collapsible) */}
            {/* {status === "success" && results && (
                <div style={{ 
                    backgroundColor: 'var(--color-bg-secondary)',
                    borderRadius: '8px',
                    border: '1px solid var(--c-light-border)',
                    padding: '16px'
                }}>
                    <div 
                        onClick={() => setShowRawResults(!showRawResults)}
                        style={{
                            display: 'flex',
                            alignItems: 'center',
                            cursor: 'pointer',
                            marginBottom: showRawResults ? '12px' : '0'
                        }}
                    >
                        <span style={{
                            color: 'var(--color-text-primary)',
                            fontSize: '1.1rem',
                            fontWeight: '500',
                            marginRight: '8px'
                        }}>
                            Raw API Response
                        </span>
                        <span style={{
                            color: 'var(--color-accent)',
                            fontSize: '1rem',
                            transition: 'transform 0.2s ease',
                            transform: showRawResults ? 'rotate(90deg)' : 'rotate(0deg)'
                        }}>
                            ▶
                        </span>
                    </div>
                    
                    {showRawResults && (
                        <div style={{ 
                            backgroundColor: 'var(--color-bg-primary)',
                            padding: '12px',
                            border: '1px solid var(--c-light-border)',
                            borderRadius: '6px',
                            overflowX: 'auto',
                            maxHeight: '400px',
                            overflowY: 'auto'
                        }}>
                            <pre style={{ 
                                margin: '0',
                                fontFamily: 'monospace',
                                fontSize: '0.85rem',
                                color: 'var(--color-text-primary)',
                                whiteSpace: 'pre-wrap',
                                wordBreak: 'break-word'
                            }}>
                                {JSON.stringify(results, null, 2)}
                            </pre>
                        </div>
                    )}
                </div>
            )} */}
        </div>
    );
};

export default PSMILESListComponent;
