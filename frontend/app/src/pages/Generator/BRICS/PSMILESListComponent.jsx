import React, { useEffect, useState, useRef, useCallback } from "react";
import { createPortal } from "react-dom";
import { motion } from "framer-motion";
import MolInputBox from "../../../components/UI/InputBox/InputBox";
import InfoBox from "../../../components/UI/InfoBox/InfoBox";
import TwoDViewer from "../../../components/UI/TwoDViewer/TwoDViewer";
import GlassyContainer from "../../../components/glassy_container/gc";
import { getBricsCandidates } from "../../../services/api/mcpToolsService";
import "./PSMILESListComponent.css";

// Standalone Portal Dropdown Menu for BRICS Input Type (adapted from PubChemGetter.jsx)
const StandaloneBRICSPortalDropdownMenu = ({
    isOpen,
    options,
    currentValue,
    onOptionClick,
    dropdownPosition,
    dropdownRef
}) => {
    if (!isOpen || !options || options.length === 0) return null;

    return createPortal(
        <div
            ref={dropdownRef}
            className="brics-search-dropdown-portal"
            style={{
                position: 'fixed',
                top: `${dropdownPosition.top}px`,
                left: `${dropdownPosition.left}px`,
                width: `${dropdownPosition.width}px`,
                zIndex: 10000,
                backgroundColor: 'var(--color-bg-primary)',
                border: '1px solid var(--color-accent)',
                borderTop: 'none',
                borderRadius: '0 0 8px 8px',
                boxShadow: '0 4px 16px var(--shadow-color)',
                maxHeight: '300px',
                overflowY: 'auto'
            }}
        >
            {options.map((option) => (
                <div
                    key={option.value}
                    className={`brics-selector-option ${option.value === currentValue ? 'selected' : ''}`}
                    onClick={() => onOptionClick(option.value)}
                    role="option"
                    aria-selected={option.value === currentValue}
                    style={{
                        display: 'flex',
                        alignItems: 'flex-start',
                        justifyContent: 'space-between',
                        padding: '12px 14px',
                        cursor: 'pointer',
                        transition: 'all 0.2s ease',
                        borderBottom: '1px solid var(--c-light-border)',
                        backgroundColor: option.value === currentValue ? 'var(--color-accent)' : 'transparent',
                        color: option.value === currentValue ? 'white' : 'var(--color-text-primary)',
                        minHeight: '48px',
                        flexShrink: 0
                    }}
                    onMouseEnter={(e) => {
                        if (option.value !== currentValue) {
                            e.target.style.backgroundColor = 'var(--glassy-color)';
                        }
                    }}
                    onMouseLeave={(e) => {
                        if (option.value !== currentValue) {
                            e.target.style.backgroundColor = 'transparent';
                        }
                    }}
                >
                    <div style={{
                        display: 'flex',
                        flexDirection: 'column',
                        gap: '2px',
                        flex: '1',
                        minWidth: 0
                    }}>
                        <span className="option-label" style={{
                            fontSize: '0.85rem',
                            fontWeight: '500',
                            color: 'inherit',
                            lineHeight: '1.2'
                        }}>
                            {option.label}
                        </span>
                        {option.description && (
                            <span className="option-description" style={{
                                fontSize: '0.75rem',
                                color: option.value === currentValue ? 'rgba(255, 255, 255, 0.8)' : 'var(--color-text-secondary)',
                                lineHeight: '1.2'
                            }}>
                                {option.description}
                            </span>
                        )}
                    </div>
                    {option.value === currentValue && (
                        <span className="check-icon" style={{
                            color: 'white',
                            fontWeight: 'bold',
                            marginLeft: '8px',
                            flexShrink: 0
                        }}>
                            <svg width="14" height="10" viewBox="0 0 14 10" fill="currentColor">
                                <path d="M13 1L5 9L1 5" stroke="currentColor" strokeWidth="2" fill="none" strokeLinecap="round" strokeLinejoin="round" />
                            </svg>
                        </span>
                    )}
                </div>
            ))}
        </div>,
        document.body
    );
};

// BRICS Input Type Selector Component with Portal
const BRICSInputTypeSelector = ({ value, onChange, options, disabled = false }) => {
    const [isOpen, setIsOpen] = useState(false);
    const [dropdownPosition, setDropdownPosition] = useState({ top: 0, left: 0, width: 0 });
    const dropdownRef = useRef(null);
    const triggerRef = useRef(null);

    // Update dropdown position when opened
    const updateDropdownPosition = useCallback(() => {
        if (isOpen && triggerRef.current) {
            const rect = triggerRef.current.getBoundingClientRect();
            setDropdownPosition({
                top: rect.bottom,
                left: rect.left,
                width: rect.width
            });
        }
    }, [isOpen]);

    // Update position when dropdown opens or on scroll/resize
    useEffect(() => {
        if (isOpen) {
            updateDropdownPosition();
            const handleScroll = () => updateDropdownPosition();
            const handleResize = () => updateDropdownPosition();
            
            window.addEventListener('scroll', handleScroll, true);
            window.addEventListener('resize', handleResize);
            
            return () => {
                window.removeEventListener('scroll', handleScroll, true);
                window.removeEventListener('resize', handleResize);
            };
        }
    }, [isOpen, updateDropdownPosition]);

    // Close dropdown when clicking outside
    useEffect(() => {
        const handleClickOutside = (event) => {
            if (dropdownRef.current && !dropdownRef.current.contains(event.target) &&
                triggerRef.current && !triggerRef.current.contains(event.target)) {
                setIsOpen(false);
            }
        };

        if (isOpen) {
            document.addEventListener('mousedown', handleClickOutside);
            return () => document.removeEventListener('mousedown', handleClickOutside);
        }
    }, [isOpen]);

    // Close dropdown on escape key
    useEffect(() => {
        const handleKeyDown = (event) => {
            if (event.key === 'Escape') {
                setIsOpen(false);
            }
        };

        if (isOpen) {
            document.addEventListener('keydown', handleKeyDown);
            return () => document.removeEventListener('keydown', handleKeyDown);
        }
    }, [isOpen]);

    const handleOptionClick = (optionValue) => {
        onChange({ target: { value: optionValue } });
        setIsOpen(false);
    };

    const selectedOption = options.find(option => option.value === value);

    return (
        <>
            <div
                className={`brics-search-selector ${disabled ? 'disabled' : ''}`}
                ref={triggerRef}
            >
                <div
                    className={`selector-trigger ${isOpen ? 'open' : ''}`}
                    onClick={() => !disabled && setIsOpen(!isOpen)}
                    role="button"
                    tabIndex={disabled ? -1 : 0}
                    onKeyDown={(e) => {
                        if (e.key === 'Enter' || e.key === ' ') {
                            e.preventDefault();
                            if (!disabled) setIsOpen(!isOpen);
                        }
                    }}
                >
                    <span className="selected-text">
                        {selectedOption ? selectedOption.label : 'Select input type'}
                    </span>
                    <span className={`dropdown-arrow ${isOpen ? 'open' : ''}`}>
                        <svg
                            width="12"
                            height="8"
                            viewBox="0 0 12 8"
                            fill="currentColor"
                            style={{
                                transition: 'transform 0.2s ease',
                                transform: isOpen ? 'rotate(180deg)' : 'rotate(0deg)'
                            }}
                        >
                            <path d="M6 8L0 0h12z" />
                        </svg>
                    </span>
                </div>
            </div>
            <StandaloneBRICSPortalDropdownMenu
                isOpen={isOpen}
                options={options}
                currentValue={value}
                onOptionClick={handleOptionClick}
                dropdownPosition={dropdownPosition}
                dropdownRef={dropdownRef}
            />
        </>
    );
};

// Input type options for the dropdown
const inputTypeOptions = [
    { 
        value: "molecule", 
        label: "Molecule", 
        description: "Small molecules/SMILES" 
    },
    { 
        value: "polymer", 
        label: "Polymer", 
        description: "Polymers/PSMILES" 
    }
];

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

    // Animation variants
    const containerVariants = {
        hidden: { opacity: 0, y: 30 },
        visible: {
            opacity: 1,
            y: 0,
            transition: {
                duration: 0.6,
                ease: "easeOut",
                staggerChildren: 0.1
            }
        }
    };

    const sectionVariants = {
        hidden: { opacity: 0, y: 20 },
        visible: {
            opacity: 1,
            y: 0,
            transition: {
                duration: 0.5,
                ease: "easeOut"
            }
        }
    };

    const candidateItemVariants = {
        hidden: { opacity: 0, scale: 0.8, y: 10 },
        visible: {
            opacity: 1,
            scale: 1,
            y: 0,
            transition: {
                duration: 0.3,
                ease: "easeOut"
            }
        }
    };

    const listContainerVariants = {
        hidden: { opacity: 0 },
        visible: {
            opacity: 1,
            transition: {
                staggerChildren: 0.05
            }
        }
    };

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
        <motion.div 
            variants={containerVariants}
            initial="hidden"
            animate="visible"
            style={{ 
                padding: '20px', 
                backgroundColor: 'var(--color-bg-primary)',
                color: 'var(--color-text-primary)',
                width: '100%',
                border: '1px solid var(--c-light-border)',
                borderRadius: '12px'
            }}
        >
            <motion.h2 
                variants={sectionVariants}
                style={{ 
                    marginBottom: '16px', 
                    color: 'var(--color-text-primary)',
                    fontSize: '1.5rem',
                    fontWeight: '600'
                }}
            >
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
            </motion.h2>

            {/* Input Section */}
            <motion.div 
                variants={sectionVariants}
                style={{ 
                    marginBottom: '16px', 
                    display: 'flex', 
                    gap: '12px', 
                    alignItems: 'flex-end',
                    flexWrap: 'wrap'
                }}
            >
                {/* Dropdown for candidate type - moved to the left */}
                <div className="brics-input-container">
                    <label className="brics-type-label">
                        Input Type
                    </label>
                    <BRICSInputTypeSelector
                        value={candidateType}
                        onChange={(e) => setCandidateType(e.target.value)}
                        options={inputTypeOptions}
                        disabled={status === "loading"}
                    />
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
            </motion.div>

            {/* List Section */}
            <motion.div 
                variants={sectionVariants}
                style={{ marginBottom: '16px' }}
            >
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
            </motion.div>

            {/* Action Section */}
            <motion.div 
                variants={sectionVariants}
                style={{ 
                    display: 'flex', 
                    gap: '12px', 
                    marginBottom: '16px',
                    flexWrap: 'wrap'
                }}
            >
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
            </motion.div>

            {/* Error Message */}
            {status === "error" && error && (
                <motion.div 
                    variants={sectionVariants}
                    initial="hidden"
                    animate="visible"
                    style={{ 
                        color: 'var(--color-alert)',
                        backgroundColor: 'var(--color-bg-secondary)',
                        padding: '12px',
                        border: '1px solid var(--color-alert)',
                        borderRadius: '8px',
                        marginBottom: '16px'
                    }}
                >
                    <p style={{ margin: '0', fontWeight: '500' }}>
                        <strong>Error:</strong> {error}
                    </p>
                </motion.div>
            )}

            {/* Results Section */}
            {status === "success" && results && (
                <motion.div 
                    variants={sectionVariants}
                    initial="hidden"
                    animate="visible"
                    style={{ 
                        backgroundColor: 'var(--color-bg-secondary)',
                        borderRadius: '8px',
                        border: '1px solid var(--c-light-border)',
                        padding: '16px',
                        marginBottom: '20px'
                    }}
                >
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
                </motion.div>
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
        </motion.div>
    );
};

export default PSMILESListComponent;
