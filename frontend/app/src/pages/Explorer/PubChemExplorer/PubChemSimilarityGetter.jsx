import React, { useState, useEffect, useCallback, useRef } from "react";
import { createPortal } from "react-dom";
import "./PubChemSimilarityGetter.css";
import SimpleInputBox from "../../../components/UI/SimpleInputBox/SimpleInputBox";
import DataViewer from "../../../components/UI/DataViewer/DataViewer";
import InfoBox from "../../../components/UI/InfoBox/InfoBox";
import TwoDViewer from "../../../components/UI/TwoDViewer/TwoDViewer";
import { motion } from "framer-motion";
import { fadeInUpVariantStatic } from "../../../components/animations/framerAnim";
import GlassyContainer from "../../../components/glassy_container/gc";
import { fastSimilarity2dSearchByCid, getCompoundByCid } from "../../../services/api/mcpToolsService";
import { FaCheckDouble } from "react-icons/fa6";
import { IoWarningOutline } from "react-icons/io5";
import { call_endpoint_async } from "../../../endpoints/caller";
import { endpoints } from "../../../endpoints/endpoints";

// Enhanced CID Validator
const validateCid = (cid) => {
    if (!cid || typeof cid !== 'string') {
        return false;
    }
    
    const trimmedCid = cid.trim();
    
    // Check if empty
    if (trimmedCid === '') {
        return false;
    }
    
    // CID must be a positive integer
    const cidNumber = parseInt(trimmedCid, 10);
    return !isNaN(cidNumber) && cidNumber > 0 && cidNumber.toString() === trimmedCid;
};

// In-house PubChem Similarity Search Type Selector Component with Portal
const PubChemSimilaritySearchTypeSelector = ({ value, onChange, options, disabled = false }) => {
    const [isOpen, setIsOpen] = useState(false);
    const [dropdownPosition, setDropdownPosition] = useState({ top: 0, left: 0, width: 0 });
    const dropdownRef = useRef(null);
    const triggerRef = useRef(null);

    // Update dropdown position when opened
    const updateDropdownPosition = useCallback(() => {
        if (triggerRef.current && isOpen) {
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
        const handleEscapeKey = (event) => {
            if (event.key === 'Escape' && isOpen) {
                setIsOpen(false);
            }
        };

        if (isOpen) {
            document.addEventListener('keydown', handleEscapeKey);
            return () => document.removeEventListener('keydown', handleEscapeKey);
        }
    }, [isOpen]);

    const handleOptionClick = (optionValue) => {
        onChange(optionValue);
        setIsOpen(false);
    };

    const selectedOption = options.find(option => option.value === value);

    // Portal dropdown component
    const PortalDropdown = () => {
        if (!isOpen) return null;

        return createPortal(
            <div
                ref={dropdownRef}
                className="pubchem-similarity-selector-options-portal"
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
                    maxHeight: '200px',
                    overflowY: 'auto'
                }}
            >
                {options.map((option) => (
                    <div
                        key={option.value}
                        className={`pubchem-similarity-selector-option ${option.value === value ? 'selected' : ''}`}
                        onClick={() => handleOptionClick(option.value)}
                        role="option"
                        aria-selected={option.value === value}
                    >
                        <span className="option-label">{option.label}</span>
                        {option.value === value && (
                            <span className="check-icon">
                                <svg width="14" height="10" viewBox="0 0 14 10" fill="currentColor">
                                    <path d="M13 1L5 9L1 5" />
                                </svg>
                            </span>
                        )}
                    </div>
                ))}
            </div>,
            document.body
        );
    };

    return (
        <>
            <div
                className={`pubchem-similarity-search-selector ${disabled ? 'disabled' : ''}`}
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
                            !disabled && setIsOpen(!isOpen);
                        }
                    }}
                >
                    <span className="selected-text">
                        {selectedOption ? selectedOption.label : 'Select search type'}
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
            <PortalDropdown />
        </>
    );
};

// Search type options for the dropdown
const similaritySearchTypeOptions = [
    { value: "cid", label: "PubChem CID" },
    { value: "smiles", label: "SMILES String" }
];

const PubChemSimilarityGetter = ({
    toolData = null,
    initialSearchValue = "",
    initialSearchType = "cid",
    initialSimilarityThreshold = 90,
    hideInputBox = false
}) => {
    const [searchValue, setSearchValue] = useState(initialSearchValue);
    const [inputValue, setInputValue] = useState(initialSearchValue);
    const [searchType, setSearchType] = useState(initialSearchType);
    const [similarityThreshold, setSimilarityThreshold] = useState(initialSimilarityThreshold);
    const [apiData, setApiData] = useState(toolData);
    const [isLoading, setIsLoading] = useState(false);
    const [error, setError] = useState(null);
    const [validationError, setValidationError] = useState(null);
    
    // Enhanced validation states
    const [isValidInput, setIsValidInput] = useState(false);
    const [validationStatus, setValidationStatus] = useState('empty'); // 'empty', 'valid', 'invalid'
    const [isValidatingSmiles, setIsValidatingSmiles] = useState(false); // For backend validation loading state
    
    // States for multiple similar molecules
    const [similarMolecules, setSimilarMolecules] = useState([]);
    const [selectedMolecule, setSelectedMolecule] = useState(null);
    const [selectedIndex, setSelectedIndex] = useState(0);

    // Ref to track if search type is being set from toolData
    const isSettingFromToolData = useRef(false);

    useEffect(() => {
        console.log('PubChem similarity tool data received >>', toolData);
        if (toolData && toolData.type && toolData.content) {
            console.log('Processing PubChem similarity tool data:', toolData);
            
            // Set flag to indicate we're setting search type from toolData
            isSettingFromToolData.current = true;
            
            // Set search type based on tool data type
            setSearchType(toolData.type);

            console.log('tool data >>', toolData);
            
            if (toolData.type === "cid" || toolData.type === "smiles") {
                // Handle similarity search results
                if (Array.isArray(toolData.content) && toolData.content.length > 0) {
                    console.log('Setting similar molecules from tool data:', toolData.content);
                    setSimilarMolecules(toolData.content);
                    setSelectedMolecule(toolData.content[0]);
                    setSelectedIndex(0);
                    setApiData(null); // Clear single molecule data for multiple results view
                    // Clear any existing errors for successful data processing
                    setError(null);
                    setValidationError(null);
                } else if (toolData.content) {
                    // If it's a single molecule, treat it as a single-item array for consistent UI
                    setSimilarMolecules([toolData.content]);
                    setSelectedMolecule(toolData.content);
                    setSelectedIndex(0);
                    setApiData(null); // Clear single molecule data
                    // Clear any existing errors for successful data processing
                    setError(null);
                    setValidationError(null);
                } else {
                    // No molecules found
                    setError(`No similar molecules found for ${toolData.type} search.`);
                    setSimilarMolecules([]);
                    setSelectedMolecule(null);
                    setApiData(null);
                }
            }
            
            // Reset flag after state updates
            setTimeout(() => {
                isSettingFromToolData.current = false;
            }, 0);
        } else if (toolData) {
            // Fallback for old format toolData without type/content structure
            setApiData(toolData);
        }
    }, [toolData]);

    useEffect(() => {
        console.log('PubChem similarity api data updated >>', apiData);
        console.log('similar molecules >>', similarMolecules);
        console.log('active similar molecules >>', selectedMolecule);
    }, [apiData]);

    // Backend SMILES validation function
    const validateSmilesWithBackend = useCallback(async (smiles) => {
        if (!smiles || smiles.trim() === '') {
            return false;
        }

        try {
            setIsValidatingSmiles(true);
            const payload = {
                type: "MOL",
                value: smiles.trim()
            };
            
            const response = await call_endpoint_async(endpoints.validate, payload);
            
            if (response.data.status === "success") {
                return response.data.valid;
            }
            return false;
        } catch (error) {
            console.log('SMILES validation error:', error);
            return false;
        } finally {
            setIsValidatingSmiles(false);
        }
    }, []);

    // Enhanced input change handler with comprehensive validation
    const handleInputChange = useCallback(async (value) => {
        console.log('handleInputChange called with:', value, 'for searchType:', searchType);
        setInputValue(value);
        setError(null);
        setValidationError(null);

        // Handle empty input
        if (value.trim() === "") {
            setValidationError(null);
            setIsValidInput(false);
            setValidationStatus('empty');
            return;
        }

        const trimmedValue = value.trim();
        let isValid = false;
        let errorMessage = "";

        // Validate based on search type
        if (searchType === "cid") {
            console.log('Validating CID:', trimmedValue);
            isValid = validateCid(trimmedValue);
            console.log('CID validation result:', isValid);
            if (!isValid) {
                errorMessage = "Invalid CID format. Must be a positive integer (e.g., 2244, 5281804).";
            }
            
            // Update validation states immediately for CID
            setIsValidInput(isValid);
            setValidationStatus(isValid ? 'valid' : 'invalid');
            setValidationError(isValid ? null : errorMessage);
        } else if (searchType === "smiles") {
            // For SMILES, use backend validation
            try {
                setValidationStatus('validating'); // Show loading state
                console.log('Validating SMILES:', trimmedValue);
                isValid = await validateSmilesWithBackend(trimmedValue);
                console.log('SMILES validation result:', isValid);
                
                if (!isValid) {
                    errorMessage = "Invalid SMILES format. Please enter a valid SMILES string (e.g., CCO, C1=CC=CC=C1).";
                }
                
                // Update validation states after backend validation
                setIsValidInput(isValid);
                setValidationStatus(isValid ? 'valid' : 'invalid');
                setValidationError(isValid ? null : errorMessage);
            } catch (error) {
                console.error('SMILES validation error:', error);
                setIsValidInput(false);
                setValidationStatus('invalid');
                setValidationError("Error validating SMILES. Please check your input.");
            }
        }
    }, [searchType, validateSmilesWithBackend]);

    // Re-validate when search type changes
    useEffect(() => {
        if (inputValue.trim() !== "") {
            // Force re-validation with the current input for the new search type
            handleInputChange(inputValue);
        } else {
            // If input is empty, ensure validation states are reset properly
            setIsValidInput(false);
            setValidationStatus('empty');
            setValidationError(null);
        }
    }, [searchType, handleInputChange]);

    const handleSubmit = useCallback(async (value) => {
        console.log('handleSubmit called with value:', value);
        const trimmedValue = value.trim();
        console.log("valid state >>", isValidInput);
        
        if (!trimmedValue) {
            setError("Please enter a search value.");
            return;
        }

        // Perform immediate validation as a fallback in case the state hasn't updated yet
        let isCurrentlyValid = isValidInput;
        
        if (searchType === "cid") {
            isCurrentlyValid = validateCid(trimmedValue);
            console.log('Immediate CID validation in submit:', isCurrentlyValid);
        } else if (searchType === "smiles") {
            // For SMILES, trust the async validation state since we can't do sync validation
            // But ensure we have a valid state
            isCurrentlyValid = isValidInput && validationStatus === 'valid';
            console.log('SMILES validation state in submit:', isCurrentlyValid, 'validationStatus:', validationStatus);
        }

        if (!isCurrentlyValid) {
            if (searchType === "cid") {
                setError("Invalid CID format. CID must be a positive integer (e.g., 2244, 5281804).");
            } else if (searchType === "smiles") {
                setError("Invalid SMILES format. Please enter a valid SMILES string (e.g., CCO, C1=CC=CC=C1).");
            }
            return;
        }

        setIsLoading(true);
        setError(null);
        setValidationError(null);

        try {
            let result;

            if (searchType === "cid") {
                // For CID-based similarity search, use the PubChem fast similarity service
                result = await fastSimilarity2dSearchByCid({
                    cid: parseInt(trimmedValue, 10),
                    threshold: similarityThreshold / 100 // Convert percentage to decimal
                });
            } else if (searchType === "smiles") {
                // For SMILES, we would need a different API call
                // Since PubChem API doesn't have direct SMILES similarity search in our service,
                // we can show a message or implement it differently
                setError("SMILES-based similarity search is not yet implemented for PubChem. Please use CID-based search.");
                return;
            }

            if (result) {
                if (result.status === "success") {
                    var processible = result.result["0"];
                    processible = JSON.parse(processible);
                    
                    console.log('PubChem similarity search processible', processible);

                    if (processible.data && Array.isArray(processible.data) && processible.data.length > 0) {
                        // For each similar CID, get the full compound data
                        const similarCompounds = [];
                        
                        for (const cidData of processible.data.slice(0, 20)) { // Limit to first 20 for performance
                            try {
                                const compoundResult = await getCompoundByCid({ cid: cidData.cid });
                                if (compoundResult && compoundResult.status === "success") {
                                    const compoundData = JSON.parse(compoundResult.result["0"]);
                                    if (compoundData.data) {
                                        // Add similarity score to the compound data
                                        compoundData.data.similarity_score = cidData.similarity;
                                        compoundData.data.similar_cid = cidData.cid;
                                        similarCompounds.push(compoundData.data);
                                    }
                                }
                            } catch (error) {
                                console.warn(`Failed to get compound data for CID ${cidData.cid}:`, error);
                            }
                        }

                        if (similarCompounds.length > 0) {
                            setSimilarMolecules(similarCompounds);
                            setSelectedMolecule(similarCompounds[0]);
                            setSelectedIndex(0);
                            setApiData(null); // Clear single molecule data
                            setSearchValue(trimmedValue);
                        } else {
                            setError(`No compound data found for similar molecules of "${trimmedValue}"`);
                            setSimilarMolecules([]);
                            setSelectedMolecule(null);
                            setApiData(null);
                        }
                    } else if (processible.data && Array.isArray(processible.data) && processible.data.length === 0) {
                        // No similar molecules found
                        setError(`No similar molecules found for CID "${trimmedValue}" with similarity threshold ${similarityThreshold}%`);
                        setSimilarMolecules([]);
                        setSelectedMolecule(null);
                        setApiData(null);
                    } else {
                        // Unexpected data format
                        setError(`No similar molecules found for "${trimmedValue}"`);
                        setSimilarMolecules([]);
                        setSelectedMolecule(null);
                        setApiData(null);
                    }
                } else if (result.status === "failed") {
                    setError(result.message || "Failed to search for similar molecules");
                    setSimilarMolecules([]);
                    setSelectedMolecule(null);
                    setApiData(null);
                }
            }
        } catch (error) {
            console.error('PubChem similarity search error:', error);
            setError(error.message || "An error occurred while searching for similar molecules");
            setSimilarMolecules([]);
            setSelectedMolecule(null);
            setApiData(null);
        } finally {
            setIsLoading(false);
        }
    }, [searchType, similarityThreshold, isValidInput, validationStatus]);

    // Effect to handle initial search value
    useEffect(() => {
        if (initialSearchValue && !toolData) {
            setInputValue(initialSearchValue);
            handleSubmit(initialSearchValue);
        }
    }, [initialSearchValue, toolData, handleSubmit]);

    const getPlaceholderText = () => {
        switch (searchType) {
            case "cid":
                return "Enter PubChem CID (e.g., 2244, 5281804)";
            case "smiles":
                return "Enter SMILES string (e.g., CCO, C1=CC=CC=C1)";
            default:
                return "Enter search value";
        }
    };

    const getHeaderText = () => {
        switch (searchType) {
            case "cid":
                return "CID-based Similarity Search";
            case "smiles":
                return "SMILES-based Similarity Search";
            default:
                return "Similarity Search";
        }
    };

    // Keyboard navigation for similar molecules
    useEffect(() => {
        const handleKeyDown = (event) => {
            if (!similarMolecules || similarMolecules.length === 0) return;
            
            if (event.key === 'ArrowDown') {
                event.preventDefault();
                const nextIndex = (selectedIndex + 1) % similarMolecules.length;
                setSelectedIndex(nextIndex);
                setSelectedMolecule(similarMolecules[nextIndex]);
            } else if (event.key === 'ArrowUp') {
                event.preventDefault();
                const prevIndex = selectedIndex === 0 ? similarMolecules.length - 1 : selectedIndex - 1;
                setSelectedIndex(prevIndex);
                setSelectedMolecule(similarMolecules[prevIndex]);
            }
        };

        window.addEventListener('keydown', handleKeyDown);
        return () => window.removeEventListener('keydown', handleKeyDown);
    }, [selectedIndex, similarMolecules]);

    // Clear error message after some time
    useEffect(() => {
        if (error) {
            const timer = setTimeout(() => {
                setError(null);
            }, 10000); // Clear error after 10 seconds
            return () => clearTimeout(timer);
        }
    }, [error]);

    useEffect(() => {
        console.log("search type: ", searchType);
    }, [searchType])

    useEffect(() => {
        console.log('Similar molecules updated:', similarMolecules);
        console.log('Selected molecule changed:', selectedMolecule);
        console.log('Selected index:', selectedIndex);
    }, [similarMolecules, selectedMolecule, selectedIndex]);

    return (
        <>
            <motion.div 
                initial="hidden"
                animate="visible"
                variants={fadeInUpVariantStatic}
                className="pubchem-similarity-getter-container"
            >
                {!hideInputBox && (
                    <div className="pubchem-similarity-getter-row-1">
                        <GlassyContainer>
                            <h3 style={{ marginBottom: '16px', fontWeight: '700' }}>
                                {getHeaderText()}
                            </h3>
                            
                            <div className="pubchem-similarity-getter-search-controls">
                                {/* Row 1: Search Type and Similarity Threshold */}
                                <div className="pubchem-similarity-getter-controls-row-1">
                                    {/* Search Type Dropdown */}
                                    <div className="pubchem-similarity-getter-dropdown-section">
                                        <h4 className="pubchem-similarity-getter-dropdown-header">Search Type:</h4>
                                        <div className="pubchem-similarity-getter-dropdown">
                                            <PubChemSimilaritySearchTypeSelector
                                                value={searchType}
                                                onChange={setSearchType}
                                                options={similaritySearchTypeOptions}
                                                disabled={isLoading}
                                            />
                                        </div>
                                    </div>

                                    {/* Similarity Threshold Input */}
                                    <div className="pubchem-similarity-getter-threshold-section">
                                        <h4 className="pubchem-similarity-getter-dropdown-header">Similarity:</h4>
                                        <input
                                            type="number"
                                            min="0"
                                            max="100"
                                            value={similarityThreshold}
                                            onChange={(e) => setSimilarityThreshold(parseInt(e.target.value) || 90)}
                                            className="pubchem-similarity-getter-threshold-input"
                                            disabled={isLoading}
                                            placeholder="90"
                                        />
                                        <span className="pubchem-similarity-getter-threshold-help">
                                            %
                                        </span>
                                    </div>
                                </div>

                                {/* Row 2: Search Input with Validation */}
                                <div className="pubchem-similarity-getter-controls-row-2">
                                    <div className="pubchem-similarity-getter-input-section">
                                        <div className="pubchem-similarity-input-row">
                                            <h4 className="pubchem-similarity-getter-dropdown-header">Search Value:</h4>
                                            <div className="pubchem-similarity-input-with-validation">
                                                <input
                                                    type="text"
                                                    value={inputValue}
                                                    onChange={(e) => handleInputChange(e.target.value)}
                                                    onKeyDown={(e) => {
                                                        if (e.key === 'Enter' && isValidInput && !isLoading) {
                                                            handleSubmit(inputValue);
                                                        }
                                                    }}
                                                    placeholder={getPlaceholderText()}
                                                    disabled={isLoading}
                                                    className={`pubchem-similarity-search-input ${validationStatus}`}
                                                />
                                                <div className="pubchem-similarity-validation-icons">
                                                    {validationStatus === 'valid' && (
                                                        <FaCheckDouble className="validation-icon success" />
                                                    )}
                                                    {validationStatus === 'invalid' && (
                                                        <IoWarningOutline className="validation-icon fail" />
                                                    )}
                                                    {validationStatus === 'validating' && (
                                                        <div className="validation-icon validating">⏳</div>
                                                    )}
                                                    {validationStatus === 'empty' && inputValue.length > 0 && (
                                                        <IoWarningOutline className="validation-icon empty" />
                                                    )}
                                                </div>
                                            </div>
                                            
                                            {/* Submit button */}
                                            <button
                                                onClick={() => handleSubmit(inputValue)}
                                                disabled={!isValidInput || isLoading}
                                                className={`pubchem-similarity-submit-btn ${!isValidInput || isLoading ? 'disabled' : 'enabled'}`}
                                            >
                                                {isLoading ? 'Searching...' : 'Search'}
                                            </button>
                                        </div>
                                        
                                        {/* Validation Error Message - now properly below the input row */}
                                        {validationError && (
                                            <div className="pubchem-similarity-validation-error">
                                                <p>{validationError}</p>
                                            </div>
                                        )}
                                    </div>
                                </div>
                            </div>

                            {/* Error Display */}
                            {error && (
                                <div className="pubchem-similarity-getter-error">
                                    <p>{error}</p>
                                </div>
                            )}

                            {/* Loading State */}
                            {isLoading && (
                                <div className="pubchem-similarity-getter-loading">
                                    <p>Searching for similar molecules...</p>
                                </div>
                            )}
                        </GlassyContainer>
                    </div>
                )}

                {/* Similar Molecules Results */}
                {((similarMolecules && similarMolecules.length > 0) || error) && (
                    <div className="pubchem-similarity-getter-row-2" style={{ marginTop: '16px' }}>
                        <GlassyContainer>
                            <h3 style={{ marginBottom: '12px', fontWeight: '700' }}>
                                {similarMolecules && similarMolecules.length > 0 
                                    ? `Similar Molecules Results (${similarMolecules.length} found)` 
                                    : 'Search Results'
                                }
                            </h3>
                            
                            {/* Show error state when there's an error and no molecules */}
                            {error && (!similarMolecules || similarMolecules.length === 0) ? (
                                <div className="pubchem-similarity-getter-error" style={{ margin: '20px 0' }}>
                                    <p>{error}</p>
                                </div>
                            ) : (
                                <div className="pubchem-similarity-results-container">
                                {/* Left Panel - Molecule List */}
                                <div className="pubchem-similarity-candidates-panel">
                                    <h4 style={{
                                        color: 'var(--color-text-primary)',
                                        fontSize: '1rem',
                                        marginBottom: '10px',
                                        fontWeight: '500'
                                    }}>
                                        Similar Molecules
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
                                        {similarMolecules.map((molecule, index) => (
                                            <div 
                                                key={index}
                                                onClick={() => {
                                                    setSelectedMolecule(molecule);
                                                    setSelectedIndex(index);
                                                }}
                                                className={`pubchem-similarity-candidate-item ${selectedMolecule === molecule ? 'selected' : ''}`}
                                            >
                                                <span className="pubchem-similarity-candidate-number">
                                                    {index + 1}.
                                                </span>
                                                <div className="pubchem-similarity-candidate-info">
                                                    <span className="pubchem-similarity-candidate-id">
                                                        CID: {molecule.similar_cid || molecule.cid || 'Unknown'}
                                                    </span>
                                                    {molecule.similarity_score && (
                                                        <span className="pubchem-similarity-candidate-score">
                                                            {(molecule.similarity_score * 100).toFixed(1)}%
                                                        </span>
                                                    )}
                                                </div>
                                            </div>
                                        ))}
                                    </div>
                                </div>

                                {/* Right Panel - Info and 2D Viewer */}
                                <div className="pubchem-similarity-details-panel">
                                    {selectedMolecule ? (
                                        <>
                                            {/* Selected molecule header */}
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
                                                    Selected ({selectedIndex + 1}/{similarMolecules.length}): <span style={{ 
                                                        fontFamily: 'monospace', 
                                                        color: 'var(--color-accent)',
                                                        fontWeight: '600' 
                                                    }}>
                                                        CID {selectedMolecule.similar_cid || selectedMolecule.cid || 'Unknown'}
                                                    </span>
                                                    {selectedMolecule.similarity_score && (
                                                        <span style={{
                                                            color: 'var(--color-success)',
                                                            fontWeight: '600',
                                                            marginLeft: '8px'
                                                        }}>
                                                            ({(selectedMolecule.similarity_score * 100).toFixed(1)}% similar)
                                                        </span>
                                                    )}
                                                </h4>
                                            </div>
                                            
                                            {/* Show visualization only if we have SMILES data */}
                                            {selectedMolecule.isomeric_smiles || selectedMolecule.canonical_smiles ? (
                                                <div className="pubchem-similarity-details-row">
                                                    {/* InfoBox */}
                                                    <div style={{ flex: '4' }}>
                                                        <InfoBox 
                                                            activeMol={selectedMolecule.isomeric_smiles || selectedMolecule.canonical_smiles} 
                                                            isValidMol={true} 
                                                            infoType={"MOL"} 
                                                        />
                                                    </div>
                                                    
                                                    {/* TwoDViewer */}
                                                    <div style={{ flex: '6' }}>
                                                        <TwoDViewer 
                                                            activeMol={selectedMolecule.isomeric_smiles || selectedMolecule.canonical_smiles} 
                                                            isValidMol={true} 
                                                            visType={"MOL"} 
                                                        />
                                                    </div>
                                                </div>
                                            ) : (
                                                <div className="pubchem-similarity-placeholder">
                                                    No structure data available for this molecule
                                                </div>
                                            )}
                                        </>
                                    ) : (
                                        <div className="pubchem-similarity-placeholder">
                                            Select a molecule from the left panel to view details
                                        </div>
                                    )}
                                </div>
                            </div>
                            )}
                        </GlassyContainer>
                    </div>
                )}

                {/* DataViewer - show selected molecule data */}
                {selectedMolecule && (
                    <div className="pubchem-similarity-getter-row-3" style={{ marginTop: '20px' }}>
                        <DataViewer
                            data={selectedMolecule}
                            title={`Similar Molecule Data${searchValue ? ` for "${searchValue}"` : ''}${
                                selectedMolecule ? 
                                ` - CID ${selectedMolecule.similar_cid || selectedMolecule.cid || 'Unknown'}` : 
                                ''
                            }`}
                            initiallyExpanded={true}
                        />
                    </div>
                )}
            </motion.div>
        </>
    );
};

export default PubChemSimilarityGetter;
