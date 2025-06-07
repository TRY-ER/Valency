import React, { useState, useEffect, useCallback, useRef } from "react";
import { createPortal } from "react-dom";
// import PropTypes from "prop-types";
import "./ChemBLSimilarityGetter.css";
import SimpleInputBox from "../../../components/UI/SimpleInputBox/SimpleInputBox";
import DataViewer from "../../../components/UI/DataViewer/DataViewer";
import InfoBox from "../../../components/UI/InfoBox/InfoBox";
import TwoDViewer from "../../../components/UI/TwoDViewer/TwoDViewer";
import { motion } from "framer-motion";
import { fadeInUpVariantStatic } from "../../../components/animations/framerAnim";
import GlassyContainer from "../../../components/glassy_container/gc";
import { findSimilarMoleculesBySmiles, findSimilarMoleculesByChemblId } from "../../../services/api/mcpToolsService";
import { FaCheckDouble } from "react-icons/fa6";
import { IoWarningOutline } from "react-icons/io5";
import { call_endpoint_async } from "../../../endpoints/caller";
import { endpoints } from "../../../endpoints/endpoints";

// Enhanced ChEMBL ID Validator
const validateChemblId = (chemblId) => {
    if (!chemblId || typeof chemblId !== 'string') {
        return false;
    }
    
    const trimmedId = chemblId.trim();
    
    // Check if empty
    if (trimmedId === '') {
        return false;
    }
    
    // ChEMBL ID must start with "CHEMBL" followed by one or more digits
    // Case insensitive check
    const chemblIdPattern = /^CHEMBL\d+$/i;
    return chemblIdPattern.test(trimmedId);
};

// In-house ChemBL Similarity Search Type Selector Component with Portal
const ChemBLSimilaritySearchTypeSelector = ({ value, onChange, options, disabled = false }) => {
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
                className="similarity-selector-options-portal"
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
                        className={`similarity-selector-option ${option.value === value ? 'selected' : ''}`}
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
                className={`chembl-similarity-search-selector ${disabled ? 'disabled' : ''}`}
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
    { value: "smiles", label: "SMILES String" },
    { value: "chembl_id", label: "ChEMBL ID" }
];

const ChemBLSimilarityGetter = ({
    toolData = null,
    initialSearchValue = "",
    initialSearchType = "smiles",
    initialSimilarityThreshold = 70,
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
        console.log('similarity tool data received >>', toolData);
        if (toolData && toolData.type && toolData.content) {
            console.log('Processing similarity tool data:', toolData);
            
            // Set flag to indicate we're setting search type from toolData
            isSettingFromToolData.current = true;
            
            // Set search type based on tool data type
            setSearchType(toolData.type);

            console.log('tool data >>', toolData);
            
            if (toolData.type === "smiles" || toolData.type === "chembl_id") {
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
        console.log('similarity api data updated >>', apiData);
        console.log('similar molecules >>', similarMolecules);
        console.log('active similar molecules >>', selectedMolecule);
    }, [apiData]);

    // Backend SMILES validation function (similar to MolExplorer)
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
        if (searchType === "chembl_id") {
            console.log('Validating ChEMBL ID:', trimmedValue);
            isValid = validateChemblId(trimmedValue);
            console.log('ChEMBL ID validation result:', isValid);
            if (!isValid) {
                errorMessage = "Invalid ChEMBL ID format. Must start with 'CHEMBL' followed by numbers (e.g., CHEMBL25, CHEMBL1234567).";
            }
            
            // Update validation states immediately for ChEMBL ID
            setIsValidInput(isValid);
            setValidationStatus(isValid ? 'valid' : 'invalid');
            setValidationError(isValid ? null : errorMessage);
        } else if (searchType === "smiles") {
            // For SMILES, use backend validation like MolExplorer
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
        
        if (searchType === "chembl_id") {
            isCurrentlyValid = validateChemblId(trimmedValue);
            console.log('Immediate ChEMBL ID validation in submit:', isCurrentlyValid);
        } else if (searchType === "smiles") {
            // For SMILES, trust the async validation state since we can't do sync validation
            // But ensure we have a valid state
            isCurrentlyValid = isValidInput && validationStatus === 'valid';
            console.log('SMILES validation state in submit:', isCurrentlyValid, 'validationStatus:', validationStatus);
        }

        if (!isCurrentlyValid) {
            if (searchType === "chembl_id") {
                const errorMsg = trimmedValue.toLowerCase().includes('chembl') 
                    ? "Invalid ChEMBL ID format. Ensure it follows the pattern: CHEMBL followed by numbers (e.g., CHEMBL25, CHEMBL1234567)."
                    : "Invalid ChEMBL ID format. ChEMBL ID must start with 'CHEMBL' followed by numbers (e.g., CHEMBL25, CHEMBL1234567).";
                setError(errorMsg);
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

            if (searchType === "smiles") {
                result = await findSimilarMoleculesBySmiles({ 
                    smiles: trimmedValue,
                    similarity_threshold: similarityThreshold
                });
            } else if (searchType === "chembl_id") {
                result = await findSimilarMoleculesByChemblId({ 
                    chembl_id: trimmedValue,
                    similarity_threshold: similarityThreshold
                });
            }

            if (result) {
                if (result.status === "success") {
                    var processible = result.result["0"];
                    processible = JSON.parse(processible);
                    
                    console.log('Similarity search processible', processible);

                    if (processible.data && Array.isArray(processible.data) && processible.data.length > 0) {
                        // Multiple similar molecules found
                        setSimilarMolecules(processible.data);
                        setSelectedMolecule(processible.data[0]);
                        setSelectedIndex(0);
                        setApiData(null); // Clear single molecule data
                        setSearchValue(trimmedValue);
                    } else if (processible.data && Array.isArray(processible.data) && processible.data.length === 0) {
                        // No similar molecules found
                        setError(`No similar molecules found for "${trimmedValue}" with similarity threshold ${similarityThreshold}%`);
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
            console.error('Similarity search error:', error);
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
            case "smiles":
                return "Enter SMILES string (e.g., CCO, C1=CC=CC=C1)";
            case "chembl_id":
                return "Enter ChEMBL ID (e.g., CHEMBL25)";
            default:
                return "Enter search value";
        }
    };

    const getHeaderText = () => {
        switch (searchType) {
            case "smiles":
                return "SMILES-based Similarity Search";
            case "chembl_id":
                return "ChEMBL ID-based Similarity Search";
            default:
                return "Similarity Search";
        }
    };

    // Reset states when search type changes (but not when set from toolData)
    // useEffect(() => {
    //     // Only reset if the search type change is user-initiated (not from toolData)
    //     if (!isSettingFromToolData.current) {
    //         setError(null);
    //         setValidationError(null);
    //         setInputValue("");
    //         setSearchValue("");
    //         setApiData(null);
    //         setSimilarMolecules([]);
    //         setSelectedMolecule(null);
    //         setSelectedIndex(0);
    //         // Reset validation states for the new search type
    //         setIsValidInput(false);
    //         setValidationStatus('empty');
    //     }
    // }, [searchType]);

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
                className="chembl-similarity-getter-container"
            >
                {!hideInputBox && (
                    <div className="chembl-similarity-getter-row-1">
                        <GlassyContainer>
                            <h3 style={{ marginBottom: '16px', fontWeight: '700' }}>
                                {getHeaderText()}
                            </h3>
                            
                            <div className="chembl-similarity-getter-search-controls">
                                {/* Row 1: Search Type and Similarity Threshold */}
                                <div className="chembl-similarity-getter-controls-row-1">
                                    {/* Search Type Dropdown */}
                                    <div className="chembl-similarity-getter-dropdown-section">
                                        <h4 className="chembl-similarity-getter-dropdown-header">Search Type:</h4>
                                        <div className="chembl-similarity-getter-dropdown">
                                            <ChemBLSimilaritySearchTypeSelector
                                                value={searchType}
                                                onChange={setSearchType}
                                                options={similaritySearchTypeOptions}
                                                disabled={isLoading}
                                            />
                                        </div>
                                    </div>

                                    {/* Similarity Threshold Input */}
                                    <div className="chembl-similarity-getter-threshold-section">
                                        <h4 className="chembl-similarity-getter-dropdown-header">Similarity:</h4>
                                        <input
                                            type="number"
                                            min="0"
                                            max="100"
                                            value={similarityThreshold}
                                            onChange={(e) => setSimilarityThreshold(parseInt(e.target.value) || 70)}
                                            className="chembl-similarity-getter-threshold-input"
                                            disabled={isLoading}
                                            placeholder="70"
                                        />
                                        <span className="chembl-similarity-getter-threshold-help">
                                            %
                                        </span>
                                    </div>
                                </div>

                                {/* Row 2: Search Input with Validation */}
                                <div className="chembl-similarity-getter-controls-row-2">
                                    <div className="chembl-similarity-getter-input-section">
                                        <div className="chembl-similarity-input-row">
                                            <h4 className="chembl-similarity-getter-dropdown-header">Search Value:</h4>
                                            <div className="chembl-similarity-input-with-validation">
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
                                                    className={`chembl-similarity-search-input ${validationStatus}`}
                                                />
                                                <div className="chembl-similarity-validation-icons">
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
                                                className={`chembl-similarity-submit-btn ${!isValidInput || isLoading ? 'disabled' : 'enabled'}`}
                                            >
                                                {isLoading ? 'Searching...' : 'Search'}
                                            </button>
                                        </div>
                                        
                                        {/* Validation Error Message - now properly below the input row */}
                                        {validationError && (
                                            <div className="chembl-similarity-validation-error">
                                                <p>{validationError}</p>
                                            </div>
                                        )}
                                    </div>
                                </div>
                            </div>

                            {/* Error Display */}
                            {error && (
                                <div className="chembl-similarity-getter-error">
                                    <p>{error}</p>
                                </div>
                            )}

                            {/* Loading State */}
                            {isLoading && (
                                <div className="chembl-similarity-getter-loading">
                                    <p>Searching for similar molecules...</p>
                                </div>
                            )}
                        </GlassyContainer>
                    </div>
                )}

                {/* Similar Molecules Results */}
                {((similarMolecules && similarMolecules.length > 0) || error) && (
                    <div className="chembl-similarity-getter-row-2" style={{ marginTop: '16px' }}>
                        <GlassyContainer>
                            <h3 style={{ marginBottom: '12px', fontWeight: '700' }}>
                                {similarMolecules && similarMolecules.length > 0 
                                    ? `Similar Molecules Results (${similarMolecules.length} found)` 
                                    : 'Search Results'
                                }
                            </h3>
                            
                            {/* Show error state when there's an error and no molecules */}
                            {error && (!similarMolecules || similarMolecules.length === 0) ? (
                                <div className="chembl-similarity-getter-error" style={{ margin: '20px 0' }}>
                                    <p>{error}</p>
                                </div>
                            ) : (
                                <div className="chembl-similarity-results-container">
                                {/* Left Panel - Molecule List */}
                                <div className="chembl-similarity-candidates-panel">
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
                                                className={`chembl-similarity-candidate-item ${selectedMolecule === molecule ? 'selected' : ''}`}
                                            >
                                                <span className="chembl-similarity-candidate-number">
                                                    {index + 1}.
                                                </span>
                                                <div className="chembl-similarity-candidate-info">
                                                    <span className="chembl-similarity-candidate-id">
                                                        {molecule.molecule_chembl_id || 'Unknown ID'}
                                                    </span>
                                                    {molecule.similarity_score && (
                                                        <span className="chembl-similarity-candidate-score">
                                                            {(molecule.similarity_score * 100).toFixed(1)}%
                                                        </span>
                                                    )}
                                                </div>
                                            </div>
                                        ))}
                                    </div>
                                </div>

                                {/* Right Panel - Info and 2D Viewer */}
                                <div className="chembl-similarity-details-panel">
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
                                                        {selectedMolecule.molecule_chembl_id || selectedMolecule.pref_name || 'Unknown'}
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
                                            {selectedMolecule.molecule_structures && selectedMolecule.molecule_structures.canonical_smiles ? (
                                                <div className="chembl-similarity-details-row">
                                                    {/* InfoBox */}
                                                    <div style={{ flex: '4' }}>
                                                        <InfoBox 
                                                            activeMol={selectedMolecule.molecule_structures.canonical_smiles} 
                                                            isValidMol={true} 
                                                            infoType={"MOL"} 
                                                        />
                                                    </div>
                                                    
                                                    {/* TwoDViewer */}
                                                    <div style={{ flex: '6' }}>
                                                        <TwoDViewer 
                                                            activeMol={selectedMolecule.molecule_structures.canonical_smiles} 
                                                            isValidMol={true} 
                                                            visType={"MOL"} 
                                                        />
                                                    </div>
                                                </div>
                                            ) : (
                                                <div className="chembl-similarity-placeholder">
                                                    No structure data available for this molecule
                                                </div>
                                            )}
                                        </>
                                    ) : (
                                        <div className="chembl-similarity-placeholder">
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
                    <div className="chembl-similarity-getter-row-3" style={{ marginTop: '20px' }}>
                        <DataViewer
                            data={selectedMolecule}
                            title={`Similar Molecule Data${searchValue ? ` for "${searchValue}"` : ''}${
                                selectedMolecule ? 
                                ` - ${selectedMolecule.molecule_chembl_id || selectedMolecule.pref_name || 'Unknown'}` : 
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

// ChemBLSimilarityGetter.propTypes = {
//     toolData: PropTypes.object,
//     initialSearchValue: PropTypes.string,
//     initialSearchType: PropTypes.oneOf(["smiles", "chembl_id"]),
//     initialSimilarityThreshold: PropTypes.number,
//     hideInputBox: PropTypes.bool,
// };

export default ChemBLSimilarityGetter;
