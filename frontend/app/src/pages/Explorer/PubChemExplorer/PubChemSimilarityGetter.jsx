import React, { useState, useEffect, useCallback, useRef } from "react";
import { createPortal } from "react-dom";
import "./PubChemSimilarityGetter.css";
import SimpleInputBox from "../../../components/UI/SimpleInputBox/SimpleInputBox";
import DataViewer from "../../../components/UI/DataViewer/DataViewer";

import TwoDViewer from "../../../components/UI/TwoDViewer/TwoDViewer";
import { motion, AnimatePresence } from "framer-motion";
import { 
    fadeInUpVariantStatic, 
    fadeInUpVariants, 
    containerVariants, 
    itemVariants,
    scaleInVariants 
} from "../../../components/animations/framerAnim";
import GlassyContainer from "../../../components/glassy_container/gc";
import { fastSimilarity2dSearchByCid, getCompoundByCid } from "../../../services/api/mcpToolsService";
import { FaCheckDouble } from "react-icons/fa6";
import { IoWarningOutline } from "react-icons/io5";
import { call_endpoint_async } from "../../../endpoints/caller";
import { endpoints } from "../../../endpoints/endpoints";

// Enhanced CID Validator
const validateCid = (cid) => {
    // Add extra safety checks
    if (cid === null || cid === undefined || typeof cid !== 'string') {
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
            <motion.div
                ref={dropdownRef}
                className="pubchem-similarity-selector-options-portal"
                initial={{ opacity: 0, y: -10, scale: 0.95 }}
                animate={{ opacity: 1, y: 0, scale: 1 }}
                exit={{ opacity: 0, y: -10, scale: 0.95 }}
                transition={{ duration: 0.2, ease: "easeOut" }}
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
                {options.map((option, index) => (
                    <motion.div
                        key={option.value}
                        className={`pubchem-similarity-selector-option ${option.value === value ? 'selected' : ''}`}
                        onClick={() => handleOptionClick(option.value)}
                        role="option"
                        aria-selected={option.value === value}
                        initial={{ opacity: 0, x: -10 }}
                        animate={{ opacity: 1, x: 0 }}
                        transition={{ 
                            delay: index * 0.05, 
                            duration: 0.2,
                            ease: "easeOut" 
                        }}
                        whileHover={{ 
                            backgroundColor: 'rgba(66, 153, 225, 0.1)',
                            x: 5,
                            transition: { duration: 0.15 }
                        }}
                    >
                        <span className="option-label">{option.label}</span>
                        {option.value === value && (
                            <motion.span 
                                className="check-icon"
                                initial={{ scale: 0, rotate: -90 }}
                                animate={{ scale: 1, rotate: 0 }}
                                transition={{ delay: 0.1, type: "spring", stiffness: 200 }}
                            >
                                <svg width="14" height="10" viewBox="0 0 14 10" fill="currentColor">
                                    <path d="M13 1L5 9L1 5" />
                                </svg>
                            </motion.span>
                        )}
                    </motion.div>
                ))}
            </motion.div>,
            document.body
        );
    };

    return (
        <>
            <div
                className={`pubchem-similarity-search-selector ${disabled ? 'disabled' : ''}`}
                ref={triggerRef}
            >
                <motion.div
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
                    whileHover={{ scale: 1.02 }}
                    whileTap={{ scale: 0.98 }}
                    transition={{ duration: 0.15 }}
                >
                    <span className="selected-text">
                        {selectedOption ? selectedOption.label : 'Select search type'}
                    </span>
                    <motion.span 
                        className={`dropdown-arrow ${isOpen ? 'open' : ''}`}
                        animate={{ rotate: isOpen ? 180 : 0 }}
                        transition={{ duration: 0.2, ease: "easeInOut" }}
                    >
                        <svg
                            width="12"
                            height="8"
                            viewBox="0 0 12 8"
                            fill="currentColor"
                        >
                            <path d="M6 8L0 0h12z" />
                        </svg>
                    </motion.span>
                </motion.div>
            </div>
            <AnimatePresence>
                <PortalDropdown />
            </AnimatePresence>
        </>
    );
};

// Search type options for the dropdown
const similaritySearchTypeOptions = [
    { value: "cid", label: "PubChem CID" },
];

const PubChemSimilarityGetter = ({
    toolData = null,
    initialSearchValue = "",
    initialSearchType = "cid",
    initialSimilarityThreshold = 90,
    hideInputBox = false
}) => {
    const [searchValue, setSearchValue] = useState(initialSearchValue || ""); // Ensure string
    const [inputValue, setInputValue] = useState(initialSearchValue || ""); // Ensure string
    const [searchType, setSearchType] = useState("cid"); // Default to CID
    const [similarityThreshold, setSimilarityThreshold] = useState(initialSimilarityThreshold);
    const [apiData, setApiData] = useState(toolData);
    const [isLoading, setIsLoading] = useState(false);
    const [error, setError] = useState(null);
    const [validationError, setValidationError] = useState(null);
    
    // Enhanced validation states
    const [isValidInput, setIsValidInput] = useState(false);
    const [validationStatus, setValidationStatus] = useState('empty'); // 'empty', 'valid', 'invalid'
    
    // States for multiple similar molecules
    const [similarMolecules, setSimilarMolecules] = useState([]);
    const [selectedMolecule, setSelectedMolecule] = useState(null);

    // Ref to track if search type is being set from toolData
    const isSettingFromToolData = useRef(false);

    useEffect(() => {
        console.log('PubChem similarity tool data received >>', toolData);
        if (toolData && toolData.type && toolData.content) {
            console.log('Processing PubChem similarity tool data:', toolData);
            
            // Set flag to indicate we're setting search type from toolData
            isSettingFromToolData.current = true;
            
            // Set search type based on tool data type
            setSearchType(toolData.type === "smiles" ? "cid" : toolData.type); // Force to CID if smiles

            console.log('tool data >>', toolData);
            
            if (toolData.type === "cid") { // Only handle CID
                // Handle similarity search results
                if (Array.isArray(toolData.content) && toolData.content.length > 0) {
                    console.log('Setting similar molecules from tool data:', toolData.content);
                    setSimilarMolecules(toolData.content);
                    setSelectedMolecule(toolData.content[0]);
                    setApiData(null); // Clear single molecule data for multiple results view
                    // Clear any existing errors for successful data processing
                    setError(null);
                    setValidationError(null);
                } else if (toolData.content) {
                    // If it's a single molecule, treat it as a single-item array for consistent UI
                    setSimilarMolecules([toolData.content]);
                    setSelectedMolecule(toolData.content);
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

    // Enhanced input change handler with comprehensive validation
    const handleInputChange = useCallback(async (value) => {
        console.log('handleInputChange called with:', value, 'for searchType:', searchType);
        
        const currentInputString = typeof value === 'string' ? value : '';
        setInputValue(currentInputString);
        setError(null);
        setValidationError(null);

        // Handle empty input
        if (currentInputString.trim() === "") {
            setValidationError(null);
            setIsValidInput(false);
            setValidationStatus('empty');
            return;
        }

        const trimmedValue = currentInputString.trim(); // Safe to trim now
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
        } 
        // else if (searchType === "smiles") { // REMOVING SMILES VALIDATION LOGIC
            // For SMILES, use backend validation
            // try {
            //     setValidationStatus('validating'); // Show loading state
            //     console.log('Validating SMILES:', trimmedValue);
            //     isValid = await validateSmilesWithBackend(trimmedValue);
            //     console.log('SMILES validation result:', isValid);
                
            //     if (!isValid) {
            //         errorMessage = "Invalid SMILES format. Please enter a valid SMILES string (e.g., CCO, C1=CC=CC=C1).";
            //     }
                
            //     // Update validation states after backend validation
            //     setIsValidInput(isValid);
            //     setValidationStatus(isValid ? 'valid' : 'invalid');
            //     setValidationError(isValid ? null : errorMessage);
            // } catch (error) {
            //     console.error('SMILES validation error:', error);
            //     setIsValidInput(false);
            //     setValidationStatus('invalid');
            //     setValidationError("Error validating SMILES. Please check your input.");
            // }
        // }
    }, [searchType]); // Removed validateSmilesWithBackend from dependencies

    // Re-validate when search type changes - This might be less relevant now or simplify
    useEffect(() => {
        const safeInputValue = inputValue || '';
        if (safeInputValue.trim() !== "") {
            // Force re-validation with the current input for the new search type
            handleInputChange(inputValue);
        } else {
            // If input is empty, ensure validation states are reset properly
            setIsValidInput(false);
            setValidationStatus('empty');
            setValidationError(null);
        }
    }, [searchType, handleInputChange, inputValue]); // Added inputValue to dependencies

    const handleSubmit = useCallback(async (value) => {
        console.log('handleSubmit called with value:', value);
        
        // Ensure value is a string before calling trim
        const valueString = value != null ? String(value) : '';
        const trimmedValue = valueString.trim(); 
        
        console.log("valid state >>", isValidInput);
        
        if (!trimmedValue) {
            setError("Please enter a search value.");
            return;
        }

        // Perform immediate validation as a fallback in case the state hasn't updated yet
        let isCurrentlyValid = isValidInput;
        console.log('validating cid >>', validateCid(trimmedValue)); 
        if (searchType === "cid") {
            isCurrentlyValid = validateCid(trimmedValue);
            console.log('Immediate CID validation in submit:', isCurrentlyValid);
        } 
        // else if (searchType === "smiles") { // REMOVING SMILES SUBMIT LOGIC
            // For SMILES, trust the async validation state since we can't do sync validation
            // But ensure we have a valid state
            // isCurrentlyValid = isValidInput && validationStatus === 'valid';
            // console.log('SMILES validation state in submit:', isCurrentlyValid, 'validationStatus:', validationStatus);
        // }

        if (!isCurrentlyValid) {
            if (searchType === "cid") {
                setError("Invalid CID format. CID must be a positive integer (e.g., 2244, 5281804).");
            } 
            // else if (searchType === "smiles") { // REMOVING SMILES ERROR
                // setError("Invalid SMILES format. Please enter a valid SMILES string (e.g., CCO, C1=CC=CC=C1).");
            // }
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
                    cid: trimmedValue, 
                });
            } 
            // else if (searchType === "smiles") { // REMOVING SMILES API CALL
                // For SMILES, we would need a different API call
                // Since PubChem API doesn't have direct SMILES similarity search in our service,
                // we can show a message or implement it differently
                // setError("SMILES-based similarity search is not yet implemented for PubChem. Please use CID-based search.");
                // setIsLoading(false); // Ensure loading is stopped
                // return;
            // }

            console.log("result >>", result )
            if (result) {
                if (result.status === "success") {
                    var processible = result.result["0"];
                    processible = JSON.parse(processible);
                    
                    console.log('PubChem similarity search processible', processible);

                    if (processible.result) {
                        // Check if result is an array (multiple similar compounds)
                        if (Array.isArray(processible.result)) {
                            if (processible.result.length > 0) {
                                // For each similar CID, get the full compound data
                                const similarCompounds = [];
                                
                                for (const cidData of processible.result.slice(0, 20)) { // Limit to first 20 for performance
                                    try {
                                        const compoundResult = await getCompoundByCid({ cid: cidData.cid });
                                        if (compoundResult && compoundResult.status === "success") {
                                            const compoundData = JSON.parse(compoundResult.result["0"]);
                                            if (compoundData.result) {
                                                // Add similarity score to the compound data
                                                compoundData.result.similarity_score = cidData.similarity;
                                                compoundData.result.similar_cid = cidData.cid;
                                                similarCompounds.push(compoundData.result);
                                            }
                                        }
                                    } catch (error) {
                                        console.warn(`Failed to get compound data for CID ${cidData.cid}:`, error);
                                    }
                                }

                                if (similarCompounds.length > 0) {
                                    setSimilarMolecules(similarCompounds);
                                    setSelectedMolecule(similarCompounds[0]);
                                    setApiData(null); // Clear single molecule data
                                    setSearchValue(trimmedValue);
                                } else {
                                    setError(`No compound data found for similar molecules of "${trimmedValue}"`);
                                    setSimilarMolecules([]);
                                    setSelectedMolecule(null);
                                    setApiData(null);
                                }
                            } else {
                                // Empty array - no similar molecules found
                                setError(`No similar molecules found for CID "${trimmedValue}" with similarity threshold ${similarityThreshold}%`);
                                setSimilarMolecules([]);
                                setSelectedMolecule(null);
                                setApiData(null);
                            }
                        } else {
                            // Result is an object - use it directly as similarity data
                            console.log('Processing object result for similarity search:', processible.result);
                            setSimilarMolecules([processible.result]);
                            setSelectedMolecule(processible.result);
                            setApiData(null); // Clear single molecule data
                            setSearchValue(trimmedValue);
                        }
                    } else {
                        // No data found
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
        } finally {
            setIsLoading(false); // Ensure loading is stopped in finally
        }
    }, [searchType, similarityThreshold, isValidInput, validationStatus]);

    // Effect to handle initial search value
    useEffect(() => {
        if (initialSearchValue && !toolData) {
            // Automatically submit if there's an initial search value and no toolData
            // This ensures that if the component is loaded with a search term, it executes
            setInputValue(initialSearchValue); // Set input value for display
            // Validate before submitting
            if (searchType === "cid") {
                if (validateCid(initialSearchValue)) {
                    setIsValidInput(true);
                    setValidationStatus('valid');
                    handleSubmit(initialSearchValue);
                } else {
                    setIsValidInput(false);
                    setValidationStatus('invalid');
                    setValidationError("Initial CID is invalid.");
                }
            }
            // Add SMILES validation here if it were still supported
        }
    }, [initialSearchValue, toolData, handleSubmit, searchType]); // Added searchType

    const getPlaceholderText = () => {
        switch (searchType) {
            case "cid":
                return "Enter PubChem CID (e.g., 2244)";
            // case "smiles": // Removed SMILES
            //     return "Enter SMILES string (e.g., CCO)";
            default:
                return "Enter search term";
        }
    };

    const getHeaderText = () => {
        switch (searchType) {
            case "cid":
                return "PubChem CID Similarity Search";
            // case "smiles": // Removed SMILES
            //     return "SMILES Similarity Search";
            default:
                return "PubChem Similarity Search";
        }
    };

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
        // This effect might be for re-validating or other actions when searchType changes.
        // If inputValue exists, re-validate it.
        if (inputValue) {
            handleInputChange(inputValue);
        }
    }, [searchType]) // Removed handleInputChange from dep array if it causes loops, ensure it's stable or manage carefully

    useEffect(() => {
        console.log('PubChem similarity states updated:');
        console.log('  similarMolecules:', similarMolecules);
        console.log('  selectedMolecule:', selectedMolecule);
        // console.log('  selectedIndex:', selectedIndex); // REMOVED
    }, [similarMolecules, selectedMolecule]);


    return (
        <GlassyContainer>
            <motion.div 
                className="pubchem-similarity-getter"
                variants={containerVariants}
                initial="hidden"
                animate="visible"
                exit="exit"
            >
                <AnimatePresence mode="wait">
                    {!hideInputBox && (
                        <motion.div 
                            variants={itemVariants} 
                            className="input-section"
                            key="input-section"
                        >
                            <motion.div 
                                className="header-title"
                                variants={fadeInUpVariants}
                                custom={0}
                            >
                                {getHeaderText()}
                            </motion.div>
                            <motion.div
                                variants={fadeInUpVariants}
                                custom={1}
                            >
                                <SimpleInputBox
                                    value={inputValue || ''}
                                    onChange={handleInputChange}
                                    onSubmit={() => handleSubmit(inputValue)}
                                    placeholder={getPlaceholderText()}
                                    buttonText="Search Similar"
                                    isLoading={isLoading}
                                    error={validationError}
                                    disabled={isLoading}
                                />
                            </motion.div>
                        </motion.div>
                    )}
                </AnimatePresence>

                <AnimatePresence mode="wait">
                    {isLoading && (
                        <motion.div
                            key="loading"
                            variants={scaleInVariants}
                            initial="hidden"
                            animate="visible"
                            exit="exit"
                            style={{
                                padding: '20px',
                                textAlign: 'center',
                                color: 'var(--color-text-primary)'
                            }}
                        >
                            <motion.p
                                animate={{ 
                                    scale: [1, 1.05, 1],
                                    opacity: [0.7, 1, 0.7]
                                }}
                                transition={{ 
                                    duration: 1.5, 
                                    repeat: Infinity,
                                    ease: "easeInOut"
                                }}
                            >
                                🔍 Searching for similar molecules...
                            </motion.p>
                        </motion.div>
                    )}
                </AnimatePresence>

                <AnimatePresence mode="wait">
                    {error && (
                        <motion.div
                            key="error"
                            variants={fadeInUpVariants}
                            initial="hidden"
                            animate="visible"
                            exit="exit"
                            style={{
                                padding: '16px 20px',
                                backgroundColor: 'rgba(220, 38, 38, 0.1)',
                                borderRadius: '8px',
                                border: '1px solid rgba(220, 38, 38, 0.3)',
                                textAlign: 'center',
                                color: 'var(--color-text-primary)',
                                margin: '20px 0'
                            }}
                        >
                            <motion.h4 
                                style={{ 
                                    marginBottom: '8px', 
                                    color: '#dc2626', 
                                    margin: '0 0 8px 0',
                                    fontWeight: '600'
                                }}
                                variants={fadeInUpVariants}
                                custom={0}
                            >
                                ⚠️ Error
                            </motion.h4>
                            <motion.p 
                                style={{ 
                                    color: 'var(--color-text-secondary)', 
                                    margin: 0,
                                    fontSize: '0.9rem'
                                }}
                                variants={fadeInUpVariants}
                                custom={1}
                            >
                                {error}
                            </motion.p>
                        </motion.div>
                    )}
                </AnimatePresence>

                <AnimatePresence mode="wait">
                    {(selectedMolecule || apiData) && !isLoading && !error && !validationError ? (
                        <motion.div 
                            key="results"
                            variants={itemVariants} 
                            className="results-display-unified"
                            initial="hidden"
                            animate="visible"
                            exit="exit"
                        >
                            <motion.div
                                variants={fadeInUpVariants}
                                custom={0}
                            >
                                <DataViewer 
                                    data={selectedMolecule || apiData} 
                                />
                            </motion.div>
                        </motion.div>
                    ) : (
                        !isLoading && !error && !validationError && similarMolecules.length === 0 && (
                            <motion.div
                                key="no-data"
                                variants={fadeInUpVariants}
                                initial="hidden"
                                animate="visible"
                                exit="exit"
                                style={{
                                    padding: '16px 20px',
                                    backgroundColor: 'rgba(66, 153, 225, 0.1)',
                                    borderRadius: '8px',
                                    border: '1px solid rgba(66, 153, 225, 0.3)',
                                    textAlign: 'center',
                                    color: 'var(--color-text-primary)',
                                    margin: '20px 0'
                                }}
                            >
                                <motion.h4 
                                    style={{ 
                                        marginBottom: '8px', 
                                        color: 'var(--color-text-primary)', 
                                        margin: '0 0 8px 0',
                                        fontWeight: '600'
                                    }}
                                    variants={fadeInUpVariants}
                                    custom={0}
                                >
                                    ℹ️ No Data Available
                                </motion.h4>
                                <motion.p 
                                    style={{ 
                                        color: 'var(--color-text-secondary)', 
                                        margin: 0,
                                        fontSize: '0.9rem'
                                    }}
                                    variants={fadeInUpVariants}
                                    custom={1}
                                >
                                    No data to display. Perform a search or check your input.
                                </motion.p>
                            </motion.div>
                        )
                    )}
                </AnimatePresence>

                {/* Similar Molecules Results Section */}
                <AnimatePresence mode="wait">
                    {similarMolecules && similarMolecules.length > 1 && !isLoading && !error && (
                        <motion.div 
                            key="similar-molecules"
                            variants={itemVariants}
                            initial="hidden"
                            animate="visible"
                            exit="exit"
                            style={{ marginTop: '20px' }}
                        >
                            <motion.div
                                variants={containerVariants}
                                style={{
                                    backgroundColor: 'var(--color-bg-secondary)',
                                    borderRadius: '12px',
                                    border: '1px solid var(--c-light-border)',
                                    padding: '20px'
                                }}
                            >
                                <motion.h3 
                                    variants={fadeInUpVariants}
                                    custom={0}
                                    style={{ 
                                        marginBottom: '16px', 
                                        fontWeight: '700',
                                        color: 'var(--color-text-primary)'
                                    }}
                                >
                                    🧪 Similar Molecules Results ({similarMolecules.length} found)
                                </motion.h3>
                                
                                <motion.p 
                                    variants={fadeInUpVariants}
                                    custom={1}
                                    style={{
                                        color: 'var(--color-text-secondary)',
                                        fontSize: '0.9rem',
                                        marginBottom: '16px',
                                        fontStyle: 'italic'
                                    }}
                                >
                                    💡 Click on any molecule to view detailed information
                                </motion.p>

                                <div
                                    style={{
                                        display: 'grid',
                                        gridTemplateColumns: 'repeat(auto-fill, minmax(280px, 1fr))',
                                        gap: '16px',
                                        marginTop: '16px'
                                    }}
                                >
                                    {similarMolecules.map((molecule, index) => (
                                        <motion.div
                                            key={`similar-${molecule.cid || molecule.similar_cid || index}`}
                                            variants={fadeInUpVariants}
                                            custom={index}
                                            whileHover={{ 
                                                scale: 1.02, 
                                                y: -5,
                                                transition: { duration: 0.2 }
                                            }}
                                            whileTap={{ scale: 0.98 }}
                                            onClick={() => setSelectedMolecule(molecule)}
                                            style={{
                                                background: selectedMolecule === molecule 
                                                    ? 'var(--color-accent)' 
                                                    : 'var(--color-bg-primary)',
                                                border: selectedMolecule === molecule 
                                                    ? '2px solid var(--color-accent-light)' 
                                                    : '1px solid var(--c-light-border)',
                                                borderRadius: '12px',
                                                padding: '16px',
                                                cursor: 'pointer',
                                                transition: 'all 0.3s ease',
                                                position: 'relative',
                                                overflow: 'hidden'
                                            }}
                                        >
                                            {/* Selection indicator */}
                                            {selectedMolecule === molecule && (
                                                <motion.div
                                                    initial={{ scale: 0, rotate: -90 }}
                                                    animate={{ scale: 1, rotate: 0 }}
                                                    transition={{ delay: 0.1, type: "spring", stiffness: 200 }}
                                                    style={{
                                                        position: 'absolute',
                                                        top: '12px',
                                                        right: '12px',
                                                        color: 'white',
                                                        fontSize: '16px',
                                                        fontWeight: 'bold'
                                                    }}
                                                >
                                                    ✓
                                                </motion.div>
                                            )}

                                            <div style={{ marginBottom: '12px' }}>
                                                <motion.div 
                                                    style={{
                                                        color: selectedMolecule === molecule 
                                                            ? 'white' 
                                                            : 'var(--color-accent)',
                                                        fontSize: '14px',
                                                        fontWeight: 'bold',
                                                        marginBottom: '4px'
                                                    }}
                                                    variants={fadeInUpVariants}
                                                    custom={0}
                                                >
                                                    CID: {molecule.cid || molecule.similar_cid || 'Unknown'}
                                                </motion.div>
                                                
                                                {molecule.similarity_score && (
                                                    <motion.div 
                                                        style={{
                                                            color: selectedMolecule === molecule 
                                                                ? 'rgba(255, 255, 255, 0.9)' 
                                                                : 'var(--color-success)',
                                                            fontSize: '12px',
                                                            fontWeight: '600'
                                                        }}
                                                        variants={fadeInUpVariants}
                                                        custom={1}
                                                    >
                                                        {(molecule.similarity_score * 100).toFixed(1)}% similar
                                                    </motion.div>
                                                )}
                                            </div>

                                            {molecule.compound_name && (
                                                <motion.div 
                                                    style={{
                                                        color: selectedMolecule === molecule 
                                                            ? 'rgba(255, 255, 255, 0.9)' 
                                                            : 'var(--color-text-primary)',
                                                        fontSize: '13px',
                                                        marginBottom: '12px',
                                                        fontWeight: '500',
                                                        lineHeight: '1.3'
                                                    }}
                                                    variants={fadeInUpVariants}
                                                    custom={2}
                                                >
                                                    {molecule.compound_name.slice(0, 80)}
                                                    {molecule.compound_name.length > 80 && '...'}
                                                </motion.div>
                                            )}

                                            {(molecule.cid || molecule.similar_cid) && (
                                                <motion.div 
                                                    style={{ marginTop: '12px' }}
                                                    variants={scaleInVariants}
                                                    custom={index}
                                                >
                                                    <img 
                                                        src={`https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/${molecule.cid || molecule.similar_cid}/PNG?record_type=2d&image_size=small`}
                                                        alt={`Compound ${molecule.cid || molecule.similar_cid}`}
                                                        style={{
                                                            width: '100%',
                                                            height: 'auto',
                                                            borderRadius: '8px',
                                                            backgroundColor: 'white',
                                                            padding: '8px',
                                                            border: selectedMolecule === molecule 
                                                                ? '2px solid rgba(255, 255, 255, 0.5)' 
                                                                : '1px solid var(--c-light-border)'
                                                        }}
                                                        onError={(e) => {
                                                            e.target.style.display = 'none';
                                                        }}
                                                    />
                                                </motion.div>
                                            )}
                                        </motion.div>
                                    ))}
                                </div>
                            </motion.div>
                        </motion.div>
                    )}
                </AnimatePresence>
            </motion.div>
        </GlassyContainer>
    );
};

export default PubChemSimilarityGetter;
