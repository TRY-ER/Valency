import React, { useState, useEffect, useCallback, useRef } from "react";
import { createPortal } from "react-dom";
// import PropTypes from "prop-types";
import "./ChemBLGetter.css";
import SimpleInputBox from "../../../components/UI/SimpleInputBox/SimpleInputBox";
import DataViewer from "../../../components/UI/DataViewer/DataViewer";
import InfoBox from "../../../components/UI/InfoBox/InfoBox";
import TwoDViewer from "../../../components/UI/TwoDViewer/TwoDViewer";
import { motion } from "framer-motion";
import { fadeInUpVariantStatic } from "../../../components/animations/framerAnim";
import GlassyContainer from "../../../components/glassy_container/gc";
import { getMoleculeByChemblId, findMoleculeByPrefName, findMoleculeBySynonym, getMoleculesByChemblIds } from "../../../services/api/mcpToolsService";

// ChEMBL ID Validator
const validateChemblId = (chemblId) => {
    // ChEMBL ID must start with "CHEMBL" followed by one or more digits
    const chemblIdPattern = /^CHEMBL\d+$/i;
    return chemblIdPattern.test(chemblId.trim());
};

// In-house ChemBL Search Type Selector Component with Portal
const ChemBLSearchTypeSelector = ({ value, onChange, options, disabled = false }) => {
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

    // Portal dropdown component
    const PortalDropdown = () => {
        if (!isOpen || disabled) return null;

        return createPortal(
            <div
                ref={dropdownRef}
                className="selector-options-portal"
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
                        className={`selector-option ${option.value === value ? 'selected' : ''}`}
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
                className={`chembl-search-selector ${disabled ? 'disabled' : ''}`}
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

// Search type options for the dropdown (moved outside component to prevent re-creation)
const searchTypeOptions = [
    { value: "chembl_id", label: "ChEMBL ID" },
    { value: "multiple_chembl_ids", label: "Multiple ChEMBL IDs" },
    { value: "pref_name", label: "Preferred Name" },
    { value: "synonym", label: "Synonym" }
];

const ChemBLGetter = ({
    toolData = null,
    initialSearchValue = "",
    initialSearchType = "chembl_id",
    hideInputBox = false
}) => {
    const [searchValue, setSearchValue] = useState(initialSearchValue);
    const [inputValue, setInputValue] = useState(initialSearchValue);
    const [searchType, setSearchType] = useState(initialSearchType);
    const [apiData, setApiData] = useState(toolData);
    const [isLoading, setIsLoading] = useState(false);
    const [error, setError] = useState(null);
    const [validationError, setValidationError] = useState(null);
    
    // States for multiple candidates (preferred name search)
    const [multipleCandidates, setMultipleCandidates] = useState([]);
    const [selectedCandidate, setSelectedCandidate] = useState(null);
    const [selectedIndex, setSelectedIndex] = useState(0);
    const [isSelectedCandidateValid, setIsSelectedCandidateValid] = useState(true);

    // States for multiple ChEMBL IDs handling
    const [currentChemblId, setCurrentChemblId] = useState("");
    const [isValidCurrentChemblId, setIsValidCurrentChemblId] = useState(false);
    const [chemblIdsList, setChemblIdsList] = useState([]);
    const [multipleChemblResults, setMultipleChemblResults] = useState([]);
    const [selectedChemblResult, setSelectedChemblResult] = useState(null);
    const [selectedChemblIndex, setSelectedChemblIndex] = useState(0);

    useEffect(() => {
        console.log('tool data received >>', toolData);
    }, [toolData]);

    useEffect(() => {
        console.log("apiData updated >>", apiData);
    }, [apiData]);

    // Handle input changes with real-time validation
    const handleInputChange = useCallback((value) => {
        setInputValue(value);

        // Clear previous validation errors
        setValidationError(null);

        // Only validate if there's a value and search type is chembl_id
        if (value && value.trim() && (searchType === "chembl_id" || searchType === "multiple_chembl_ids")) {
            if (!validateChemblId(value.trim())) {
                setValidationError("Invalid ChEMBL ID format. Must start with 'CHEMBL' followed by numbers.");
            }
        }

        // For multiple ChEMBL IDs, update current ChEMBL ID state
        if (searchType === "multiple_chembl_ids") {
            setCurrentChemblId(value);
            setIsValidCurrentChemblId(value && value.trim() && validateChemblId(value.trim()));
        }
    }, [searchType]);

    const handleSubmit = useCallback(async (value) => {
        // For multiple ChEMBL IDs, handle different logic
        if (searchType === "multiple_chembl_ids") {
            if (chemblIdsList.length === 0) {
                setError("Please add at least one ChEMBL ID to the list.");
                return;
            }

            setIsLoading(true);
            setError(null);
            setValidationError(null);

            try {
                const result = await getMoleculesByChemblIds({ chembl_ids: chemblIdsList });
                
                if (result && result.status === "success") {
                    const processible = JSON.parse(result.result["0"]);
                    console.log('Multiple ChEMBL results:', processible);

                    if (processible.result && Array.isArray(processible.result) && processible.result.length > 0) {
                        setMultipleChemblResults(processible.result);
                        setSelectedChemblResult(processible.result[0]);
                        setSelectedChemblIndex(0);
                        setApiData(null); // Clear single molecule data
                        setMultipleCandidates([]);
                        setSelectedCandidate(null);
                    } else {
                        setError("No molecules found for the provided ChEMBL IDs.");
                        setMultipleChemblResults([]);
                        setSelectedChemblResult(null);
                    }
                } else {
                    setError(result?.message || "Failed to fetch molecules for the provided ChEMBL IDs.");
                    setMultipleChemblResults([]);
                    setSelectedChemblResult(null);
                }
            } catch (error) {
                console.error('Error fetching multiple ChEMBL molecules:', error);
                setError(error.message || "Failed to fetch ChEMBL molecules. Please try again.");
                setMultipleChemblResults([]);
                setSelectedChemblResult(null);
            } finally {
                setIsLoading(false);
            }
            return;
        }

        // Original logic for other search types
        if (!value || value.trim() === "") {
            setError("Please enter a search value.");
            return;
        }

        const trimmedValue = value.trim();

        // Validate ChEMBL ID format if search type is chembl_id
        if (searchType === "chembl_id") {
            if (!validateChemblId(trimmedValue)) {
                setError("Invalid ChEMBL ID format. ChEMBL ID must start with 'CHEMBL' followed by numbers (e.g., CHEMBL25, CHEMBL1234567).");
                return;
            }
        }

        setIsLoading(true);
        setError(null);
        setValidationError(null);

        try {
            let result;

            if (searchType === "chembl_id") {
                result = await getMoleculeByChemblId({ chembl_id: trimmedValue });
            } else if (searchType === "pref_name") {
                result = await findMoleculeByPrefName({ pref_name: trimmedValue });
            } else if (searchType === "synonym") {
                result = await findMoleculeBySynonym({ synonym: trimmedValue });
            }

            if (result) {
                if (result.status === "success") {
                    var processible = result.result["0"]
                    processible = JSON.parse(processible);
                    
                    console.log('Processible', processible);

                    if (searchType === "chembl_id") {
                        // For ChEMBL ID search, we expect a single molecule
                        if (processible.result){
                            setApiData(processible.result);
                            setSearchValue(trimmedValue);
                            // Clear multiple candidates for single molecule view
                            setMultipleCandidates([]);
                            setSelectedCandidate(null);
                            setMultipleChemblResults([]);
                            setSelectedChemblResult(null);
                        } else {
                            // No data found for ChEMBL ID
                            setError(`No molecule found for ChEMBL ID "${trimmedValue}"`);
                            setApiData(null);
                            setMultipleCandidates([]);
                            setSelectedCandidate(null);
                            setMultipleChemblResults([]);
                            setSelectedChemblResult(null);
                        }
                    } else if (searchType === "pref_name") {
                        // For preferred name search, we might get multiple candidates
                        if (processible.result && Array.isArray(processible.result) && processible.result.length > 0) {
                            // Multiple candidates found
                            setMultipleCandidates(processible.result);
                            setSelectedCandidate(processible.result[0]);
                            setSelectedIndex(0);
                            setIsSelectedCandidateValid(true);
                            setApiData(null); // Clear single molecule data
                            setSearchValue(trimmedValue);
                            setMultipleChemblResults([]);
                            setSelectedChemblResult(null);
                        } else if (processible.result && Array.isArray(processible.result) && processible.result.length === 0) {
                            // Empty array - no records found
                            console.log(`No records found for preference search "${trimmedValue}"`);
                            setError(`No records found with preference search "${trimmedValue}"`);
                            setMultipleCandidates([]);
                            setSelectedCandidate(null);
                            setApiData(null);
                            setMultipleChemblResults([]);
                            setSelectedChemblResult(null);
                            // Don't set searchValue when there's an error to ensure error displays properly
                        } else if (processible.result) {
                            // Single candidate found (non-array data)
                            setApiData(processible.result);
                            setSearchValue(trimmedValue);
                            setMultipleCandidates([]);
                            setSelectedCandidate(null);
                            setMultipleChemblResults([]);
                            setSelectedChemblResult(null);
                        } else {
                            // No data found
                            setError(`No records found with preference search "${trimmedValue}"`);
                            setMultipleCandidates([]);
                            setSelectedCandidate(null);
                            setApiData(null);
                            setMultipleChemblResults([]);
                            setSelectedChemblResult(null);
                            // Don't set searchValue when there's an error to ensure error displays properly
                        }
                    } else if (searchType === "synonym") {
                        // For synonym search, we might get multiple candidates
                        if (processible.result && Array.isArray(processible.result) && processible.result.length > 0) {
                            // Multiple candidates found
                            setMultipleCandidates(processible.result);
                            setSelectedCandidate(processible.result[0]);
                            setSelectedIndex(0);
                            setIsSelectedCandidateValid(true);
                            setApiData(null); // Clear single molecule data
                            setSearchValue(trimmedValue);
                            setMultipleChemblResults([]);
                            setSelectedChemblResult(null);
                        } else if (processible.result && Array.isArray(processible.result) && processible.result.length === 0) {
                            // Empty array - no records found
                            console.log(`No records found for synonym search "${trimmedValue}"`);
                            setError(`No records found with synonym search "${trimmedValue}"`);
                            setMultipleCandidates([]);
                            setSelectedCandidate(null);
                            setApiData(null);
                            setMultipleChemblResults([]);
                            setSelectedChemblResult(null);
                            // Don't set searchValue when there's an error to ensure error displays properly
                        } else if (processible.result) {
                            // Single candidate found (non-array data)
                            setApiData(processible.result);
                            setSearchValue(trimmedValue);
                            setMultipleCandidates([]);
                            setSelectedCandidate(null);
                            setMultipleChemblResults([]);
                            setSelectedChemblResult(null);
                        } else {
                            // No data found
                            setError(`No records found with synonym search "${trimmedValue}"`);
                            setMultipleCandidates([]);
                            setSelectedCandidate(null);
                            setApiData(null);
                            setMultipleChemblResults([]);
                            setSelectedChemblResult(null);
                            // Don't set searchValue when there's an error to ensure error displays properly
                        }
                    }
                } else {
                    // API returned non-success status
                    setError(result.message || `Search failed for "${trimmedValue}". Please try again.`);
                    setApiData(null);
                    setMultipleCandidates([]);
                    setSelectedCandidate(null);
                    setMultipleChemblResults([]);
                    setSelectedChemblResult(null);
                }
            } else {
                setError("No data found for the provided search value.");
                setApiData(null);
                setMultipleCandidates([]);
                setSelectedCandidate(null);
                setMultipleChemblResults([]);
                setSelectedChemblResult(null);
            }
        } catch (error) {
            console.error('Error fetching ChEMBL data:', error);
            setError(error.message || "Failed to fetch ChEMBL data. Please try again.");
        } finally {
            setIsLoading(false);
        }
    }, [searchType, chemblIdsList]);

    // Effect to handle initial data passed as props
    useEffect(() => {
        if (toolData && toolData.type && toolData.content) {
            console.log('Processing tool data:', toolData);
            
            // Set search type based on tool data type
            setSearchType(toolData.type);
            
            if (toolData.type === "pref_name" || toolData.type === "synonym") {
                // Handle preferred name and synonym search results
                if (Array.isArray(toolData.content) && toolData.content.length > 0) {
                    setMultipleCandidates(toolData.content);
                    setSelectedCandidate(toolData.content[0]);
                    setSelectedIndex(0);
                    setIsSelectedCandidateValid(true);
                    setApiData(null); // Clear single molecule data
                    setMultipleChemblResults([]);
                    setSelectedChemblResult(null);
                } else if (toolData.content) {
                    // If it's a single candidate, treat it as a single-item array for consistent UI
                    setMultipleCandidates([toolData.content]);
                    setSelectedCandidate(toolData.content);
                    setSelectedIndex(0);
                    setIsSelectedCandidateValid(true);
                    setApiData(null); // Clear single molecule data
                    setMultipleChemblResults([]);
                    setSelectedChemblResult(null);
                } else {
                    // No candidates found
                    setError(`No records found with ${toolData.type} search.`);
                }
            } else if (toolData.type === "chembl_id") {
                // Handle single ChEMBL ID search result
                if (toolData.content) {
                    setApiData(toolData.content);
                    setMultipleCandidates([]);
                    setSelectedCandidate(null);
                    setMultipleChemblResults([]);
                    setSelectedChemblResult(null);
                } else {
                    setError("No molecule found for the provided ChEMBL ID.");
                }
            } else if (toolData.type === "multiple_chembl_ids") {
                // Handle multiple ChEMBL IDs search result
                if (Array.isArray(toolData.content) && toolData.content.length > 0) {
                    setMultipleChemblResults(toolData.content);
                    setSelectedChemblResult(toolData.content[0]);
                    setSelectedChemblIndex(0);
                    setApiData(null);
                    setMultipleCandidates([]);
                    setSelectedCandidate(null);
                } else {
                    setError("No molecules found for the provided ChEMBL IDs.");
                }
            }
            
            // Clear any existing errors
            setError(null);
            setValidationError(null);
        }
    }, [toolData]);

    // Effect to handle initial search value
    useEffect(() => {
        if (initialSearchValue && initialSearchValue.trim() !== "") {
            setSearchValue(initialSearchValue);
            setInputValue(initialSearchValue);
            if (!toolData) {
                handleSubmit(initialSearchValue);
            }
        }
    }, [initialSearchValue, toolData, handleSubmit]);

    const getPlaceholderText = () => {
        switch (searchType) {
            case "chembl_id":
                return "Format: CHEMBL + numbers (e.g., CHEMBL25)";
            case "multiple_chembl_ids":
                return "Add ChEMBL IDs to list (e.g., CHEMBL25)";
            case "pref_name":
                return "Try name: Aspirin";
            case "synonym":
                return "Try synonym: ACIDUM ACETYLSALICYLICUM";
            default:
                return "Enter search value";
        }
    };

    const getHeaderText = () => {
        switch (searchType) {
            case "chembl_id":
                return "Enter ChEMBL ID";
            case "multiple_chembl_ids":
                return "Add ChEMBL ID to List";
            case "pref_name":
                return "Enter Molecule Preferred Name";
            case "synonym":
                return "Enter Molecule Synonym";
            default:
                return "Enter Search Value";
        }
    };

    useEffect(() => {
        console.log("Search type changed:", searchType);
        // Clear validation errors when search type changes, but keep existing errors
        setValidationError(null);
        // Only clear error if user is actively changing search type (not on initial mount)
    }, [searchType])

    // Keyboard navigation for multiple candidates
    useEffect(() => {
        const handleKeyDown = (event) => {
            if (multipleCandidates.length > 0) {
                if (event.key === 'ArrowUp') {
                    event.preventDefault();
                    setSelectedIndex(prev => {
                        const newIndex = prev > 0 ? prev - 1 : multipleCandidates.length - 1;
                        setSelectedCandidate(multipleCandidates[newIndex]);
                        return newIndex;
                    });
                } else if (event.key === 'ArrowDown') {
                    event.preventDefault();
                    setSelectedIndex(prev => {
                        const newIndex = prev < multipleCandidates.length - 1 ? prev + 1 : 0;
                        setSelectedCandidate(multipleCandidates[newIndex]);
                        return newIndex;
                    });
                }
            }
            
            // Keyboard navigation for multiple ChEMBL results
            if (multipleChemblResults.length > 0) {
                if (event.key === 'ArrowUp') {
                    event.preventDefault();
                    setSelectedChemblIndex(prev => {
                        const newIndex = prev > 0 ? prev - 1 : multipleChemblResults.length - 1;
                        setSelectedChemblResult(multipleChemblResults[newIndex]);
                        return newIndex;
                    });
                } else if (event.key === 'ArrowDown') {
                    event.preventDefault();
                    setSelectedChemblIndex(prev => {
                        const newIndex = prev < multipleChemblResults.length - 1 ? prev + 1 : 0;
                        setSelectedChemblResult(multipleChemblResults[newIndex]);
                        return newIndex;
                    });
                }
            }
        };

        window.addEventListener('keydown', handleKeyDown);
        return () => window.removeEventListener('keydown', handleKeyDown);
    }, [multipleCandidates, multipleChemblResults]);

    // Handlers for multiple ChEMBL IDs functionality
    const handleAddChemblId = () => {
        if (isValidCurrentChemblId && currentChemblId.trim() !== "") {
            const trimmedId = currentChemblId.trim();
            if (!chemblIdsList.includes(trimmedId)) {
                setChemblIdsList(prevList => [...prevList, trimmedId]);
                setCurrentChemblId("");
                setIsValidCurrentChemblId(false);
                setInputValue(""); // Clear the input box
            } else {
                setValidationError("This ChEMBL ID is already in the list.");
            }
        } else {
            setValidationError("Please enter a valid ChEMBL ID before adding to the list.");
        }
    };

    const handleRemoveChemblId = (indexToRemove) => {
        setChemblIdsList(prevList => prevList.filter((_, index) => index !== indexToRemove));
    };

    const handleProcessChemblList = async () => {
        if (chemblIdsList.length === 0) {
            setError("Please add at least one ChEMBL ID to the list.");
            return;
        }
        
        // Trigger the handleSubmit with a dummy value to process the list
        await handleSubmit("process_list");
    };

    const handleResetChemblList = () => {
        setCurrentChemblId("");
        setIsValidCurrentChemblId(false);
        setChemblIdsList([]);
        setMultipleChemblResults([]);
        setSelectedChemblResult(null);
        setSelectedChemblIndex(0);
        setInputValue("");
        setError(null);
        setValidationError(null);
    };

    useEffect(() => {
        console.log("Error state changed:", error);    
    }, [error])

    return (
        <>
            <motion.div
                initial="hidden"
                animate="visible"
                variants={fadeInUpVariantStatic}
                className="chembl-getter-container"
            >
                {!hideInputBox && (
                    <div className="chembl-getter-row-1">
                        <GlassyContainer>
                            <div className="chembl-getter-search-controls">
                                <div className="chembl-getter-dropdown-section">
                                    <h3 className="chembl-getter-dropdown-header">Search Type</h3>
                                    <ChemBLSearchTypeSelector
                                        value={searchType}
                                        onChange={(e) => setSearchType(e.target.value)}
                                        options={searchTypeOptions}
                                        disabled={isLoading}
                                    />
                                </div>

                                <div className="chembl-getter-input-section">
                                    {searchType === "multiple_chembl_ids" ? (
                                        <div style={{ display: 'flex', flexDirection: 'column', gap: '16px' }}>
                                            <SimpleInputBox
                                                value={inputValue}
                                                onChange={handleInputChange}
                                                onSubmit={(value) => {
                                                    if (value && value.trim() && validateChemblId(value.trim())) {
                                                        handleAddChemblId();
                                                    } else {
                                                        setValidationError("Please enter a valid ChEMBL ID.");
                                                    }
                                                }}
                                                header={getHeaderText()}
                                                placeholder={getPlaceholderText()}
                                                buttonText="Add to List"
                                                isLoading={false}
                                            />
                                            
                                            {/* ChEMBL IDs List Display */}
                                            {chemblIdsList.length > 0 && (
                                                <div style={{
                                                    padding: '16px',
                                                    border: '1px solid var(--c-light-border)',
                                                    borderRadius: '8px',
                                                    backgroundColor: 'var(--color-bg-secondary)'
                                                }}>
                                                    <h4 style={{
                                                        margin: '0 0 12px 0',
                                                        fontSize: '0.9rem',
                                                        color: 'var(--color-text-primary)',
                                                        fontWeight: '500'
                                                    }}>
                                                        ChEMBL IDs List ({chemblIdsList.length})
                                                    </h4>
                                                    <div style={{
                                                        display: 'flex',
                                                        flexWrap: 'wrap',
                                                        gap: '8px',
                                                        marginBottom: '12px'
                                                    }}>
                                                        {chemblIdsList.map((id, index) => (
                                                            <div
                                                                key={index}
                                                                style={{
                                                                    display: 'flex',
                                                                    alignItems: 'center',
                                                                    padding: '6px 12px',
                                                                    backgroundColor: 'var(--color-accent)',
                                                                    color: 'white',
                                                                    borderRadius: '20px',
                                                                    fontSize: '0.8rem',
                                                                    fontFamily: 'monospace',
                                                                    gap: '8px'
                                                                }}
                                                            >
                                                                <span>{id}</span>
                                                                <button
                                                                    onClick={() => handleRemoveChemblId(index)}
                                                                    style={{
                                                                        background: 'none',
                                                                        border: 'none',
                                                                        color: 'white',
                                                                        cursor: 'pointer',
                                                                        fontSize: '12px',
                                                                        padding: '0',
                                                                        width: '16px',
                                                                        height: '16px',
                                                                        borderRadius: '50%',
                                                                        display: 'flex',
                                                                        alignItems: 'center',
                                                                        justifyContent: 'center'
                                                                    }}
                                                                    title="Remove"
                                                                >
                                                                    ×
                                                                </button>
                                                            </div>
                                                        ))}
                                                    </div>
                                                    <div style={{ display: 'flex', gap: '8px' }}>
                                                        <button
                                                            onClick={handleProcessChemblList}
                                                            disabled={isLoading || chemblIdsList.length === 0}
                                                            style={{
                                                                padding: '8px 16px',
                                                                backgroundColor: 'var(--color-accent)',
                                                                color: 'white',
                                                                border: 'none',
                                                                borderRadius: '6px',
                                                                cursor: isLoading || chemblIdsList.length === 0 ? 'not-allowed' : 'pointer',
                                                                fontSize: '0.9rem',
                                                                fontWeight: '500',
                                                                opacity: isLoading || chemblIdsList.length === 0 ? 0.6 : 1
                                                            }}
                                                        >
                                                            {isLoading ? 'Processing...' : 'Fetch All Molecules'}
                                                        </button>
                                                        <button
                                                            onClick={handleResetChemblList}
                                                            disabled={isLoading}
                                                            style={{
                                                                padding: '8px 16px',
                                                                backgroundColor: 'var(--color-bg-tertiary)',
                                                                color: 'var(--color-text-primary)',
                                                                border: '1px solid var(--c-light-border)',
                                                                borderRadius: '6px',
                                                                cursor: isLoading ? 'not-allowed' : 'pointer',
                                                                fontSize: '0.9rem',
                                                                fontWeight: '500',
                                                                opacity: isLoading ? 0.6 : 1
                                                            }}
                                                        >
                                                            Reset List
                                                        </button>
                                                    </div>
                                                </div>
                                            )}
                                        </div>
                                    ) : (
                                        <SimpleInputBox
                                            value={inputValue}
                                            onChange={handleInputChange}
                                            onSubmit={handleSubmit}
                                            header={getHeaderText()}
                                            placeholder={getPlaceholderText()}
                                            buttonText="Search Molecule"
                                            isLoading={isLoading}
                                            // error={error || validationError}
                                        />
                                    )}
                                </div>
                            </div>
                        </GlassyContainer>
                    </div>
                )}
                {(error || validationError) && !isLoading && (
                    <div className="chembl-getter-error-section" style={{ marginTop: '20px' }}>
                        <GlassyContainer>
                            <div style={{
                                padding: '16px',
                                borderRadius: '8px',
                                backgroundColor: 'rgba(255, 59, 48, 0.1)',
                                border: '1px solid rgba(255, 59, 48, 0.3)',
                                display: 'flex',
                                alignItems: 'center',
                                gap: '12px'
                            }}>
                                <div style={{
                                    color: '#ff3b30',
                                    fontSize: '20px',
                                    minWidth: '20px'
                                }}>
                                    ⚠️
                                </div>
                                <div>
                                    <h4 style={{
                                        color: '#ff3b30',
                                        margin: '0 0 4px 0',
                                        fontSize: '1rem',
                                        fontWeight: '600'
                                    }}>
                                        Search Error
                                    </h4>
                                    <p style={{
                                        color: '#ff3b30',
                                        margin: '0',
                                        fontSize: '0.9rem',
                                        lineHeight: '1.4'
                                    }}>
                                        {error || validationError}
                                    </p>
                                </div>
                            </div>
                        </GlassyContainer>
                    </div>
                )}

                {isLoading && (
                    <div className="chembl-getter-loading" style={{ marginTop: '20px', textAlign: 'center' }}>
                        <GlassyContainer>
                            <div style={{ padding: '20px' }}>
                                <div className="chembl-getter-loading-spinner" style={{
                                    border: '4px solid #f3f3f3',
                                    borderTop: '4px solid rgb(3, 196, 3)',
                                    borderRadius: '50%',
                                    width: '40px',
                                    height: '40px',
                                    animation: 'spin 2s linear infinite',
                                    margin: '0 auto'
                                }}></div>
                                <p style={{ marginTop: '10px' }}>Searching ChEMBL database...</p>
                            </div>
                        </GlassyContainer>
                    </div>
                )}

                {/* Single molecule visualization for ChEMBL ID search */}
                {apiData && apiData.molecule_structures && apiData.molecule_structures.canonical_smiles && searchType === "chembl_id" && (
                    <div className="chembl-getter-row-2" style={{ marginTop: '20px' }}>
                        <GlassyContainer>
                            <h3 style={{ marginBottom: '15px', fontWeight: '700' }}>Molecule Structure Visualization</h3>
                            <div style={{ display: 'flex', gap: '20px' }}>
                                <div style={{ flex: '1' }}>
                                    <InfoBox 
                                        activeMol={apiData.molecule_structures.canonical_smiles} 
                                        isValidMol={true} 
                                        infoType={"MOL"} 
                                    />
                                </div>
                                <div style={{ flex: '1' }}>
                                    <TwoDViewer 
                                        activeMol={apiData.molecule_structures.canonical_smiles} 
                                        isValidMol={true} 
                                        visType={"MOL"} 
                                    />
                                </div>
                            </div>
                        </GlassyContainer>
                    </div>
                )}

                {/* Multiple candidates view for preferred name and synonym search */}
                {multipleCandidates && multipleCandidates.length > 0 && (searchType === "pref_name" || searchType === "synonym") && (
                    <div className="chembl-getter-row-2" style={{ marginTop: '20px' }}>
                        <GlassyContainer>
                            <h3 style={{ marginBottom: '15px', fontWeight: '700' }}>
                                Molecule Candidates {multipleCandidates.length > 0 ? `(${multipleCandidates.length} found)` : ''}
                            </h3>
                            
                            <div className="chembl-results-container">
                                {/* Left Panel - Candidate List */}
                                <div className="chembl-candidates-panel">
                                    <h4 style={{
                                        color: 'var(--color-text-primary)',
                                        fontSize: '1rem',
                                        marginBottom: '10px',
                                        fontWeight: '500'
                                    }}>
                                        Candidate Molecules
                                    </h4>
                                    <p style={{
                                        color: 'var(--color-text-secondary)',
                                        fontSize: '0.8rem',
                                        marginBottom: '12px',
                                        fontStyle: 'italic'
                                    }}>
                                        💡 Use ↑↓ arrow keys to navigate
                                    </p>
                                    <div className="chembl-candidates-list">
                                        {multipleCandidates.map((candidate, index) => (
                                            <div 
                                                key={index}
                                                onClick={() => {
                                                    setSelectedCandidate(candidate);
                                                    setIsSelectedCandidateValid(true);
                                                    setSelectedIndex(index);
                                                }}
                                                className={`chembl-candidate-item ${selectedCandidate === candidate ? 'selected' : ''}`}
                                            >
                                                <span className="chembl-candidate-number">
                                                    {index + 1}.
                                                </span>
                                                <span className="chembl-candidate-smiles">
                                                    {candidate.pref_name || candidate.molecule_chembl_id || 'Unknown'}
                                                </span>
                                            </div>
                                        ))}
                                    </div>
                                </div>

                                {/* Right Panel - Info and 2D Viewer */}
                                <div className="chembl-details-panel">
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
                                                    Selected ({selectedIndex + 1}/{multipleCandidates.length}): <span style={{ 
                                                        fontFamily: 'monospace', 
                                                        color: 'var(--color-accent)',
                                                        fontWeight: '600' 
                                                    }}>
                                                        {selectedCandidate.pref_name || selectedCandidate.molecule_chembl_id || 'Unknown'}
                                                    </span>
                                                </h4>
                                            </div>
                                            
                                            {/* Show visualization only if we have SMILES data */}
                                            {selectedCandidate.molecule_structures && selectedCandidate.molecule_structures.canonical_smiles ? (
                                                <div className="chembl-details-row">
                                                    {/* InfoBox */}
                                                    <div style={{ flex: '4' }}>
                                                        <InfoBox 
                                                            activeMol={selectedCandidate.molecule_structures.canonical_smiles} 
                                                            isValidMol={isSelectedCandidateValid} 
                                                            infoType={"MOL"} 
                                                        />
                                                    </div>
                                                    
                                                    {/* TwoDViewer */}
                                                    <div style={{ flex: '6' }}>
                                                        <TwoDViewer 
                                                            activeMol={selectedCandidate.molecule_structures.canonical_smiles} 
                                                            isValidMol={isSelectedCandidateValid} 
                                                            visType={"MOL"} 
                                                        />
                                                    </div>
                                                </div>
                                            ) : (
                                                <div className="chembl-placeholder">
                                                    No structure data available for this candidate
                                                </div>
                                            )}
                                        </>
                                    ) : (
                                        <div className="chembl-placeholder">
                                            Select a candidate molecule from the left panel to view details
                                        </div>
                                    )}
                                </div>
                            </div>
                        </GlassyContainer>
                    </div>
                )}

                {/* Multiple ChEMBL results view for multiple_chembl_ids search */}
                {multipleChemblResults && multipleChemblResults.length > 0 && searchType === "multiple_chembl_ids" && (
                    <div className="chembl-getter-row-2" style={{ marginTop: '20px' }}>
                        <GlassyContainer>
                            <h3 style={{ marginBottom: '15px', fontWeight: '700' }}>
                                ChEMBL Molecules Results {multipleChemblResults.length > 0 ? `(${multipleChemblResults.length} found)` : ''}
                            </h3>
                            
                            <div className="chembl-results-container">
                                {/* Left Panel - Molecule List */}
                                <div className="chembl-candidates-panel">
                                    <h4 style={{
                                        color: 'var(--color-text-primary)',
                                        fontSize: '1rem',
                                        marginBottom: '10px',
                                        fontWeight: '500'
                                    }}>
                                        Found Molecules
                                    </h4>
                                    <p style={{
                                        color: 'var(--color-text-secondary)',
                                        fontSize: '0.8rem',
                                        marginBottom: '12px',
                                        fontStyle: 'italic'
                                    }}>
                                        💡 Use ↑↓ arrow keys to navigate
                                    </p>
                                    <div className="chembl-candidates-list">
                                        {multipleChemblResults.map((molecule, index) => (
                                            <div 
                                                key={index}
                                                onClick={() => {
                                                    setSelectedChemblResult(molecule);
                                                    setSelectedChemblIndex(index);
                                                }}
                                                className={`chembl-candidate-item ${selectedChemblResult === molecule ? 'selected' : ''}`}
                                            >
                                                <span className="chembl-candidate-number">
                                                    {index + 1}.
                                                </span>
                                                <span className="chembl-candidate-smiles">
                                                    {molecule.molecule_chembl_id || molecule.pref_name || 'Unknown'}
                                                </span>
                                            </div>
                                        ))}
                                    </div>
                                </div>

                                {/* Right Panel - Info and 2D Viewer */}
                                <div className="chembl-details-panel">
                                    {selectedChemblResult ? (
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
                                                    Selected ({selectedChemblIndex + 1}/{multipleChemblResults.length}): <span style={{ 
                                                        fontFamily: 'monospace', 
                                                        color: 'var(--color-accent)',
                                                        fontWeight: '600' 
                                                    }}>
                                                        {selectedChemblResult.molecule_chembl_id || selectedChemblResult.pref_name || 'Unknown'}
                                                    </span>
                                                </h4>
                                            </div>
                                            
                                            {/* Show visualization only if we have SMILES data */}
                                            {selectedChemblResult.molecule_structures && selectedChemblResult.molecule_structures.canonical_smiles ? (
                                                <div className="chembl-details-row">
                                                    {/* InfoBox */}
                                                    <div style={{ flex: '4' }}>
                                                        <InfoBox 
                                                            activeMol={selectedChemblResult.molecule_structures.canonical_smiles} 
                                                            isValidMol={true} 
                                                            infoType={"MOL"} 
                                                        />
                                                    </div>
                                                    
                                                    {/* TwoDViewer */}
                                                    <div style={{ flex: '6' }}>
                                                        <TwoDViewer 
                                                            activeMol={selectedChemblResult.molecule_structures.canonical_smiles} 
                                                            isValidMol={true} 
                                                            visType={"MOL"} 
                                                        />
                                                    </div>
                                                </div>
                                            ) : (
                                                <div className="chembl-placeholder">
                                                    No structure data available for this molecule
                                                </div>
                                            )}
                                        </>
                                    ) : (
                                        <div className="chembl-placeholder">
                                            Select a molecule from the left panel to view details
                                        </div>
                                    )}
                                </div>
                            </div>
                        </GlassyContainer>
                    </div>
                )}

                {/* DataViewer - show single molecule data or selected candidate data */}
                {(apiData || selectedCandidate || selectedChemblResult) && (
                    <div className="chembl-getter-row-3" style={{ marginTop: '20px' }}>
                        <DataViewer
                            data={
                                searchType === "chembl_id" ? apiData : 
                                searchType === "multiple_chembl_ids" ? selectedChemblResult :
                                selectedCandidate
                            }
                            title={`ChEMBL Molecule Data${searchValue ? ` for "${searchValue}"` : ''}${
                                searchType === "pref_name" && selectedCandidate ? 
                                ` - ${selectedCandidate.pref_name || selectedCandidate.molecule_chembl_id || 'Unknown'}` : 
                                searchType === "multiple_chembl_ids" && selectedChemblResult ?
                                ` - ${selectedChemblResult.molecule_chembl_id || selectedChemblResult.pref_name || 'Unknown'}` :
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

// ChemBLGetter.propTypes = {
//     toolData: PropTypes.object,
//     initialSearchValue: PropTypes.string,
//     initialSearchType: PropTypes.oneOf(["chembl_id", "multiple_chembl_ids", "pref_name", "synonym"]),
//     hideInputBox: PropTypes.bool,
// };

export default ChemBLGetter;
