import React, { useState, useEffect, useCallback, useRef } from "react";
import { createPortal } from "react-dom";
import "./PubChemGetter.css";
import SimpleInputBox from "../../../components/UI/SimpleInputBox/SimpleInputBox";
import DataViewer from "../../../components/UI/DataViewer/DataViewer";
import InfoBox from "../../../components/UI/InfoBox/InfoBox";
import TwoDViewer from "../../../components/UI/TwoDViewer/TwoDViewer";
import { motion } from "framer-motion";
import { fadeInUpVariantStatic } from "../../../components/animations/framerAnim";
import GlassyContainer from "../../../components/glassy_container/gc";
import { 
    getCompoundByCid, 
    getCidsByName, 
    getCidsBySmiles, 
    getCidsByInchikey,
    getCompoundProperties 
} from "../../../services/api/mcpToolsService";

// PubChem CID Validator
const validateCid = (cid) => {
    // CID must be a positive integer
    const cidPattern = /^\d+$/;
    return cidPattern.test(cid.trim()) && parseInt(cid.trim()) > 0;
};

// In-house PubChem Search Type Selector Component with Portal
const PubChemSearchTypeSelector = ({ value, onChange, options, disabled = false }) => {
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
                className={`pubchem-search-selector ${disabled ? 'disabled' : ''}`}
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
const searchTypeOptions = [
    { value: "cid", label: "PubChem CID" },
    { value: "name", label: "Compound Name" },
    { value: "smiles", label: "SMILES String" },
    { value: "inchikey", label: "InChIKey" }
];

const PubChemGetter = ({
    toolData = null,
    initialSearchValue = "",
    initialSearchType = "cid",
    hideInputBox = false
}) => {
    const [searchValue, setSearchValue] = useState(initialSearchValue);
    const [inputValue, setInputValue] = useState(initialSearchValue);
    const [searchType, setSearchType] = useState(initialSearchType);
    const [apiData, setApiData] = useState(toolData);
    const [isLoading, setIsLoading] = useState(false);
    const [error, setError] = useState(null);
    const [validationError, setValidationError] = useState(null);
    
    // States for multiple candidates (name/SMILES/InChIKey search that might return multiple CIDs)
    const [multipleCandidates, setMultipleCandidates] = useState([]);
    const [selectedCandidate, setSelectedCandidate] = useState(null);
    const [selectedIndex, setSelectedIndex] = useState(0);
    const [isSelectedCandidateValid, setIsSelectedCandidateValid] = useState(true);

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

        // Only validate if there's a value and search type is CID
        if (value && value.trim() && searchType === "cid") {
            if (!validateCid(value.trim())) {
                setValidationError("Invalid CID format. Must be a positive integer.");
            }
        }
    }, [searchType]);

    const handleSubmit = useCallback(async (value) => {
        // Original logic for search types
        if (!value || value.trim() === "") {
            setError("Please enter a search value.");
            return;
        }

        const trimmedValue = value.trim();

        // Validate CID format if search type is CID
        if (searchType === "cid") {
            if (!validateCid(trimmedValue)) {
                setError("Invalid CID format. CID must be a positive integer (e.g., 2244, 5988).");
                return;
            }
        }

        setIsLoading(true);
        setError(null);
        setValidationError(null);

        try {
            let result;

            if (searchType === "cid") {
                result = await getCompoundByCid({ cid: trimmedValue });
            } else if (searchType === "name") {
                result = await getCidsByName({ name: trimmedValue });
            } else if (searchType === "smiles") {
                result = await getCidsBySmiles({ smiles: trimmedValue });
            } else if (searchType === "inchikey") {
                result = await getCidsByInchikey({ inchikey: trimmedValue });
            }

            if (result) {
                if (result.status === "success") {
                    var processible = result.result["0"]
                    processible = JSON.parse(processible);
                    
                    console.log('Processible', processible);

                    if (searchType === "cid") {
                        // For CID search, we expect a single compound
                        if (processible.data){
                            setApiData(processible.data);
                            setSearchValue(trimmedValue);
                            // Clear multiple candidates for single compound view
                            setMultipleCandidates([]);
                            setSelectedCandidate(null);
                        } else {
                            // No data found for CID
                            setError(`No compound found for CID "${trimmedValue}"`);
                            setApiData(null);
                            setMultipleCandidates([]);
                            setSelectedCandidate(null);
                        }
                    } else {
                        // For name/SMILES/InChIKey search, we might get multiple CIDs
                        if (processible.data && Array.isArray(processible.data) && processible.data.length > 0) {
                            // Multiple CIDs found - we need to fetch compound data for each
                            const cids = processible.data.slice(0, 10); // Limit to first 10 results
                            const compounds = [];
                            
                            for (const cid of cids) {
                                try {
                                    const compoundResult = await getCompoundByCid({ cid: cid.toString() });
                                    if (compoundResult && compoundResult.status === "success") {
                                        const compoundData = JSON.parse(compoundResult.result["0"]);
                                        if (compoundData.data) {
                                            compounds.push({
                                                ...compoundData.data,
                                                cid: cid
                                            });
                                        }
                                    }
                                } catch (err) {
                                    console.error(`Error fetching compound for CID ${cid}:`, err);
                                }
                            }
                            
                            if (compounds.length > 0) {
                                setMultipleCandidates(compounds);
                                setSelectedCandidate(compounds[0]);
                                setSelectedIndex(0);
                                setIsSelectedCandidateValid(true);
                                setApiData(null); // Clear single compound data
                                setSearchValue(trimmedValue);
                            } else {
                                setError(`No compound data found for search "${trimmedValue}"`);
                                setMultipleCandidates([]);
                                setSelectedCandidate(null);
                                setApiData(null);
                            }
                        } else if (processible.data && Array.isArray(processible.data) && processible.data.length === 0) {
                            // Empty array - no records found
                            console.log(`No records found for ${searchType} search "${trimmedValue}"`);
                            setError(`No records found with ${searchType} search "${trimmedValue}"`);
                            setMultipleCandidates([]);
                            setSelectedCandidate(null);
                            setApiData(null);
                        } else if (processible.data) {
                            // Single CID found
                            try {
                                const compoundResult = await getCompoundByCid({ cid: processible.data.toString() });
                                if (compoundResult && compoundResult.status === "success") {
                                    const compoundData = JSON.parse(compoundResult.result["0"]);
                                    if (compoundData.data) {
                                        setApiData({
                                            ...compoundData.data,
                                            cid: processible.data
                                        });
                                        setSearchValue(trimmedValue);
                                        setMultipleCandidates([]);
                                        setSelectedCandidate(null);
                                    } else {
                                        setError(`No compound data found for search "${trimmedValue}"`);
                                        setApiData(null);
                                        setMultipleCandidates([]);
                                        setSelectedCandidate(null);
                                    }
                                } else {
                                    setError(`Failed to fetch compound data for search "${trimmedValue}"`);
                                    setApiData(null);
                                    setMultipleCandidates([]);
                                    setSelectedCandidate(null);
                                }
                            } catch (err) {
                                console.error(`Error fetching compound:`, err);
                                setError(`Error fetching compound data for search "${trimmedValue}"`);
                                setApiData(null);
                                setMultipleCandidates([]);
                                setSelectedCandidate(null);
                            }
                        } else {
                            // No data found
                            setError(`No records found with ${searchType} search "${trimmedValue}"`);
                            setMultipleCandidates([]);
                            setSelectedCandidate(null);
                            setApiData(null);
                        }
                    }
                } else {
                    // API returned non-success status
                    setError(result.message || `Search failed for "${trimmedValue}". Please try again.`);
                    setApiData(null);
                    setMultipleCandidates([]);
                    setSelectedCandidate(null);
                }
            } else {
                setError("No data found for the provided search value.");
                setApiData(null);
                setMultipleCandidates([]);
                setSelectedCandidate(null);
            }
        } catch (error) {
            console.error('Error fetching PubChem data:', error);
            setError(error.message || "Failed to fetch PubChem data. Please try again.");
        } finally {
            setIsLoading(false);
        }
    }, [searchType]);

    // Effect to handle initial data passed as props
    useEffect(() => {
        if (toolData && toolData.type && toolData.content) {
            console.log('Processing tool data:', toolData);
            
            // Set search type based on tool data type
            setSearchType(toolData.type);
            
            if (toolData.type === "name" || toolData.type === "smiles" || toolData.type === "inchikey") {
                // Handle search results that might return multiple compounds
                if (Array.isArray(toolData.content) && toolData.content.length > 0) {
                    setMultipleCandidates(toolData.content);
                    setSelectedCandidate(toolData.content[0]);
                    setSelectedIndex(0);
                    setIsSelectedCandidateValid(true);
                    setApiData(null); // Clear single compound data
                } else if (toolData.content) {
                    // If it's a single candidate, treat it as a single-item array for consistent UI
                    setMultipleCandidates([toolData.content]);
                    setSelectedCandidate(toolData.content);
                    setSelectedIndex(0);
                    setIsSelectedCandidateValid(true);
                    setApiData(null); // Clear single compound data
                } else {
                    // No candidates found
                    setError(`No records found with ${toolData.type} search.`);
                }
            } else if (toolData.type === "cid") {
                // Handle single CID search result
                if (toolData.content) {
                    setApiData(toolData.content);
                    setMultipleCandidates([]);
                    setSelectedCandidate(null);
                } else {
                    setError("No compound found for the provided CID.");
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
            case "cid":
                return "Format: positive integer (e.g., 2244, 5988)";
            case "name":
                return "Try name: Aspirin";
            case "smiles":
                return "Try SMILES: CC(=O)OC1=CC=CC=C1C(=O)O";
            case "inchikey":
                return "Try InChIKey: BSYNRYMUTXBXSQ-UHFFFAOYSA-N";
            default:
                return "Enter search value";
        }
    };

    const getHeaderText = () => {
        switch (searchType) {
            case "cid":
                return "Enter PubChem CID";
            case "name":
                return "Enter Compound Name";
            case "smiles":
                return "Enter SMILES String";
            case "inchikey":
                return "Enter InChIKey";
            default:
                return "Enter Search Value";
        }
    };

    useEffect(() => {
        console.log("Search type changed:", searchType);
        // Clear validation errors when search type changes, but keep existing errors
        setValidationError(null);
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
        };

        window.addEventListener('keydown', handleKeyDown);
        return () => window.removeEventListener('keydown', handleKeyDown);
    }, [multipleCandidates]);

    useEffect(() => {
        console.log("Error state changed:", error);    
    }, [error])

    return (
        <>
            <motion.div
                initial="hidden"
                animate="visible"
                variants={fadeInUpVariantStatic}
                className="pubchem-getter-container"
            >
                {!hideInputBox && (
                    <div className="pubchem-getter-row-1">
                        <GlassyContainer>
                            <div className="pubchem-getter-search-controls">
                                <div className="pubchem-getter-dropdown-section">
                                    <h3 className="pubchem-getter-dropdown-header">Search Type</h3>
                                    <PubChemSearchTypeSelector
                                        value={searchType}
                                        onChange={(e) => setSearchType(e.target.value)}
                                        options={searchTypeOptions}
                                        disabled={isLoading}
                                    />
                                </div>

                                <div className="pubchem-getter-input-section">
                                    <SimpleInputBox
                                        value={inputValue}
                                        onChange={handleInputChange}
                                        onSubmit={handleSubmit}
                                        header={getHeaderText()}
                                        placeholder={getPlaceholderText()}
                                        buttonText="Search"
                                        isLoading={isLoading}
                                        disabled={isLoading}
                                    />
                                </div>
                            </div>

                            {/* Display validation or general errors */}
                            {(validationError || error) && (
                                <div className="pubchem-getter-error" style={{ 
                                    color: 'var(--color-alert)', 
                                    marginTop: '8px', 
                                    fontSize: '0.9rem' 
                                }}>
                                    {validationError || error}
                                </div>
                            )}
                        </GlassyContainer>
                    </div>
                )}

                {/* Single compound visualization for CID search */}
                {apiData && apiData.atoms && searchType === "cid" && (
                    <div className="pubchem-getter-row-2" style={{ marginTop: '20px' }}>
                        <GlassyContainer>
                            <h3 style={{ marginBottom: '15px', fontWeight: '700' }}>Compound Structure Visualization</h3>
                            <div style={{ display: 'flex', gap: '20px' }}>
                                <div style={{ flex: '1' }}>
                                    <InfoBox 
                                        activeMol={apiData.canonical_smiles || ""} 
                                        isValidMol={true} 
                                        infoType={"MOL"} 
                                    />
                                </div>
                                <div style={{ flex: '1' }}>
                                    <TwoDViewer 
                                        activeMol={apiData.canonical_smiles || ""} 
                                        isValidMol={true} 
                                        visType={"MOL"} 
                                    />
                                </div>
                            </div>
                        </GlassyContainer>
                    </div>
                )}

                {/* Multiple candidates view for name/SMILES/InChIKey search */}
                {multipleCandidates && multipleCandidates.length > 0 && (searchType === "name" || searchType === "smiles" || searchType === "inchikey") && (
                    <div className="pubchem-getter-row-2" style={{ marginTop: '20px' }}>
                        <GlassyContainer>
                            <h3 style={{ marginBottom: '15px', fontWeight: '700' }}>
                                Compound Candidates {multipleCandidates.length > 0 ? `(${multipleCandidates.length} found)` : ''}
                            </h3>
                            
                            <div className="pubchem-results-container">
                                {/* Left Panel - Candidate List */}
                                <div className="pubchem-candidates-panel">
                                    <h4 style={{
                                        color: 'var(--color-text-primary)',
                                        fontSize: '1rem',
                                        marginBottom: '10px',
                                        fontWeight: '500'
                                    }}>
                                        Candidate Compounds
                                    </h4>
                                    <p style={{
                                        color: 'var(--color-text-secondary)',
                                        fontSize: '0.8rem',
                                        marginBottom: '12px',
                                        fontStyle: 'italic'
                                    }}>
                                        💡 Use ↑↓ arrow keys to navigate
                                    </p>
                                    <div className="pubchem-candidates-list">
                                        {multipleCandidates.map((candidate, index) => (
                                            <div 
                                                key={index}
                                                onClick={() => {
                                                    setSelectedCandidate(candidate);
                                                    setIsSelectedCandidateValid(true);
                                                    setSelectedIndex(index);
                                                }}
                                                className={`pubchem-candidate-item ${selectedCandidate === candidate ? 'selected' : ''}`}
                                            >
                                                <span className="pubchem-candidate-number">
                                                    {index + 1}.
                                                </span>
                                                <span className="pubchem-candidate-smiles">
                                                    {candidate.iupac_name || candidate.title || `CID: ${candidate.cid}` || 'Unknown'}
                                                </span>
                                            </div>
                                        ))}
                                    </div>
                                </div>

                                {/* Right Panel - Info and 2D Viewer */}
                                <div className="pubchem-details-panel">
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
                                                        {selectedCandidate.iupac_name || selectedCandidate.title || `CID: ${selectedCandidate.cid}` || 'Unknown'}
                                                    </span>
                                                </h4>
                                            </div>
                                            
                                            {/* Show visualization only if we have SMILES data */}
                                            {selectedCandidate.canonical_smiles ? (
                                                <div className="pubchem-details-row">
                                                    {/* InfoBox */}
                                                    <div style={{ flex: '4' }}>
                                                        <InfoBox 
                                                            activeMol={selectedCandidate.canonical_smiles} 
                                                            isValidMol={isSelectedCandidateValid} 
                                                            infoType={"MOL"} 
                                                        />
                                                    </div>
                                                    
                                                    {/* TwoDViewer */}
                                                    <div style={{ flex: '6' }}>
                                                        <TwoDViewer 
                                                            activeMol={selectedCandidate.canonical_smiles} 
                                                            isValidMol={isSelectedCandidateValid} 
                                                            visType={"MOL"} 
                                                        />
                                                    </div>
                                                </div>
                                            ) : (
                                                <div className="pubchem-placeholder">
                                                    No structure data available for this candidate
                                                </div>
                                            )}
                                        </>
                                    ) : (
                                        <div className="pubchem-placeholder">
                                            Select a candidate compound from the left panel to view details
                                        </div>
                                    )}
                                </div>
                            </div>
                        </GlassyContainer>
                    </div>
                )}

                {/* DataViewer - show single compound data or selected candidate data */}
                {(apiData || selectedCandidate) && (
                    <div className="pubchem-getter-row-3" style={{ marginTop: '20px' }}>
                        <DataViewer
                            data={
                                searchType === "cid" ? apiData : 
                                selectedCandidate
                            }
                            title={`PubChem Compound Data${searchValue ? ` for "${searchValue}"` : ''}${
                                searchType !== "cid" && selectedCandidate ? 
                                ` - ${selectedCandidate.iupac_name || selectedCandidate.title || `CID: ${selectedCandidate.cid}` || 'Unknown'}` : 
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

export default PubChemGetter;
