import React, { useState, useEffect, useCallback, useRef } from 'react';
import { createPortal } from 'react-dom';
import PropTypes from 'prop-types';
import {
    rcsbStructureSimilarityByEntryId,
    rcsbStructureSimilarityByFileUrl
} from '../../../services/api/mcpToolsService'; // Adjust path as needed
import DataViewer from '../../../components/UI/DataViewer/DataViewer'; // Adjusted path
import './RCSBPDBExplorer.css'; // Import the same CSS file for consistent styling

// Debounce function
const debounce = (func, delay) => {
    let timeout;
    return function(...args) {
        const context = this;
        clearTimeout(timeout);
        timeout = setTimeout(() => func.apply(context, args), delay);
    };
};

// Moved PortalDropdown to be a standalone component
const StandalonePortalDropdownMenu = ({
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
            className="activity-search-dropdown-portal"
            style={{
                position: 'fixed',
                top: `${dropdownPosition.top}px`,
                left: `${dropdownPosition.left}px`,
                width: 'auto',
                minWidth: `${dropdownPosition.width}px`,
                maxWidth: '400px',
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
                    className={`activity-selector-option ${option.value === currentValue ? 'selected' : ''}`}
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
                            e.currentTarget.style.backgroundColor = 'var(--color-bg-secondary)';
                        }
                    }}
                    onMouseLeave={(e) => {
                        if (option.value !== currentValue) {
                            e.currentTarget.style.backgroundColor = 'transparent';
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
                                color: option.value === currentValue ? 'rgba(255,255,255,0.8)' : 'var(--color-text-secondary)',
                                opacity: 0.8,
                                lineHeight: '1.3',
                                wordWrap: 'break-word'
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

// Activity Search Type Selector Component with Portal
const ActivitySearchTypeSelector = ({ value, onChange, options = [], disabled = false, placeholder = "Select an option" }) => {
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

    const selectedOption = options?.find(option => option.value === value);

    return (
        <>
            <div className="dropdown-field-container">
                <label className="dropdown-field-label">{placeholder}</label>
                <div
                    className={`activity-search-selector ${disabled ? 'disabled' : ''}`}
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
                        style={{
                            display: 'flex',
                            alignItems: 'center',
                            justifyContent: 'space-between',
                            padding: '10px 15px',
                            backgroundColor: 'var(--glassy-color)',
                            border: `1px solid ${isOpen ? 'var(--color-accent)' : 'var(--c-light-border)'}`,
                            borderRadius: isOpen ? '15px 15px 0 0' : '15px',
                            borderBottom: isOpen ? 'none' : `1px solid var(--c-light-border)`,
                            cursor: disabled ? 'not-allowed' : 'pointer',
                            transition: 'all 0.2s ease',
                            minHeight: '50px',
                            minWidth: '200px',
                            maxWidth: '300px',
                            width: 'auto',
                            opacity: disabled ? '0.6' : '1',
                            color: 'var(--color-text-primary)',
                            fontSize: '16px',
                            fontWeight: '600'
                        }}
                        onMouseEnter={(e) => {
                            if (!disabled && !isOpen) {
                                e.target.style.borderColor = 'var(--color-accent)';
                                e.target.style.backgroundColor = 'var(--color-bg-secondary)';
                            }
                        }}
                        onMouseLeave={(e) => {
                            if (!disabled && !isOpen) {
                                e.target.style.borderColor = 'var(--c-light-border)';
                                e.target.style.backgroundColor = 'var(--glassy-color)';
                            }
                        }}
                    >
                        <span className="selected-text" style={{
                            color: 'var(--color-text-primary)',
                            fontSize: '16px',
                            fontWeight: '600'
                        }}>
                            {selectedOption ? selectedOption.label : 'Select option'}
                        </span>
                        <span className={`dropdown-arrow ${isOpen ? 'open' : ''}`} style={{
                            color: 'var(--color-text-secondary)',
                            display: 'flex',
                            alignItems: 'center',
                            justifyContent: 'center',
                            transition: 'transform 0.2s ease',
                            transform: isOpen ? 'rotate(180deg)' : 'rotate(0deg)'
                        }}>
                            <svg
                                width="12"
                                height="8"
                                viewBox="0 0 12 8"
                                fill="currentColor"
                            >
                                <path d="M6 8L0 0h12z" />
                            </svg>
                        </span>
                    </div>
                </div>
            </div>
            <StandalonePortalDropdownMenu
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

// Generic Dropdown Component using the exact same structure as ActivitySearchTypeSelector
const GenericDropdownSelector = ({ value, onChange, options = [], disabled = false, placeholder = "Select an option" }) => {
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

    const selectedOption = options?.find(option => option.value === value);

    return (
        <>
            <div className="dropdown-field-container">
                <label className="dropdown-field-label">{placeholder}</label>
                <div
                    className={`activity-search-selector ${disabled ? 'disabled' : ''}`}
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
                        style={{
                            display: 'flex',
                            alignItems: 'center',
                            justifyContent: 'space-between',
                            padding: '10px 15px',
                            backgroundColor: 'var(--glassy-color)',
                            border: `1px solid ${isOpen ? 'var(--color-accent)' : 'var(--c-light-border)'}`,
                            borderRadius: isOpen ? '15px 15px 0 0' : '15px',
                            borderBottom: isOpen ? 'none' : `1px solid var(--c-light-border)`,
                            cursor: disabled ? 'not-allowed' : 'pointer',
                            transition: 'all 0.2s ease',
                            minHeight: '50px',
                            minWidth: '200px',
                            maxWidth: '300px',
                            width: 'auto',
                            opacity: disabled ? '0.6' : '1',
                            color: 'var(--color-text-primary)',
                            fontSize: '16px',
                            fontWeight: '600'
                        }}
                        onMouseEnter={(e) => {
                            if (!disabled && !isOpen) {
                                e.target.style.borderColor = 'var(--color-accent)';
                                e.target.style.backgroundColor = 'var(--color-bg-secondary)';
                            }
                        }}
                        onMouseLeave={(e) => {
                            if (!disabled && !isOpen) {
                                e.target.style.borderColor = 'var(--c-light-border)';
                                e.target.style.backgroundColor = 'var(--glassy-color)';
                            }
                        }}
                    >
                        <span className="selected-text" style={{
                            color: 'var(--color-text-primary)',
                            fontSize: '16px',
                            fontWeight: '600'
                        }}>
                            {selectedOption ? selectedOption.label : 'Select option'}
                        </span>
                        <span className={`dropdown-arrow ${isOpen ? 'open' : ''}`} style={{
                            color: 'var(--color-text-secondary)',
                            display: 'flex',
                            alignItems: 'center',
                            justifyContent: 'center',
                            transition: 'transform 0.2s ease',
                            transform: isOpen ? 'rotate(180deg)' : 'rotate(0deg)'
                        }}>
                            <svg
                                width="12"
                                height="8"
                                viewBox="0 0 12 8"
                                fill="currentColor"
                            >
                                <path d="M6 8L0 0h12z" />
                            </svg>
                        </span>
                    </div>
                </div>
            </div>
            <StandalonePortalDropdownMenu
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


const RCSBStructureSimilaritySearch = ({ initialSearchType = "entry_id_search", toolData = null }) => {
    const [searchType, setSearchType] = useState(initialSearchType);
    const [apiData, setApiData] = useState(toolData);
    const [isLoading, setIsLoading] = useState(false);
    const [error, setError] = useState(null);

    // Input states for each search type
    const [entryId, setEntryId] = useState("");
    const [assemblyId, setAssemblyId] = useState("1");
    const [operator, setOperator] = useState("strict_shape_match");
    const [targetSearchSpace, setTargetSearchSpace] = useState("assembly");
    const [returnType, setReturnType] = useState("assembly");
    const [fileUrl, setFileUrl] = useState("");
    const [fileFormat, setFileFormat] = useState("cif");

    // Options for dropdowns
    const searchTypeOptions = [
        { 
            value: "entry_id_search", 
            label: "PDB Entry ID Search",
            description: "Find structures similar to a given PDB entry ID"
        },
        { 
            value: "file_url_search", 
            label: "File URL Search",
            description: "Find structures similar to a structure file provided via URL"
        },
    ];

    const operatorOptions = [
        { 
            value: "strict_shape_match", 
            label: "Strict Shape Match",
            description: "More stringent similarity matching"
        },
        { 
            value: "relaxed_shape_match", 
            label: "Relaxed Shape Match",
            description: "More lenient similarity matching"
        },
    ];
    
    const targetSearchSpaceOptions = [
        { 
            value: "assembly", 
            label: "Assembly",
            description: "Compare against assemblies"
        },
        { 
            value: "polymer_entity_instance", 
            label: "Polymer Entity Instance",
            description: "Compare against polymer entity instances"
        },
    ];

    const returnTypeOptions = [
        { 
            value: "assembly", 
            label: "Assembly",
            description: "Return assembly identifiers"
        },
        { 
            value: "polymer_entity", 
            label: "Polymer Entity",
            description: "Return polymer entity identifiers"
        },
    ];

    const fileFormatOptions = [
        { value: "cif", label: "CIF" },
        { value: "bcif", label: "BCIF" },
        { value: "pdb", label: "PDB" },
        { value: "cif.gz", label: "CIF (Gzipped)" },
        { value: "pdb.gz", label: "PDB (Gzipped)" },
    ];


    const handleApiCall = async (args) => {
        setIsLoading(true);
        setError(null);
        setApiData(null);
        try {
            let result;
            switch (searchType) {
                case "entry_id_search":
                    result = await rcsbStructureSimilarityByEntryId({ 
                        entry_id: entryId, 
                        assembly_id: assemblyId,
                        operator,
                        target_search_space: targetSearchSpace,
                        return_type: returnType,
                        ...args 
                    });
                    break;
                case "file_url_search":
                    result = await rcsbStructureSimilarityByFileUrl({ 
                        file_url: fileUrl,
                        file_format: fileFormat,
                        operator,
                        target_search_space: targetSearchSpace,
                        return_type: returnType,
                        ...args 
                    });
                    break;
                default:
                    throw new Error("Invalid search type");
            }
            setApiData(result);
        } catch (err) {
            setError(err.message || "An error occurred during the structure similarity search.");
            console.error("RCSB PDB Structure Similarity Search Error:", err);
        } finally {
            setIsLoading(false);
        }
    };
    
    const debouncedApiCall = useCallback(debounce(handleApiCall, 500), [
        searchType, entryId, assemblyId, operator, targetSearchSpace, returnType, 
        fileUrl, fileFormat
    ]);

    useEffect(() => {
        if (toolData) {
            setApiData(toolData);
        }
    }, [toolData]);

    const handleSubmit = (e) => {
        e.preventDefault();
        // Validate inputs before calling API
        let valid = true;
        switch (searchType) {
            case "entry_id_search":
                if (!entryId.trim()) {
                    setError("PDB Entry ID cannot be empty.");
                    valid = false;
                }
                break;
            case "file_url_search":
                if (!fileUrl.trim()) {
                    setError("File URL cannot be empty.");
                    valid = false;
                }
                if (!fileFormat.trim()) {
                    setError("File format must be selected.");
                    valid = false;
                }
                break;
            default:
                break;
        }
        if (valid) {
            setError(null); // Clear previous errors
            handleApiCall({}); // Pass empty object if no extra args needed immediately
        }
    };

    const renderSearchInputs = () => {
        switch (searchType) {
            case "entry_id_search":
                return (
                    <div className="search-inputs-container">
                        <div className="input-field-container">
                            <label className="input-field-label">PDB Entry ID:</label>
                            <input
                                type="text"
                                value={entryId}
                                onChange={(e) => setEntryId(e.target.value)}
                                placeholder="Enter PDB ID (e.g., 4HHB)"
                                className="rcsb-input-field"
                            />
                        </div>
                        <div className="input-field-container">
                            <label className="input-field-label">Assembly ID:</label>
                            <input
                                type="text"
                                value={assemblyId}
                                onChange={(e) => setAssemblyId(e.target.value)}
                                placeholder="Assembly ID (default: 1)"
                                className="rcsb-input-field"
                            />
                        </div>
                        <GenericDropdownSelector
                            value={operator}
                            onChange={setOperator}
                            options={operatorOptions}
                            placeholder="Select Similarity Operator"
                        />
                        <GenericDropdownSelector
                            value={targetSearchSpace}
                            onChange={setTargetSearchSpace}
                            options={targetSearchSpaceOptions}
                            placeholder="Select Target Search Space"
                        />
                        <GenericDropdownSelector
                            value={returnType}
                            onChange={setReturnType}
                            options={returnTypeOptions}
                            placeholder="Select Return Type"
                        />
                    </div>
                );
            case "file_url_search":
                return (
                    <div className="search-inputs-container">
                        <div className="input-field-container">
                            <label className="input-field-label">File URL:</label>
                            <input
                                type="text"
                                value={fileUrl}
                                onChange={(e) => setFileUrl(e.target.value)}
                                placeholder="URL to structure file (e.g., https://files.rcsb.org/view/4HHB.cif)"
                                className="rcsb-input-field"
                            />
                        </div>
                        <GenericDropdownSelector
                            value={fileFormat}
                            onChange={setFileFormat}
                            options={fileFormatOptions}
                            placeholder="Select File Format"
                        />
                        <GenericDropdownSelector
                            value={operator}
                            onChange={setOperator}
                            options={operatorOptions}
                            placeholder="Select Similarity Operator"
                        />
                        <GenericDropdownSelector
                            value={targetSearchSpace}
                            onChange={setTargetSearchSpace}
                            options={targetSearchSpaceOptions}
                            placeholder="Select Target Search Space"
                        />
                        <GenericDropdownSelector
                            value={returnType}
                            onChange={setReturnType}
                            options={returnTypeOptions}
                            placeholder="Select Return Type"
                        />
                    </div>
                );
            default:
                return <p>Select a search type.</p>;
        }
    };

    return (
        <div className="rcsb-pdb-explorer-container">
            <div className="rcsb-search-type-selector-section">
                <ActivitySearchTypeSelector
                    value={searchType}
                    onChange={setSearchType}
                    options={searchTypeOptions}
                    placeholder="Select Search Type"
                />
            </div>

            <form onSubmit={handleSubmit} className="rcsb-search-form">
                {renderSearchInputs()}
                <button type="submit" disabled={isLoading} className="rcsb-submit-button">
                    {isLoading ? 'Searching...' : 'Find Similar Structures'}
                </button>
            </form>

            {error && <p style={{ color: 'red', marginTop: '10px' }}>Error: {error}</p>}

            <DataViewer data={apiData} />
        </div>
    );
};

RCSBStructureSimilaritySearch.propTypes = {
    initialSearchType: PropTypes.oneOf([
        "entry_id_search",
        "file_url_search",
    ]),
    toolData: PropTypes.object,
};

export default RCSBStructureSimilaritySearch;

// Basic CSS for PortalDropdown (can be moved to a separate CSS file)
const styleSheet = document.createElement("style");
styleSheet.type = "text/css";
styleSheet.innerText = `
  .portal-dropdown .dropdown-option:hover {
    background-color: #f0f0f0;
  }
  .portal-dropdown .dropdown-option.selected {
    background-color: #e0e0e0;
    font-weight: bold;
  }
`;
document.head.appendChild(styleSheet);
