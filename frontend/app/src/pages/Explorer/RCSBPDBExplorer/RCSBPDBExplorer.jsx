import React, { useState, useEffect, useCallback, useRef } from 'react';
import { createPortal } from 'react-dom';
import PropTypes from 'prop-types';
import {
    rcsbTextSearchPdb,
    rcsbAttributeSearchPdb,
    rcsbCombinedTextAndAttributeSearch,
    rcsbSequenceIdentitySearch,
    rcsbSequenceMotifSearch,
    rcsbGetProteinDetailsByIdPypdb
} from '../../../services/api/mcpToolsService'; // Adjust path as needed
import DataViewer from '../../../components/UI/DataViewer/DataViewer'; // Adjusted path
import './RCSBPDBExplorer.css'; // Import the new CSS file

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
            {/* Use the standalone portal dropdown component */}
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
            {/* Use the standalone portal dropdown component */}
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


const RCSBPDBExplorer = ({ initialSearchType = "text_search", toolData = null }) => {
    const [searchType, setSearchType] = useState(initialSearchType);
    const [apiData, setApiData] = useState(toolData);
    const [isLoading, setIsLoading] = useState(false);
    const [error, setError] = useState(null);

    // Input states for each search type
    const [textQuery, setTextQuery] = useState("");
    const [pdbId, setPdbId] = useState(""); // New state for PDB ID search
    const [attributePath, setAttributePath] = useState("");
    const [operator, setOperator] = useState("exact_match"); // Default operator
    const [value, setValue] = useState("");
    const [combinedTextQuery, setCombinedTextQuery] = useState("");
    const [attributeFilters, setAttributeFilters] = useState([{ attribute_path: "", operator: "exact_match", value: "" }]);
    const [sequence, setSequence] = useState("");
    const [identityCutoff, setIdentityCutoff] = useState(0.9);
    const [eValueCutoff, setEValueCutoff] = useState(1.0);
    const [sequenceType, setSequenceType] = useState("protein");
    const [motifPattern, setMotifPattern] = useState("");
    const [patternType, setPatternType] = useState("prosite");

    // Options for dropdowns
    const searchTypeOptions = [
        { value: "text_search", label: "Text Search" },
        { value: "pdb_id_search", label: "PDB ID Search" },
        { value: "attribute_search", label: "Attribute Search" },
        { value: "combined_search", label: "Combined Text & Attribute Search" },
        { value: "sequence_identity_search", label: "Sequence Identity Search" },
        { value: "sequence_motif_search", label: "Sequence Motif Search" },
    ];

    const operatorOptions = [
        { value: "exact_match", label: "Exact Match" },
        { value: "greater", label: "Greater Than" },
        { value: "greater_or_equal", label: "Greater Than or Equal To" },
        { value: "less", label: "Less Than" },
        { value: "less_or_equal", label: "Less Than or Equal To" },
        { value: "in", label: "In (comma-separated list)" },
        { value: "contains_phrase", label: "Contains Phrase" },
        { value: "contains_words", label: "Contains Words" },
    ];
    
    const sequenceTypeOptions = [
        { value: "protein", label: "Protein" },
        { value: "dna", label: "DNA" },
        { value: "rna", label: "RNA" },
    ];

    const patternTypeOptions = [
        { value: "prosite", label: "PROSITE" },
        { value: "regex", label: "Regex" },
        { value: "simple", label: "Simple" },
    ];


    const handleApiCall = async (args) => {
        setIsLoading(true);
        setError(null);
        setApiData(null);
        try {
            let result;
            switch (searchType) {
                case "text_search":
                    result = await rcsbTextSearchPdb({ query_string: textQuery, ...args });
                    break;
                case "pdb_id_search":
                    result = await rcsbGetProteinDetailsByIdPypdb({ pdb_id_string: pdbId, ...args });
                    break;
                case "attribute_search":
                    const val = operator === 'in' ? value.split(',').map(s => s.trim()) : value;
                    result = await rcsbAttributeSearchPdb({ attribute_path: attributePath, operator, value: val, ...args });
                    break;
                case "combined_search":
                    const processedFilters = attributeFilters.map(f => ({
                        ...f,
                        value: f.operator === 'in' ? f.value.split(',').map(s => s.trim()) : f.value
                    })).filter(f => f.attribute_path && f.value); // Ensure filters are valid
                    if (processedFilters.length === 0) {
                        setError("Please provide at least one valid attribute filter for combined search.");
                        setIsLoading(false);
                        return;
                    }
                    result = await rcsbCombinedTextAndAttributeSearch({ text_query_string: combinedTextQuery, attribute_filters: processedFilters, ...args });
                    break;
                case "sequence_identity_search":
                    result = await rcsbSequenceIdentitySearch({ sequence, identity_cutoff: parseFloat(identityCutoff), e_value_cutoff: parseFloat(eValueCutoff), sequence_type: sequenceType, ...args });
                    break;
                case "sequence_motif_search":
                    result = await rcsbSequenceMotifSearch({ motif_pattern: motifPattern, pattern_type: patternType, sequence_type: sequenceType, ...args });
                    break;
                default:
                    throw new Error("Invalid search type");
            }
            setApiData(result);
        } catch (err) {
            setError(err.message || "An error occurred during the search.");
            console.error("RCSB PDB Search Error:", err);
        } finally {
            setIsLoading(false);
        }
    };
    
    const debouncedApiCall = useCallback(debounce(handleApiCall, 500), [searchType, textQuery, pdbId, attributePath, operator, value, combinedTextQuery, attributeFilters, sequence, identityCutoff, eValueCutoff, sequenceType, motifPattern, patternType]);

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
            case "text_search":
                if (!textQuery.trim()) {
                    setError("Text query cannot be empty.");
                    valid = false;
                }
                break;
            case "pdb_id_search":
                if (!pdbId.trim()) {
                    setError("PDB ID cannot be empty.");
                    valid = false;
                }
                break;
            case "attribute_search":
                if (!attributePath.trim() || !operator.trim() || (value === '' || value === undefined)) {
                    setError("Attribute path, operator, and value are required for attribute search.");
                    valid = false;
                }
                break;
            case "combined_search":
                if (!combinedTextQuery.trim()) {
                     setError("Text query is required for combined search.");
                     valid = false;
                }
                const validFilters = attributeFilters.filter(f => f.attribute_path && f.value);
                if (validFilters.length === 0) {
                    setError("At least one valid attribute filter is required for combined search.");
                    valid = false;
                }
                break;
            case "sequence_identity_search":
                if (!sequence.trim()) {
                    setError("Sequence cannot be empty.");
                    valid = false;
                }
                if (isNaN(parseFloat(identityCutoff)) || parseFloat(identityCutoff) < 0 || parseFloat(identityCutoff) > 1) {
                    setError("Identity cutoff must be a number between 0.0 and 1.0.");
                    valid = false;
                }
                if (isNaN(parseFloat(eValueCutoff))) {
                    setError("E-value cutoff must be a number.");
                    valid = false;
                }
                break;
            case "sequence_motif_search":
                if (!motifPattern.trim()) {
                    setError("Motif pattern cannot be empty.");
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
    
    const handleAddAttributeFilter = () => {
        setAttributeFilters([...attributeFilters, { attribute_path: "", operator: "exact_match", value: "" }]);
    };

    const handleRemoveAttributeFilter = (index) => {
        const newFilters = attributeFilters.filter((_, i) => i !== index);
        setAttributeFilters(newFilters);
    };

    const handleAttributeFilterChange = (index, field, val) => {
        const newFilters = attributeFilters.map((filter, i) => {
            if (i === index) {
                return { ...filter, [field]: val };
            }
            return filter;
        });
        setAttributeFilters(newFilters);
    };


    const renderSearchInputs = () => {
        switch (searchType) {
            case "text_search":
                return (
                    <div className="search-inputs-container">
                        <div className="input-field-container">
                            <label className="input-field-label">Text Query:</label>
                            <input
                                type="text"
                                value={textQuery}
                                onChange={(e) => setTextQuery(e.target.value)}
                                placeholder="Enter text query (e.g., hemoglobin)"
                                className="rcsb-input-field"
                            />
                        </div>
                    </div>
                );
            case "pdb_id_search":
                return (
                    <div className="search-inputs-container">
                        <div className="input-field-container">
                            <label className="input-field-label">PDB ID:</label>
                            <input
                                type="text"
                                value={pdbId}
                                onChange={(e) => setPdbId(e.target.value)}
                                placeholder="Enter PDB ID (e.g., 6M0J, 1TIM)"
                                className="rcsb-input-field"
                            />
                        </div>
                    </div>
                );
            case "attribute_search":
                return (
                    <div className="search-inputs-container">
                        <div className="input-field-container">
                            <label className="input-field-label">Attribute Path:</label>
                            <input
                                type="text"
                                value={attributePath}
                                onChange={(e) => setAttributePath(e.target.value)}
                                placeholder="Attribute Path (e.g., exptl.method)"
                                className="rcsb-input-field"
                            />
                        </div>
                        <GenericDropdownSelector
                            value={operator}
                            onChange={setOperator}
                            options={operatorOptions}
                            placeholder="Select Operator"
                        />
                        <div className="input-field-container">
                            <label className="input-field-label">Value:</label>
                            <input
                                type="text"
                                value={value}
                                onChange={(e) => setValue(e.target.value)}
                                placeholder="Value (e.g., X-RAY DIFFRACTION or list for 'in')"
                                className="rcsb-input-field"
                            />
                        </div>
                    </div>
                );
            case "combined_search":
                return (
                    <div className="search-inputs-container">
                        <div className="input-field-container">
                            <label className="input-field-label">Main Text Query:</label>
                            <input
                                type="text"
                                value={combinedTextQuery}
                                onChange={(e) => setCombinedTextQuery(e.target.value)}
                                placeholder="Main text query"
                                className="rcsb-input-field"
                            />
                        </div>
                        {attributeFilters.map((filter, index) => (
                            <div key={index} className="attribute-filter-group">
                                <h4>Attribute Filter {index + 1}</h4>
                                <div className="input-field-container">
                                    <label className="input-field-label">Attribute Path:</label>
                                    <input
                                        type="text"
                                        value={filter.attribute_path}
                                        onChange={(e) => handleAttributeFilterChange(index, 'attribute_path', e.target.value)}
                                        placeholder="Attribute Path"
                                        className="rcsb-input-field"
                                    />
                                </div>
                                <GenericDropdownSelector
                                    value={filter.operator}
                                    onChange={(val) => handleAttributeFilterChange(index, 'operator', val)}
                                    options={operatorOptions}
                                    placeholder="Select Operator"
                                />
                                <div className="input-field-container">
                                    <label className="input-field-label">Value:</label>
                                    <input
                                        type="text"
                                        value={filter.value}
                                        onChange={(e) => handleAttributeFilterChange(index, 'value', e.target.value)}
                                        placeholder="Value"
                                        className="rcsb-input-field"
                                    />
                                </div>
                                {attributeFilters.length > 1 && (
                                    <button type="button" onClick={() => handleRemoveAttributeFilter(index)} className="remove-filter-btn">Remove Filter</button>
                                )}
                            </div>
                        ))}
                        <button type="button" onClick={handleAddAttributeFilter} className="add-filter-btn">Add Attribute Filter</button>
                    </div>
                );
            case "sequence_identity_search":
                return (
                    <div className="search-inputs-container">
                        <div className="input-field-container">
                            <label className="input-field-label">Sequence:</label>
                            <textarea
                                value={sequence}
                                onChange={(e) => setSequence(e.target.value)}
                                placeholder="Protein, DNA, or RNA sequence string"
                                rows={4}
                                className="rcsb-textarea-field"
                            />
                        </div>
                        <GenericDropdownSelector
                            value={sequenceType}
                            onChange={setSequenceType}
                            options={sequenceTypeOptions}
                            placeholder="Select Sequence Type"
                        />
                        <div className="input-field-container">
                            <label className="input-field-label">Identity Cutoff:</label>
                            <input
                                type="number"
                                value={identityCutoff}
                                onChange={(e) => setIdentityCutoff(e.target.value)}
                                placeholder="Identity Cutoff (0.0 - 1.0, default 0.9)"
                                step="0.01"
                                min="0"
                                max="1"
                                className="rcsb-input-field"
                            />
                        </div>
                        <div className="input-field-container">
                            <label className="input-field-label">E-value Cutoff:</label>
                            <input
                                type="number"
                                value={eValueCutoff}
                                onChange={(e) => setEValueCutoff(e.target.value)}
                                placeholder="E-value Cutoff (default 1.0)"
                                step="0.1"
                                className="rcsb-input-field"
                            />
                        </div>
                    </div>
                );
            case "sequence_motif_search":
                return (
                    <div className="search-inputs-container">
                        <div className="input-field-container">
                            <label className="input-field-label">Motif Pattern:</label>
                            <input
                                type="text"
                                value={motifPattern}
                                onChange={(e) => setMotifPattern(e.target.value)}
                                placeholder="Motif pattern (e.g., PROSITE format)"
                                className="rcsb-input-field"
                            />
                        </div>
                        <GenericDropdownSelector
                            value={patternType}
                            onChange={setPatternType}
                            options={patternTypeOptions}
                            placeholder="Select Pattern Type"
                        />
                         <GenericDropdownSelector
                            value={sequenceType} // Reusing sequenceType state and options
                            onChange={setSequenceType}
                            options={sequenceTypeOptions}
                            placeholder="Select Sequence Type"
                        />
                    </div>
                );
            default:
                return <p>Select a search type.</p>;
        }
    };

    return (
        <div className="rcsb-pdb-explorer-container">
            {/* Header or Title if any */}
            {/* ... */}

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
                    {isLoading ? 'Searching...' : 'Search PDB'}
                </button>
            </form>

            {error && <p style={{ color: 'red', marginTop: '10px' }}>Error: {error}</p>}

            <DataViewer data={apiData} /> {/* Using DataViewer */}
        </div>
    );
};

RCSBPDBExplorer.propTypes = {
    initialSearchType: PropTypes.oneOf([
        "text_search",
        "pdb_id_search",
        "attribute_search",
        "combined_search",
        "sequence_identity_search",
        "sequence_motif_search",
    ]),
    toolData: PropTypes.object, // For pre-filling data if needed
};

export default RCSBPDBExplorer;

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
