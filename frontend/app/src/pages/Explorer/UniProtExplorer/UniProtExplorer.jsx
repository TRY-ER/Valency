import React, { useState, useEffect, useRef, useCallback } from 'react';
import { createPortal } from 'react-dom';
import PropTypes from 'prop-types';
import {
    uniprotSearchUniprotkb,
    uniprotGetUniprotkbEntry
} from '../../../services/api/mcpToolsService'; // Adjust path as needed
import DataViewer from '../../../components/UI/DataViewer/DataViewer'; // Adjusted path
import './UniProtExplorer.css'; // Import the new CSS file

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
                top: rect.bottom + window.scrollY,
                left: rect.left + window.scrollX,
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
            if (isOpen && 
                !triggerRef.current?.contains(event.target) && 
                !dropdownRef.current?.contains(event.target)) {
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
        const handleEscape = (event) => {
            if (event.key === 'Escape' && isOpen) {
                setIsOpen(false);
            }
        };

        if (isOpen) {
            document.addEventListener('keydown', handleEscape);
            return () => document.removeEventListener('keydown', handleEscape);
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
                            if (!disabled && (e.key === 'Enter' || e.key === ' ')) {
                                e.preventDefault();
                                setIsOpen(!isOpen);
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
                                e.target.style.backgroundColor = 'var(--color-bg-secondary)';
                                e.target.style.borderColor = 'var(--color-accent)';
                            }
                        }}
                        onMouseLeave={(e) => {
                            if (!disabled && !isOpen) {
                                e.target.style.backgroundColor = 'var(--glassy-color)';
                                e.target.style.borderColor = 'var(--c-light-border)';
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
                top: rect.bottom + window.scrollY,
                left: rect.left + window.scrollX,
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
            if (isOpen && 
                !triggerRef.current?.contains(event.target) && 
                !dropdownRef.current?.contains(event.target)) {
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
        const handleEscape = (event) => {
            if (event.key === 'Escape' && isOpen) {
                setIsOpen(false);
            }
        };

        if (isOpen) {
            document.addEventListener('keydown', handleEscape);
            return () => document.removeEventListener('keydown', handleEscape);
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
                            if (!disabled && (e.key === 'Enter' || e.key === ' ')) {
                                e.preventDefault();
                                setIsOpen(!isOpen);
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
                                e.target.style.backgroundColor = 'var(--color-bg-secondary)';
                                e.target.style.borderColor = 'var(--color-accent)';
                            }
                        }}
                        onMouseLeave={(e) => {
                            if (!disabled && !isOpen) {
                                e.target.style.backgroundColor = 'var(--glassy-color)';
                                e.target.style.borderColor = 'var(--c-light-border)';
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

const UniProtExplorer = ({ initialSearchType = "search_uniprotkb", toolData = null }) => {
    const [searchType, setSearchType] = useState(initialSearchType);
    const [loading, setLoading] = useState(false);
    const [results, setResults] = useState(null);
    const [error, setError] = useState(null);
    
    // Search parameters
    const [queryString, setQueryString] = useState('');
    const [uniprotId, setUniprotId] = useState('');
    const [resultFormat, setResultFormat] = useState('json'); // Fixed to JSON for frontend parsing
    const [size, setSize] = useState(30);
    const [includeIsoform, setIncludeIsoform] = useState(false);

    // Search type options
    const searchTypeOptions = [
        {
            value: "search_uniprotkb",
            label: "Search UniProtKB",
            description: "Search the UniProtKB database with query strings and filters"
        },
        {
            value: "get_uniprotkb_entry",
            label: "Get UniProtKB Entry",
            description: "Retrieve a specific UniProtKB entry by its ID"
        }
    ];

    // Fixed format to JSON for frontend parsing

    // Effect to handle initial data passed as props
    useEffect(() => {
        if (toolData) {
            console.log('Processing UniProt tool data:', toolData);
            
            try {
                let processedData = null;
                
                // Handle different possible toolData structures
                if (toolData.result && toolData.result["0"]) {
                    // If toolData comes from API response format
                    const parsedResult = JSON.parse(toolData.result["0"]);
                    if (parsedResult.data) {
                        processedData = parsedResult.data;
                    } else if (parsedResult.results) {
                        processedData = parsedResult.results;
                    } else {
                        processedData = parsedResult;
                    }
                } else if (toolData.data) {
                    // Direct data structure
                    processedData = toolData.data;
                } else if (toolData.results) {
                    // Results structure
                    processedData = toolData.results;
                } else if (Array.isArray(toolData) || (typeof toolData === 'object' && toolData !== null)) {
                    // Direct array or object data
                    processedData = toolData;
                }
                
                if (processedData) {
                    setResults(processedData);
                    setError(null);
                    console.log('Set initial UniProt results from toolData:', processedData);
                } else {
                    console.warn('Could not extract results from toolData:', toolData);
                }
            } catch (err) {
                console.error('Error processing UniProt toolData:', err);
                setError('Error processing initial data');
            }
        }
    }, [toolData]);

    const debouncedSearch = useCallback(
        debounce(async () => {
            if (!queryString.trim() && searchType === "search_uniprotkb") return;
            if (!uniprotId.trim() && searchType === "get_uniprotkb_entry") return;
            
            setLoading(true);
            setError(null);
            
            try {
                let result;
                if (searchType === "search_uniprotkb") {
                    result = await uniprotSearchUniprotkb({
                        query_string: queryString,
                        result_format: resultFormat,
                        size: size,
                        include_isoform: includeIsoform
                    });
                } else if (searchType === "get_uniprotkb_entry") {
                    result = await uniprotGetUniprotkbEntry({
                        uniprot_id: uniprotId,
                        result_format: resultFormat
                    });
                }
                if (result.status === "success") {
                    if (result.result){
                        var processible = result.result["0"];
                        processible = JSON.parse(processible);
                        if (processible.data){
                            setResults(processible.data);
                        }
                        else if (processible.results){
                            setResults(processible.results);
                        }
                        else {
                            setResults(processible);
                        }
                    }
                    else {
                        setResults(result);
                    }
                }
            } catch (err) {
                setError(err.message || 'An error occurred during the search');
                setResults(null);
            } finally {
                setLoading(false);
            }
        }, 500),
        [searchType, queryString, uniprotId, resultFormat, size, includeIsoform]
    );

    const handleSearch = () => {
        debouncedSearch();
    };

    const handleClear = () => {
        setQueryString('');
        setUniprotId('');
        setResults(null);
        setError(null);
        setSize(30);
        setIncludeIsoform(false);
    };

    const renderSearchFields = () => {
        switch (searchType) {
            case "search_uniprotkb":
                return (
                    <div className="search-fields">
                        <div className="field-group">
                            <label>Query String *</label>
                            <div className="input-container">
                                <input
                                    type="text"
                                    value={queryString}
                                    onChange={(e) => setQueryString(e.target.value)}
                                    placeholder="e.g., gene:BRCA1, organism:human, protein kinase"
                                    className="search-input"
                                />
                                <small className="field-hint">
                                    Enter a UniProt query string (gene names, protein names, organisms, etc.)
                                </small>
                            </div>
                        </div>
                        
                        <div className="field-group">
                            <label>Number of Results</label>
                            <div className="input-container">
                                <input
                                    type="number"
                                    value={size}
                                    onChange={(e) => setSize(parseInt(e.target.value) || 30)}
                                    min="1"
                                    max="500"
                                    className="search-input"
                                />
                            </div>
                        </div>
                        
                        <div className="field-group checkbox-group">
                            <label>Include Isoforms</label>
                            <div className="input-container">
                                <input
                                    type="checkbox"
                                    checked={includeIsoform}
                                    onChange={(e) => setIncludeIsoform(e.target.checked)}
                                />
                            </div>
                        </div>
                    </div>
                );
                
            case "get_uniprotkb_entry":
                return (
                    <div className="search-fields">
                        <div className="field-group">
                            <label>UniProt ID *</label>
                            <div className="input-container">
                                <input
                                    type="text"
                                    value={uniprotId}
                                    onChange={(e) => setUniprotId(e.target.value)}
                                    placeholder="e.g., P12345, SPIKE_SARS2"
                                    className="search-input"
                                />
                                <small className="field-hint">
                                    Enter a UniProtKB ID (accession number or entry name)
                                </small>
                            </div>
                        </div>
                    </div>
                );
                
            default:
                return null;
        }
    };

    return (
        <div className="uniprot-explorer-container">
            <div className="uniprot-search-type-selector-section">
                <ActivitySearchTypeSelector
                    value={searchType}
                    onChange={setSearchType}
                    options={searchTypeOptions}
                    placeholder="Search Type"
                />
            </div>

            <div className="search-section">
                {renderSearchFields()}
                
                <div className="search-buttons">
                    <button 
                        className={`search-button primary ${loading ? 'loading' : ''}`}
                        onClick={handleSearch}
                        disabled={loading || 
                            (searchType === "search_uniprotkb" && !queryString.trim()) ||
                            (searchType === "get_uniprotkb_entry" && !uniprotId.trim())
                        }
                    >
                        {loading ? 'Searching...' : 'Search'}
                    </button>
                    <button 
                        className="search-button secondary"
                        onClick={handleClear}
                        disabled={loading}
                    >
                        Clear
                    </button>
                </div>
            </div>

            {error && (
                <div className="error-section">
                    <div className="error-message">
                        <strong>Error:</strong> {error}
                    </div>
                </div>
            )}

            {results && (
                <div className="results-section">
                    <DataViewer 
                        data={results} 
                        title={`UniProt ${searchType === "search_uniprotkb" ? "Search" : "Entry"} Results`}
                    />
                </div>
            )}
        </div>
    );
};

UniProtExplorer.propTypes = {
    initialSearchType: PropTypes.oneOf([
        "search_uniprotkb",
        "get_uniprotkb_entry"
    ]),
    toolData: PropTypes.object,
};

export default UniProtExplorer;

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
