import React, { useState, useEffect, useCallback, useRef } from "react";
import { createPortal } from "react-dom";
import "./ApprovedDrugsViewer.css";
import SimpleInputBox from "../../../components/UI/SimpleInputBox/SimpleInputBox";
import DataViewer from "../../../components/UI/DataViewer/DataViewer";
import InfoBox from "../../../components/UI/InfoBox/InfoBox";
import TwoDViewer from "../../../components/UI/TwoDViewer/TwoDViewer";
import { motion } from "framer-motion";
import { fadeInUpVariantStatic } from "../../../components/animations/framerAnim";
import GlassyContainer from "../../../components/glassy_container/gc";
import { getApprovedDrugs } from "../../../services/api/mcpToolsService";

// In-house Sort Option Selector Component with Portal
const SortOptionSelector = ({ value, onChange, options, disabled = false }) => {
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
            if (isOpen && 
                dropdownRef.current && 
                !dropdownRef.current.contains(event.target) &&
                triggerRef.current && 
                !triggerRef.current.contains(event.target)) {
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

    const selectedOption = options.find(option => option.value === value);

    // Portal dropdown component
    const PortalDropdown = () => {
        if (!isOpen) return null;

        return createPortal(
            <div
                ref={dropdownRef}
                className="approved-drugs-dropdown-portal"
                style={{
                    position: 'fixed',
                    top: `${dropdownPosition.top}px`,
                    left: `${dropdownPosition.left}px`,
                    width: `${dropdownPosition.width}px`,
                    zIndex: 9999
                }}
            >
                <div className="approved-drugs-dropdown-options">
                    {options.map((option) => (
                        <div
                            key={option.value}
                            className={`approved-drugs-dropdown-option ${option.value === value ? 'selected' : ''}`}
                            onClick={() => handleOptionClick(option.value)}
                        >
                            {option.label}
                        </div>
                    ))}
                </div>
            </div>,
            document.body
        );
    };

    return (
        <>
            <div
                className={`approved-drugs-sort-selector ${disabled ? 'disabled' : ''}`}
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
                        {selectedOption ? selectedOption.label : 'Select sort option'}
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

// Sort options for the dropdown
const sortOptions = [
    { value: "default", label: "Default Order" },
    { value: "molecular_weight", label: "By Molecular Weight" }
];

const ApprovedDrugsViewer = ({
    toolData = null,
    initialSortType = "default",
    hideControls = false
}) => {
    const [sortType, setSortType] = useState(initialSortType);
    const [apiData, setApiData] = useState(toolData);
    const [isLoading, setIsLoading] = useState(false);
    const [error, setError] = useState(null);
    
    // States for drug list handling
    const [drugsList, setDrugsList] = useState([]);
    const [selectedDrug, setSelectedDrug] = useState(null);
    const [selectedIndex, setSelectedIndex] = useState(0);
    const [searchFilter, setSearchFilter] = useState("");
    const [filteredDrugs, setFilteredDrugs] = useState([]);

    useEffect(() => {
        console.log('tool data received >>', toolData);
    }, [toolData]);

    useEffect(() => {
        console.log("apiData updated >>", apiData);
    }, [apiData]);

    // Handle fetching approved drugs
    const handleFetchDrugs = useCallback(async () => {
        setIsLoading(true);
        setError(null);

        try {
            const args = sortType === "molecular_weight" ? { order_by_mw: true } : {};
            const result = await getApprovedDrugs(args);
            console.log('API result:', result); 
            if (result && result.status === "success") {
                const processible = JSON.parse(result.result["0"]);
                console.log('Approved drugs results:', processible);

                if (processible.result && Array.isArray(processible.result) && processible.result.length > 0) {
                    setDrugsList(processible.result);
                    setSelectedDrug(processible.result[0]);
                    setSelectedIndex(0);
                    setApiData(processible.result[0]);
                } else {
                    setError("No approved drugs found.");
                    setDrugsList([]);
                    setSelectedDrug(null);
                    setApiData(null);
                }
            } else {
                setError(result?.message || "Failed to fetch approved drugs.");
                setDrugsList([]);
                setSelectedDrug(null);
                setApiData(null);
            }
        } catch (error) {
            console.error('Error fetching approved drugs:', error);
            setError(error.message || "Failed to fetch approved drugs. Please try again.");
            setDrugsList([]);
            setSelectedDrug(null);
            setApiData(null);
        } finally {
            setIsLoading(false);
        }
    }, [sortType]);

    // Effect to handle initial data passed as props
    useEffect(() => {
        if (toolData) {
            setApiData(toolData);
        }
    }, [toolData]);

    // Auto-fetch on component mount or sort type change
    useEffect(() => {
        if (!toolData) {
            handleFetchDrugs();
        }
    }, [sortType, toolData, handleFetchDrugs]);

    // Handle search filter
    useEffect(() => {
        if (!searchFilter.trim()) {
            setFilteredDrugs(drugsList);
        } else {
            const filtered = drugsList.filter(drug => {
                const searchTerm = searchFilter.toLowerCase();
                return (
                    (drug.pref_name && drug.pref_name.toLowerCase().includes(searchTerm)) ||
                    (drug.molecule_chembl_id && drug.molecule_chembl_id.toLowerCase().includes(searchTerm)) ||
                    (drug.molecule_synonyms && drug.molecule_synonyms.some(syn => 
                        syn.molecule_synonym && syn.molecule_synonym.toLowerCase().includes(searchTerm)
                    ))
                );
            });
            setFilteredDrugs(filtered);
        }
    }, [drugsList, searchFilter]);

    // Keyboard navigation for drug selection
    useEffect(() => {
        const handleKeyDown = (event) => {
            if (filteredDrugs.length > 0) {
                if (event.key === 'ArrowDown') {
                    event.preventDefault();
                    const newIndex = (selectedIndex + 1) % filteredDrugs.length;
                    setSelectedIndex(newIndex);
                    setSelectedDrug(filteredDrugs[newIndex]);
                    setApiData(filteredDrugs[newIndex]);
                } else if (event.key === 'ArrowUp') {
                    event.preventDefault();
                    const newIndex = selectedIndex === 0 ? filteredDrugs.length - 1 : selectedIndex - 1;
                    setSelectedIndex(newIndex);
                    setSelectedDrug(filteredDrugs[newIndex]);
                    setApiData(filteredDrugs[newIndex]);
                }
            }
        };

        if (drugsList.length > 0) {
            window.addEventListener('keydown', handleKeyDown);
            return () => window.removeEventListener('keydown', handleKeyDown);
        }
    }, [filteredDrugs, selectedIndex]);

    const handleDrugSelect = (drug, index) => {
        setSelectedDrug(drug);
        setSelectedIndex(index);
        setApiData(drug);
    };

    const getHeaderText = () => {
        if (drugsList.length === 0) return "ChEMBL Approved Drugs";
        return `ChEMBL Approved Drugs (${drugsList.length} total${filteredDrugs.length !== drugsList.length ? `, ${filteredDrugs.length} filtered` : ''})`;
    };

    useEffect(() => {
        if (error) {
            const timer = setTimeout(() => {
                setError(null);
            }, 5000);
            return () => clearTimeout(timer);
        }
    }, [error]);

    return (
        <>
            <motion.div
                className="approved-drugs-container"
                variants={fadeInUpVariantStatic}
                initial="initial"
                animate="animate"
            >
                {/* Header and Controls */}
                {!hideControls && (
                    <div className="approved-drugs-row-1">
                        <GlassyContainer>
                            <div className="approved-drugs-search-controls">
                                <h2 className="approved-drugs-header">{getHeaderText()}</h2>
                                
                                <div className="approved-drugs-controls-row">
                                    {/* Sort Options */}
                                    <div className="approved-drugs-dropdown-section">
                                        <h3 className="approved-drugs-dropdown-header">Sort Options</h3>
                                        <div className="approved-drugs-dropdown">
                                            <SortOptionSelector
                                                value={sortType}
                                                onChange={setSortType}
                                                options={sortOptions}
                                                disabled={isLoading}
                                            />
                                        </div>
                                    </div>

                                    {/* Search Filter */}
                                    <div className="approved-drugs-search-section">
                                        <h3 className="approved-drugs-dropdown-header">Search Filter</h3>
                                        <SimpleInputBox
                                            placeholder="Filter by name, ChEMBL ID, or synonym..."
                                            value={searchFilter}
                                            onChange={setSearchFilter}
                                            disabled={isLoading || drugsList.length === 0}
                                        />
                                    </div>

                                    {/* Refresh Button */}
                                    <div className="approved-drugs-action-section">
                                        <button
                                            className="approved-drugs-refresh-btn"
                                            onClick={handleFetchDrugs}
                                            disabled={isLoading}
                                        >
                                            {isLoading ? 'Loading...' : 'Refresh Data'}
                                        </button>
                                    </div>
                                </div>
                            </div>
                        </GlassyContainer>
                    </div>
                )}

                {/* Loading State */}
                {isLoading && (
                    <div className="approved-drugs-loading">
                        <div className="loading-spinner"></div>
                        <p>Fetching approved drugs data...</p>
                    </div>
                )}

                {/* Error State */}
                {error && (
                    <div className="approved-drugs-error">
                        <p>⚠️ {error}</p>
                    </div>
                )}

                {/* Drugs List Display */}
                {filteredDrugs.length > 0 && (
                    <div className="approved-drugs-row-2">
                        <GlassyContainer>
                            <div className="approved-drugs-selection-container">
                                {/* Left Panel - Drug List */}
                                <div className="approved-drugs-list-panel">
                                    <div className="approved-drugs-list-header">
                                        <h4>Select a Drug ({filteredDrugs.length} results)</h4>
                                        <p className="approved-drugs-list-instruction">
                                            Click to select • Use ↑↓ arrow keys to navigate
                                        </p>
                                    </div>
                                    <div className="approved-drugs-list-container">
                                        {filteredDrugs.map((drug, index) => (
                                            <div
                                                key={drug.molecule_chembl_id || index}
                                                className={`approved-drugs-item ${
                                                    selectedDrug && selectedDrug.molecule_chembl_id === drug.molecule_chembl_id ? 'selected' : ''
                                                }`}
                                                onClick={() => handleDrugSelect(drug, index)}
                                            >
                                                <span className="approved-drugs-item-number">
                                                    {index + 1}.
                                                </span>
                                                <div className="approved-drugs-item-content">
                                                    <span className="approved-drugs-item-name">
                                                        {drug.pref_name || 'Unnamed Drug'}
                                                    </span>
                                                    <span className="approved-drugs-item-id">
                                                        {drug.molecule_chembl_id}
                                                    </span>
                                                    {drug.molecule_properties && drug.molecule_properties.mw_freebase && typeof drug.molecule_properties.mw_freebase === 'number' && (
                                                        <span className="approved-drugs-item-mw">
                                                            MW: {drug.molecule_properties.mw_freebase.toFixed(2)}
                                                        </span>
                                                    )}
                                                </div>
                                            </div>
                                        ))}
                                    </div>
                                </div>

                                {/* Right Panel - Info and 2D Viewer */}
                                <div className="approved-drugs-details-panel">
                                    {selectedDrug ? (
                                        <>
                                            {/* Selected drug header */}
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
                                                    Selected ({selectedIndex + 1}/{filteredDrugs.length}): <span style={{ 
                                                        fontFamily: 'monospace', 
                                                        color: 'var(--color-accent)',
                                                        fontWeight: '600' 
                                                    }}>
                                                        {selectedDrug.pref_name || selectedDrug.molecule_chembl_id || 'Unknown'}
                                                    </span>
                                                </h4>
                                            </div>
                                            
                                            {/* Show visualization only if we have SMILES data */}
                                            {selectedDrug.molecule_structures && selectedDrug.molecule_structures.canonical_smiles ? (
                                                <div className="approved-drugs-details-row">
                                                    {/* InfoBox */}
                                                    <div style={{ flex: '4' }}>
                                                        <InfoBox 
                                                            activeMol={selectedDrug.molecule_structures.canonical_smiles} 
                                                            isValidMol={true} 
                                                            infoType={"MOL"} 
                                                        />
                                                    </div>
                                                    
                                                    {/* TwoDViewer */}
                                                    <div style={{ flex: '6' }}>
                                                        <TwoDViewer 
                                                            activeMol={selectedDrug.molecule_structures.canonical_smiles} 
                                                            isValidMol={true} 
                                                            visType={"MOL"} 
                                                        />
                                                    </div>
                                                </div>
                                            ) : (
                                                <div className="approved-drugs-placeholder">
                                                    No structure data available for this drug
                                                </div>
                                            )}
                                        </>
                                    ) : (
                                        <div className="approved-drugs-placeholder">
                                            Select a drug from the left panel to view details
                                        </div>
                                    )}
                                </div>
                            </div>
                        </GlassyContainer>
                    </div>
                )}

                {/* DataViewer - show selected drug data */}
                {selectedDrug && (
                    <div className="approved-drugs-row-3" style={{ marginTop: '20px' }}>
                        <DataViewer
                            data={selectedDrug}
                            title={`Approved Drug Data - ${selectedDrug.pref_name || selectedDrug.molecule_chembl_id || 'Unknown'}`}
                            initiallyExpanded={true}
                        />
                    </div>
                )}
            </motion.div>
        </>
    );
};

export default ApprovedDrugsViewer;
