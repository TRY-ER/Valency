import React, { useState, useEffect, useCallback, useRef } from "react";
import { createPortal } from "react-dom";
import "./ChemBLActivityFetcher.css";
import SimpleInputBox from "../../../components/UI/SimpleInputBox/SimpleInputBox";
import DataViewer from "../../../components/UI/DataViewer/DataViewer";
import InfoBox from "../../../components/UI/InfoBox/InfoBox";
import TwoDViewer from "../../../components/UI/TwoDViewer/TwoDViewer";
import { motion } from "framer-motion";
import { fadeInUpVariantStatic } from "../../../components/animations/framerAnim";
import GlassyContainer from "../../../components/glassy_container/gc";
import { getActivitiesForTarget, getActivitiesForMolecule } from "../../../services/api/mcpToolsService";

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
const ActivitySearchTypeSelector = ({ value, onChange, options = [], disabled = false }) => {
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
                        padding: '8px 12px',
                        backgroundColor: 'var(--color-bg-primary)',
                        border: `1px solid ${isOpen ? 'var(--color-accent)' : 'var(--c-light-border)'}`,
                        borderRadius: isOpen ? '6px 6px 0 0' : '6px',
                        borderBottom: isOpen ? 'none' : `1px solid var(--c-light-border)`,
                        cursor: disabled ? 'not-allowed' : 'pointer',
                        transition: 'all 0.2s ease',
                        minHeight: '40px',
                        opacity: disabled ? '0.6' : '1'
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
                            e.target.style.backgroundColor = 'var(--color-bg-primary)';
                        }
                    }}
                >
                    <span className="selected-text" style={{
                        color: 'var(--color-text-primary)',
                        fontSize: '0.85rem',
                        fontWeight: '500'
                    }}>
                        {selectedOption ? selectedOption.label : 'Select activity search type'}
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

// Activity search type options
const activitySearchTypeOptions = [
    {
        value: "target",
        label: "Target Activities",
        description: "Get activities for a biological target using ChEMBL target ID"
    },
    {
        value: "molecule",
        label: "Molecule Activities",
        description: "Get activities for a molecule using ChEMBL molecule ID"
    }
];

// Standard type options for target activities
const standardTypeOptions = [
    { value: "", label: "All", description: "Fetch all standard types" },
    { value: "IC50", label: "IC50", description: "Half maximal inhibitory concentration" },
    { value: "Ki", label: "Ki", description: "Inhibition constant" },
    { value: "Kd", label: "Kd", description: "Dissociation constant" },
    { value: "EC50", label: "EC50", description: "Half maximal effective concentration" },
    { value: "AC50", label: "AC50", description: "Half maximal activity concentration" },
    { value: "CC50", label: "CC50", description: "Half maximal cytotoxic concentration" },
    { value: "GI50", label: "GI50", description: "Half maximal growth inhibition" },
    { value: "LC50", label: "LC50", description: "Half lethal concentration" },
    { value: "LD50", label: "LD50", description: "Half lethal dose" },
    { value: "MIC", label: "MIC", description: "Minimum inhibitory concentration" },
    { value: "Potency", label: "Potency", description: "General potency measurement" },
    { value: "Activity", label: "Activity", description: "General activity measurement" }
];

const ChemBLActivityFetcher = ({
    toolData = null,
    initialSearchValue = "",
    initialSearchType = "target",
    hideInputBox = false
}) => {
    const [searchValue, setSearchValue] = useState(initialSearchValue);
    const [inputValue, setInputValue] = useState(initialSearchValue);
    const [searchType, setSearchType] = useState(initialSearchType);
    const [standardType, setStandardType] = useState("IC50");
    const [pchemblValueExists, setPchemblValueExists] = useState(true);
    const [apiData, setApiData] = useState(toolData);
    const [isLoading, setIsLoading] = useState(false);
    const [error, setError] = useState(null);
    const [validationError, setValidationError] = useState(null);
    
    // States for target search with activity molecules
    const [targetData, setTargetData] = useState(null);
    const [activityMolecules, setActivityMolecules] = useState([]);
    const [selectedActivityMolecule, setSelectedActivityMolecule] = useState(null);
    const [selectedActivityIndex, setSelectedActivityIndex] = useState(0);
    
    // States for molecule search with activity data
    const [moleculeActivities, setMoleculeActivities] = useState([]);
    const [selectedMoleculeActivity, setSelectedMoleculeActivity] = useState(null);
    const [selectedMoleculeIndex, setSelectedMoleculeIndex] = useState(0);

    // Ref to track if search type is being set from toolData
    const isSettingFromToolData = useRef(false);

    useEffect(() => {
        console.log('activity tool data received >>', toolData);
        if (toolData && toolData.type && toolData.content) {
            console.log('Processing activity tool data:', toolData);
            
            // Set flag to indicate we're setting search type from toolData
            isSettingFromToolData.current = true;
            
            // Set search type based on tool data type
            setSearchType(toolData.type);

            console.log('activity tool data >>', toolData);
            
            if (toolData.type === "target") {
                // Handle target search results - expecting structure with target_data and activity_data
                if (toolData.content && typeof toolData.content === 'object' && 
                    toolData.content.target_data && toolData.content.activity_data) {
                    console.log('Setting target activity data from tool data:', toolData.content);
                    setTargetData(toolData.content.target_data);
                    setActivityMolecules(Array.isArray(toolData.content.activity_data) ? toolData.content.activity_data : []);
                    setSelectedActivityMolecule(Array.isArray(toolData.content.activity_data) && toolData.content.activity_data.length > 0 ? toolData.content.activity_data[0] : null);
                    setSelectedActivityIndex(0);
                    setApiData(null); // Clear regular data display
                    setMoleculeActivities([]);
                    setSelectedMoleculeActivity(null);
                    setSelectedMoleculeIndex(0);
                    // Clear any existing errors for successful data processing
                    setError(null);
                } else {
                    // No target activities found
                    setError(`No target activities found for ${toolData.type} search.`);
                    setTargetData(null);
                    setActivityMolecules([]);
                    setSelectedActivityMolecule(null);
                    setApiData(null);
                    setMoleculeActivities([]);
                    setSelectedMoleculeActivity(null);
                }
            } else if (toolData.type === "molecule") {
                // Handle molecule search results - expecting array of activity data
                if (Array.isArray(toolData.content) && toolData.content.length > 0) {
                    console.log('Setting molecule activity data from tool data:', toolData.content);
                    setMoleculeActivities(toolData.content);
                    setSelectedMoleculeActivity(toolData.content[0]);
                    setSelectedMoleculeIndex(0);
                    setApiData(null); // Clear regular data display
                    setTargetData(null);
                    setActivityMolecules([]);
                    setSelectedActivityMolecule(null);
                    setSelectedActivityIndex(0);
                    // Clear any existing errors for successful data processing
                    setError(null);
                } else if (toolData.content) {
                    // If it's a single activity, treat it as a single-item array for consistent UI
                    setMoleculeActivities([toolData.content]);
                    setSelectedMoleculeActivity(toolData.content);
                    setSelectedMoleculeIndex(0);
                    setApiData(null); // Clear regular data display
                    setTargetData(null);
                    setActivityMolecules([]);
                    setSelectedActivityMolecule(null);
                    setSelectedActivityIndex(0);
                    // Clear any existing errors for successful data processing
                    setError(null);
                } else {
                    // No molecule activities found
                    setError(`No molecule activities found for ${toolData.type} search.`);
                    setMoleculeActivities([]);
                    setSelectedMoleculeActivity(null);
                    setApiData(null);
                    setTargetData(null);
                    setActivityMolecules([]);
                    setSelectedActivityMolecule(null);
                }
            }
            
            // Reset flag after state updates
            setTimeout(() => {
                isSettingFromToolData.current = false;
            }, 0);
        } else if (toolData) {
            // Fallback for old format toolData without type/content structure
            console.log('Processing legacy activity tool data format:', toolData);
            
            // Extract the data array from the nested toolData structure if it exists
            // Structure: toolData.results.0.data (array)
            let extractedData = null;
            if (toolData.result && toolData.result["0"]) {
                extractedData = JSON.parse(toolData.result["0"])
                extractedData = extractedData.data || extractedData; // Fallback to data if available
            } else {
                // Fallback to original toolData if structure is different
                extractedData = toolData;
            }
            
            // Check if this is target search data with target_data and activity_data
            if (extractedData && typeof extractedData === 'object' && 
                extractedData.target_data && extractedData.activity_data) {
                // Handle target search structure
                setTargetData(extractedData.target_data);
                setActivityMolecules(Array.isArray(extractedData.activity_data) ? extractedData.activity_data : []);
                setSelectedActivityMolecule(Array.isArray(extractedData.activity_data) && extractedData.activity_data.length > 0 ? extractedData.activity_data[0] : null);
                setSelectedActivityIndex(0);
                setApiData(null); // Clear regular data display
                setMoleculeActivities([]);
                setSelectedMoleculeActivity(null);
                setSelectedMoleculeIndex(0);
            } else if (Array.isArray(extractedData) && extractedData.length > 0) {
                // Handle molecule search structure (array of activity data)
                setMoleculeActivities(extractedData);
                setSelectedMoleculeActivity(extractedData[0]);
                setSelectedMoleculeIndex(0);
                setApiData(null); // Clear regular data display
                setTargetData(null);
                setActivityMolecules([]);
                setSelectedActivityMolecule(null);
                setSelectedActivityIndex(0);
            } else {
                // Handle regular single data structure
                setApiData(extractedData);
                setTargetData(null);
                setActivityMolecules([]);
                setSelectedActivityMolecule(null);
                setSelectedActivityIndex(0);
                setMoleculeActivities([]);
                setSelectedMoleculeActivity(null);
                setSelectedMoleculeIndex(0);
            }
            
            setError(null);
        }
    }, [toolData]);

    useEffect(() => {
        console.log('activity api data updated >>', apiData);
        console.log('target data >>', targetData);
        console.log('activity molecules >>', activityMolecules);
        console.log('molecule activities >>', moleculeActivities);
    }, [apiData, targetData, activityMolecules, moleculeActivities]);

    // Handle input changes with real-time validation
    const handleInputChange = useCallback((value) => {
        setInputValue(value);
        setValidationError(null);

        if (value.trim() === '') {
            setValidationError(null);
            return;
        }

        const trimmedValue = value.trim();

        // Enhanced validation for molecule search (ChEMBL ID required)
        if (searchType === "molecule") {
            if (!validateChemblId(trimmedValue)) {
                const errorMessage = trimmedValue.toLowerCase().includes('chembl') 
                    ? "Invalid ChEMBL ID format. Ensure it follows the pattern: CHEMBL followed by numbers (e.g., CHEMBL25, CHEMBL1234567)."
                    : "Invalid ChEMBL ID format. ChEMBL ID must start with 'CHEMBL' followed by numbers (e.g., CHEMBL25, CHEMBL1234567).";
                setValidationError(errorMessage);
            }
        }
        // For target search, we don't validate format since it's just a target name
    }, [searchType]);

    const handleSubmit = useCallback(async (value) => {
        console.log('Submitting search for:', value, 'with type:', searchType);
        const trimmedValue = value.trim();

        if (!trimmedValue) {
            if (searchType === "molecule") {
                setError('Please enter a ChEMBL ID');
            } else {
                setError('Please enter a target name');
            }
            return;
        }

        if (searchType === "molecule") {
            if (!validateChemblId(trimmedValue)) {
                const errorMessage = trimmedValue.toLowerCase().includes('chembl') 
                    ? "Invalid ChEMBL ID format. Ensure it follows the pattern: CHEMBL followed by numbers (e.g., CHEMBL25, CHEMBL1234567)."
                    : "Invalid ChEMBL ID format. ChEMBL ID must start with 'CHEMBL' followed by numbers (e.g., CHEMBL25, CHEMBL1234567).";
                setValidationError(errorMessage);
                return;
            }
        }

        setIsLoading(true);
        setError(null);
        setApiData(null);
        setTargetData(null);
        setActivityMolecules([]);
        setSelectedActivityMolecule(null);
        setSelectedActivityIndex(0);
        setMoleculeActivities([]);
        setSelectedMoleculeActivity(null);
        setSelectedMoleculeIndex(0);
        setSearchValue(trimmedValue);

        try {
            let result;

            if (searchType === "target") {
                result = await getActivitiesForTarget({
                    biological_target_name: trimmedValue,
                    strandard_type: standardType
                });
            } else if (searchType === "molecule") {
                result = await getActivitiesForMolecule({
                    molecule_chembl_id: trimmedValue,
                    pchembl_value_exists: pchemblValueExists
                });
            }

            // Extract the data array from the nested response structure
            // Response structure: result.results.0.data (array)
            let extractedData = null;
            console.log('API result:', JSON.parse(result.result["0"]));
            if (result && result.result && result.result["0"]) {
                console.log("this part is getting triggered")
                extractedData = JSON.parse(result.result["0"]);
                console.log('Extracted data:', extractedData);
                extractedData = extractedData.data || extractedData; // Fallback to data if available
            } else {
                // Fallback to original result if structure is different
                extractedData = result;
            }

            // Check if this is target search data with target_data and activity_data
            if (extractedData && typeof extractedData === 'object' && 
                extractedData.target_data && extractedData.activity_data) {
                // Handle target search structure
                setTargetData(extractedData.target_data);
                setActivityMolecules(Array.isArray(extractedData.activity_data) ? extractedData.activity_data : []);
                setSelectedActivityMolecule(Array.isArray(extractedData.activity_data) && extractedData.activity_data.length > 0 ? extractedData.activity_data[0] : null);
                setSelectedActivityIndex(0);
                setApiData(null); // Clear regular data display
                setMoleculeActivities([]);
                setSelectedMoleculeActivity(null);
                setSelectedMoleculeIndex(0);
            } else if (Array.isArray(extractedData) && extractedData.length > 0) {
                // Handle molecule search structure (array of activity data)
                setMoleculeActivities(extractedData);
                setSelectedMoleculeActivity(extractedData[0]);
                setSelectedMoleculeIndex(0);
                setApiData(null); // Clear regular data display
                setTargetData(null);
                setActivityMolecules([]);
                setSelectedActivityMolecule(null);
                setSelectedActivityIndex(0);
            } else {
                // Handle regular single data structure
                setApiData(extractedData);
                setTargetData(null);
                setActivityMolecules([]);
                setSelectedActivityMolecule(null);
                setSelectedActivityIndex(0);
                setMoleculeActivities([]);
                setSelectedMoleculeActivity(null);
                setSelectedMoleculeIndex(0);
            }

            setError(null);
        } catch (err) {
            console.error('Activity fetch error:', err);
            setError(err.message || 'Failed to fetch activity data');
            setApiData(null);
            setTargetData(null);
            setActivityMolecules([]);
            setSelectedActivityMolecule(null);
            setSelectedActivityIndex(0);
            setMoleculeActivities([]);
            setSelectedMoleculeActivity(null);
            setSelectedMoleculeIndex(0);
        } finally {
            setIsLoading(false);
        }
    }, [searchType, standardType, pchemblValueExists]);

    // Keyboard navigation for activity molecules and molecule activities
    useEffect(() => {
        const handleKeyDown = (event) => {
            // Handle activity molecules navigation (target search)
            if (activityMolecules && activityMolecules.length > 0) {
                if (event.key === 'ArrowDown') {
                    event.preventDefault();
                    const nextIndex = (selectedActivityIndex + 1) % activityMolecules.length;
                    setSelectedActivityIndex(nextIndex);
                    setSelectedActivityMolecule(activityMolecules[nextIndex]);
                } else if (event.key === 'ArrowUp') {
                    event.preventDefault();
                    const prevIndex = selectedActivityIndex === 0 ? activityMolecules.length - 1 : selectedActivityIndex - 1;
                    setSelectedActivityIndex(prevIndex);
                    setSelectedActivityMolecule(activityMolecules[prevIndex]);
                }
                return;
            }
            
            // Handle molecule activities navigation (molecule search)
            if (moleculeActivities && moleculeActivities.length > 0) {
                if (event.key === 'ArrowDown') {
                    event.preventDefault();
                    const nextIndex = (selectedMoleculeIndex + 1) % moleculeActivities.length;
                    setSelectedMoleculeIndex(nextIndex);
                    setSelectedMoleculeActivity(moleculeActivities[nextIndex]);
                } else if (event.key === 'ArrowUp') {
                    event.preventDefault();
                    const prevIndex = selectedMoleculeIndex === 0 ? moleculeActivities.length - 1 : selectedMoleculeIndex - 1;
                    setSelectedMoleculeIndex(prevIndex);
                    setSelectedMoleculeActivity(moleculeActivities[prevIndex]);
                }
            }
        };

        window.addEventListener('keydown', handleKeyDown);
        return () => window.removeEventListener('keydown', handleKeyDown);
    }, [selectedActivityIndex, activityMolecules, selectedMoleculeIndex, moleculeActivities]);

    // Effect to handle initial search value
    useEffect(() => {
        if (initialSearchValue && !toolData) {
            handleSubmit(initialSearchValue);
        }
    }, [initialSearchValue, toolData, handleSubmit]);

    // Clear error message after some time
    useEffect(() => {
        if (error) {
            const timer = setTimeout(() => {
                setError(null);
            }, 10000); // Clear error after 10 seconds
            return () => clearTimeout(timer);
        }
    }, [error]);

    const getPlaceholderText = () => {
        if (searchType === "target") {
            return "Enter target name (e.g., hERG)";
        } else {
            return "Enter molecule ChEMBL ID (e.g., CHEMBL25)";
        }
    };

    const getHeaderText = () => {
        if (searchType === "target") {
            return "ChEMBL Target Activity Fetcher";
        } else {
            return "ChEMBL Molecule Activity Fetcher";
        }
    };

    // Reset validation error when search type changes (but not when set from toolData)
    // useEffect(() => {
    //     // Only reset if the search type change is user-initiated (not from toolData)
    //     if (!isSettingFromToolData.current) {
    //         setValidationError(null);
    //         setError(null);
    //         setApiData(null);
    //         setTargetData(null);
    //         setActivityMolecules([]);
    //         setSelectedActivityMolecule(null);
    //         setSelectedActivityIndex(0);
    //         setMoleculeActivities([]);
    //         setSelectedMoleculeActivity(null);
    //         setSelectedMoleculeIndex(0);
    //     }
    // }, [searchType]);

    return (
        <motion.div
            className="activity-fetcher-container"
            variants={fadeInUpVariantStatic}
            initial="hidden"
            animate="visible"
        >
            <GlassyContainer>
                <div className="activity-fetcher-row-1">
                    <h3 style={{ marginBottom: '16px', fontWeight: '700' }}>
                        {getHeaderText()}
                    </h3>

                    <div className="activity-fetcher-search-controls">
                        {/* Row 1: Activity Search Type and Additional Options */}
                        <div className="activity-fetcher-controls-row-1">
                            {/* Activity Search Type Dropdown */}
                            <div className="activity-fetcher-dropdown-section">
                                <h4 className="activity-fetcher-dropdown-header">Search Type:</h4>
                                <div className="activity-fetcher-dropdown">
                                    <ActivitySearchTypeSelector
                                        value={searchType}
                                        onChange={setSearchType}
                                        options={activitySearchTypeOptions}
                                        disabled={isLoading}
                                    />
                                </div>
                            </div>

                            {/* Standard Type for Target Activities */}
                            {searchType === "target" && (
                                <div className="activity-fetcher-dropdown-section">
                                    <h4 className="activity-fetcher-dropdown-header">Standard Type:</h4>
                                    <div className="activity-fetcher-dropdown">
                                        <ActivitySearchTypeSelector
                                            value={standardType}
                                            onChange={setStandardType}
                                            options={standardTypeOptions}
                                            disabled={isLoading}
                                        />
                                    </div>
                                </div>
                            )}

                            {/* pChEMBL Checkbox for Molecule Activities */}
                            {searchType === "molecule" && (
                                <div className="activity-fetcher-checkbox-section">
                                    <h4 className="activity-fetcher-dropdown-header">Options:</h4>
                                    <label className="checkbox-label">
                                        <input
                                            type="checkbox"
                                            checked={pchemblValueExists}
                                            onChange={(e) => setPchemblValueExists(e.target.checked)}
                                            disabled={isLoading}
                                        />
                                        <span className="checkbox-text">Only activities with pChEMBL values</span>
                                    </label>
                                </div>
                            )}
                        </div>

                        {/* Row 2: Search Input */}
                        {!hideInputBox && (
                            <div className="activity-fetcher-controls-row-2">
                                <div className="activity-fetcher-input-section">
                                    <h4 className="activity-fetcher-dropdown-header">Search Value:</h4>
                                    <SimpleInputBox
                                        value={inputValue}
                                        onChange={handleInputChange}
                                        onSubmit={handleSubmit}
                                        header=""
                                        placeholder={getPlaceholderText()}
                                        buttonText={isLoading ? 'Fetching...' : 'Fetch Activities'}
                                        isLoading={isLoading}
                                        error={validationError}
                                        disabled={isLoading}
                                    />
                                </div>
                            </div>
                        )}
                    </div>
                </div>

                <div className="activity-fetcher-row-2">
                    {isLoading && (
                        <div className="activity-fetcher-loading" style={{ marginTop: '20px', textAlign: 'center' }}>
                            <GlassyContainer>
                                <div style={{ padding: '20px' }}>
                                    <div style={{
                                        border: '4px solid #f3f3f3',
                                        borderTop: '4px solid rgb(3, 196, 3)',
                                        borderRadius: '50%',
                                        width: '40px',
                                        height: '40px',
                                        animation: 'spin 2s linear infinite',
                                        margin: '0 auto'
                                    }}></div>
                                    <p style={{ marginTop: '10px', color: 'var(--color-text-secondary)' }}>
                                        Fetching activity data...
                                    </p>
                                </div>
                            </GlassyContainer>
                        </div>
                    )}

                    {error && (
                        <div className="activity-fetcher-error" style={{ marginTop: '20px' }}>
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
                                            Fetch Error
                                        </h4>
                                        <p style={{
                                            color: '#ff3b30',
                                            margin: '0',
                                            fontSize: '0.9rem',
                                            lineHeight: '1.4'
                                        }}>
                                            {error}
                                        </p>
                                    </div>
                                </div>
                            </GlassyContainer>
                        </div>
                    )}

                    {/* Target Data Display - show target information when target search is performed */}
                    {targetData && !isLoading && !error && (
                        <div className="activity-fetcher-target-data" style={{ marginTop: '20px' }}>
                            <DataViewer
                                data={targetData}
                                title={`Target Information${searchValue ? ` for "${searchValue}"` : ''}${standardType ? ` (${standardType})` : ''}`}
                                initiallyExpanded={false}
                            />
                        </div>
                    )}

                    {/* Activity Molecules Results - show activity molecules with toggleable selection */}
                    {activityMolecules && activityMolecules.length > 0 && !isLoading && !error && (
                        <div className="activity-fetcher-activity-molecules" style={{ marginTop: '20px' }}>
                            <GlassyContainer>
                                <h3 style={{ marginBottom: '12px', fontWeight: '700' }}>
                                    Activity Molecules Results ({activityMolecules.length} found)
                                </h3>
                                
                                <div className="activity-molecules-results-container">
                                    {/* Left Panel - Activity Molecules List */}
                                    <div className="activity-molecules-candidates-panel">
                                        <h4 style={{
                                            color: 'var(--color-text-primary)',
                                            fontSize: '1rem',
                                            marginBottom: '10px',
                                            fontWeight: '500'
                                        }}>
                                            Activity Molecules
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
                                            {activityMolecules.map((molecule, index) => (
                                                <div 
                                                    key={index}
                                                    onClick={() => {
                                                        setSelectedActivityMolecule(molecule);
                                                        setSelectedActivityIndex(index);
                                                    }}
                                                    className={`activity-molecule-candidate-item ${selectedActivityMolecule === molecule ? 'selected' : ''}`}
                                                >
                                                    <span className="activity-molecule-candidate-number">
                                                        {index + 1}.
                                                    </span>
                                                    <div className="activity-molecule-candidate-info">
                                                        <span className="activity-molecule-candidate-id">
                                                            {molecule.molecule_chembl_id || molecule.canonical_smiles || 'Unknown ID'}
                                                        </span>
                                                        {molecule.standard_value && molecule.standard_units && (
                                                            <span className="activity-molecule-candidate-activity">
                                                                {molecule.standard_value} {molecule.standard_units}
                                                            </span>
                                                        )}
                                                    </div>
                                                </div>
                                            ))}
                                        </div>
                                    </div>

                                    {/* Right Panel - Info and 2D Viewer */}
                                    <div className="activity-molecules-details-panel">
                                        {selectedActivityMolecule ? (
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
                                                        Selected ({selectedActivityIndex + 1}/{activityMolecules.length}): <span style={{ 
                                                            fontFamily: 'monospace', 
                                                            color: 'var(--color-accent)',
                                                            fontWeight: '600' 
                                                        }}>
                                                            {selectedActivityMolecule.molecule_chembl_id || selectedActivityMolecule.canonical_smiles || 'Unknown'}
                                                        </span>
                                                        {selectedActivityMolecule.standard_value && selectedActivityMolecule.standard_units && (
                                                            <span style={{
                                                                color: 'var(--color-success)',
                                                                fontWeight: '600',
                                                                marginLeft: '8px'
                                                            }}>
                                                                ({selectedActivityMolecule.standard_value} {selectedActivityMolecule.standard_units})
                                                            </span>
                                                        )}
                                                    </h4>
                                                </div>
                                                
                                                {/* Show visualization only if we have SMILES data */}
                                                {selectedActivityMolecule.canonical_smiles ? (
                                                    <div className="activity-molecules-details-row">
                                                        {/* InfoBox */}
                                                        <div style={{ flex: '4' }}>
                                                            <InfoBox 
                                                                key={`activity-infobox-${selectedActivityIndex}`}
                                                                activeMol={selectedActivityMolecule.canonical_smiles} 
                                                                isValidMol={true} 
                                                                infoType={"MOL"} 
                                                            />
                                                        </div>
                                                        
                                                        {/* TwoDViewer */}
                                                        <div style={{ flex: '6' }}>
                                                            <TwoDViewer 
                                                                key={`activity-2dviewer-${selectedActivityIndex}`}
                                                                activeMol={selectedActivityMolecule.canonical_smiles} 
                                                                isValidMol={true} 
                                                                visType={"MOL"} 
                                                            />
                                                        </div>
                                                    </div>
                                                ) : (
                                                    <div className="activity-molecules-placeholder">
                                                        No SMILES structure data available for this molecule
                                                    </div>
                                                )}
                                            </>
                                        ) : (
                                            <div className="activity-molecules-placeholder">
                                                Select an activity molecule from the left panel to view details
                                            </div>
                                        )}
                                    </div>
                                </div>
                            </GlassyContainer>
                        </div>
                    )}

                    {/* Selected Activity Molecule Data Viewer */}
                    {selectedActivityMolecule && (
                        <div className="activity-fetcher-selected-molecule-data" style={{ marginTop: '20px' }}>
                            <DataViewer
                                key={`activity-molecule-data-${selectedActivityIndex}`}
                                data={selectedActivityMolecule}
                                title={`Activity Molecule Data${searchValue ? ` for "${searchValue}"` : ''}${
                                    selectedActivityMolecule ? 
                                    ` - ${selectedActivityMolecule.molecule_chembl_id || selectedActivityMolecule.canonical_smiles || 'Unknown'}` : 
                                    ''
                                }`}
                                initiallyExpanded={false}
                            />
                        </div>
                    )}

                    {/* Molecule Activities Results - show molecule activities with toggleable selection */}
                    {moleculeActivities && moleculeActivities.length > 0 && !isLoading && !error && (
                        <div className="activity-fetcher-molecule-activities" style={{ marginTop: '20px' }}>
                            <GlassyContainer>
                                <h3 style={{ marginBottom: '12px', fontWeight: '700' }}>
                                    Molecule Activities Results ({moleculeActivities.length} found)
                                </h3>
                                
                                <div className="activity-molecules-results-container">
                                    {/* Left Panel - Molecule Activities List */}
                                    <div className="activity-molecules-candidates-panel">
                                        <h4 style={{
                                            color: 'var(--color-text-primary)',
                                            fontSize: '1rem',
                                            marginBottom: '10px',
                                            fontWeight: '500'
                                        }}>
                                            Activity Targets 
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
                                            {moleculeActivities.map((activity, index) => (
                                                <div 
                                                    key={index}
                                                    onClick={() => {
                                                        setSelectedMoleculeActivity(activity);
                                                        setSelectedMoleculeIndex(index);
                                                    }}
                                                    className={`activity-molecule-candidate-item ${selectedMoleculeActivity === activity ? 'selected' : ''}`}
                                                >
                                                    <span className="activity-molecule-candidate-number">
                                                        {index + 1}.
                                                    </span>
                                                    <div className="activity-molecule-candidate-info">
                                                        <span className="activity-molecule-candidate-id">
                                                            {activity.target_chembl_id || activity.target_pref_name || 'Unknown Target'}
                                                        </span>
                                                        {activity.standard_value && activity.standard_units && activity.standard_type && (
                                                            <span className="activity-molecule-candidate-activity">
                                                                {activity.standard_type}: {activity.standard_value} {activity.standard_units}
                                                            </span>
                                                        )}
                                                    </div>
                                                </div>
                                            ))}
                                        </div>
                                    </div>

                                    {/* Right Panel - Info and 2D Viewer */}
                                    <div className="activity-molecules-details-panel">
                                        {selectedMoleculeActivity ? (
                                            <>
                                                {/* Selected activity header */}
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
                                                        Selected ({selectedMoleculeIndex + 1}/{moleculeActivities.length}): <span style={{ 
                                                            fontFamily: 'monospace', 
                                                            color: 'var(--color-accent)',
                                                            fontWeight: '600' 
                                                        }}>
                                                            {selectedMoleculeActivity.target_chembl_id || selectedMoleculeActivity.target_pref_name || 'Unknown'}
                                                        </span>
                                                        {selectedMoleculeActivity.standard_value && selectedMoleculeActivity.standard_units && selectedMoleculeActivity.standard_type && (
                                                            <span style={{
                                                                color: 'var(--color-success)',
                                                                fontWeight: '600',
                                                                marginLeft: '8px'
                                                            }}>
                                                                ({selectedMoleculeActivity.standard_type}: {selectedMoleculeActivity.standard_value} {selectedMoleculeActivity.standard_units})
                                                            </span>
                                                        )}
                                                    </h4>
                                                </div>
                                                
                                                {/* Show visualization only if we have SMILES data */}
                                                {selectedMoleculeActivity.canonical_smiles ? (
                                                    <div className="activity-molecules-details-row">
                                                        {/* InfoBox */}
                                                        <div style={{ flex: '4' }}>
                                                            <InfoBox 
                                                                key={`molecule-infobox-${selectedMoleculeIndex}`}
                                                                activeMol={selectedMoleculeActivity.canonical_smiles} 
                                                                isValidMol={true} 
                                                                infoType={"MOL"} 
                                                            />
                                                        </div>
                                                        
                                                        {/* TwoDViewer */}
                                                        <div style={{ flex: '6' }}>
                                                            <TwoDViewer 
                                                                key={`molecule-2dviewer-${selectedMoleculeIndex}`}
                                                                activeMol={selectedMoleculeActivity.canonical_smiles} 
                                                                isValidMol={true} 
                                                                visType={"MOL"} 
                                                            />
                                                        </div>
                                                    </div>
                                                ) : (
                                                    <div className="activity-molecules-placeholder">
                                                        No SMILES structure data available for this activity record
                                                    </div>
                                                )}
                                            </>
                                        ) : (
                                            <div className="activity-molecules-placeholder">
                                                Select an activity record from the left panel to view details
                                            </div>
                                        )}
                                    </div>
                                </div>
                            </GlassyContainer>
                        </div>
                    )}

                    {/* Selected Molecule Activity Data Viewer */}
                    {selectedMoleculeActivity && (
                        <div className="activity-fetcher-selected-activity-data" style={{ marginTop: '20px' }}>
                            <DataViewer
                                key={`molecule-activity-data-${selectedMoleculeIndex}`}
                                data={selectedMoleculeActivity}
                                title={`Molecule Activity Data${searchValue ? ` for "${searchValue}"` : ''}${
                                    selectedMoleculeActivity ? 
                                    ` - ${selectedMoleculeActivity.target_chembl_id || selectedMoleculeActivity.target_pref_name || 'Unknown'}` : 
                                    ''
                                }`}
                                initiallyExpanded={false}
                            />
                        </div>
                    )}

                    {/* DataViewer - show fetched activity data */}
                    {apiData && !isLoading && !error && (
                        <div className="activity-fetcher-row-3" style={{ marginTop: '20px' }}>
                            <DataViewer
                                data={apiData}
                                title={`${getHeaderText()} Results${searchValue ? ` for "${searchValue}"` : ''}${searchType === "target" ? ` (${standardType})` :
                                        searchType === "molecule" ? ` (pChEMBL ${pchemblValueExists ? 'enabled' : 'disabled'})` :
                                            ''
                                    }`}
                                initiallyExpanded={false}
                            />
                        </div>
                    )}
                </div>
            </GlassyContainer>
        </motion.div>
    );
};

export default ChemBLActivityFetcher;
