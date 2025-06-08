import React, { useState, useEffect } from "react";
import { motion } from "framer-motion";
import { fadeInUpVariantStatic } from "../../../../components/animations/framerAnim";
import { getCidsByName, getCidsBySmiles, getCidsByInchikey, getCompoundByCid } from "../../../../services/api/mcpToolsService";
import GlassyContainer from "../../../../components/glassy_container/gc";
import { FaCheckDouble } from "react-icons/fa6";
import { IoWarningOutline } from "react-icons/io5";

const CompoundSearch = ({ toolData = null, initialSearchValue = "", initialSearchType = "name" }) => {
    const [searchValue, setSearchValue] = useState(initialSearchValue);
    const [searchType, setSearchType] = useState(initialSearchType);
    const [isLoading, setIsLoading] = useState(false);
    const [results, setResults] = useState(toolData);
    const [error, setError] = useState(null);
    const [isValidInput, setIsValidInput] = useState(false);

    // Handle initial results from props
    useEffect(() => {
        if (toolData) {
            setResults(toolData);
        }
    }, [toolData]);

    // Handle initial search value from props
    useEffect(() => {
        if (initialSearchValue) {
            setSearchValue(initialSearchValue);
            setIsValidInput(true);
        }
    }, [initialSearchValue]);

    // Validate input in real-time
    const validateInput = (value) => {
        if (!value || value.trim() === '') return false;
        
        switch (searchType) {
            case "name":
                return value.trim().length > 1;
            case "smiles":
                // Basic SMILES validation - at least contains some structure characters
                return /[A-Za-z0-9()=\[\]#@+-]/.test(value.trim());
            case "inchikey":
                // InChI Key format validation
                return /^[A-Z]{14}-[A-Z]{10}-[A-Z]$/.test(value.trim());
            default:
                return false;
        }
    };

    const handleInputChange = (value) => {
        setSearchValue(value);
        setIsValidInput(validateInput(value));
        if (error) setError(null);
    };    const handleFormSubmit = async () => {
        if (!searchValue.trim()) {
            setError('Please enter a search value');
            return;
        }
        
        if (!isValidInput) {
            setError('Please enter a valid search value');
            return;
        }
        
        setIsLoading(true);
        setError(null);
        setResults(null);
        
        try {
            let response;
            const trimmedValue = searchValue.trim();
            
            switch (searchType) {
                case "name":
                    response = await getCidsByName({ name: trimmedValue });
                    break;
                case "smiles":
                    response = await getCidsBySmiles({ smiles: trimmedValue });
                    break;
                case "inchikey":
                    response = await getCidsByInchikey({ inchikey: trimmedValue });
                    break;
                default:
                    throw new Error("Invalid search type");
            }

            if (response && response.status === "success") {
                const data = JSON.parse(response.result["0"]);
                
                if (data.result && Array.isArray(data.result) && data.result.length > 0) {
                    // Get detailed compound information for the first few CIDs
                    const cidList = data.result.slice(0, 10); // Limit to first 10 for performance
                    const compoundDetails = [];
                    
                    for (const cid of cidList) {
                        try {
                            const compoundResponse = await getCompoundByCid({ cid: parseInt(cid, 10) });
                            if (compoundResponse && compoundResponse.status === "success") {
                                const compoundData = JSON.parse(compoundResponse.result["0"]);
                                if (compoundData.result) {
                                    compoundDetails.push({
                                        cid: cid,
                                        ...compoundData.result
                                    });
                                }
                            }
                        } catch (compoundError) {
                            console.warn(`Failed to get details for CID ${cid}:`, compoundError);
                        }
                    }
                    
                    setResults({
                        search_type: searchType,
                        search_value: trimmedValue,
                        total_found: data.result.length,
                        cids: data.result,
                        compound_details: compoundDetails
                    });
                } else {
                    setError(`No compounds found for ${searchType}: "${trimmedValue}"`);
                }
            } else {
                setError(response.message || `Failed to search for ${searchType}`);
            }
            
        } catch (err) {
            console.error('Error searching compounds:', err);
            setError(err.message || 'Failed to search compounds');
        } finally {
            setIsLoading(false);
        }
    };

    const getPlaceholderText = () => {
        switch (searchType) {
            case "name":
                return "Enter compound name (e.g., aspirin, caffeine)";
            case "smiles":
                return "Enter SMILES string (e.g., CCO, C1=CC=CC=C1)";
            case "inchikey":
                return "Enter InChI Key (e.g., BSYNRYMUTXBXSQ-UHFFFAOYSA-N)";
            default:
                return "Enter search value";
        }
    };

    const getHeaderText = () => {
        switch (searchType) {
            case "name":
                return "Search by Compound Name";
            case "smiles":
                return "Search by SMILES";
            case "inchikey":
                return "Search by InChI Key";
            default:
                return "Compound Search";
        }
    };

    const handleReset = () => {
        setResults(null);
        setSearchValue("");
        setError(null);
        setIsValidInput(false);
    };

    return (
        <motion.div
            variants={fadeInUpVariantStatic}
            initial="initial"
            animate="animate"
            className="utility-component-container"
        >
            <div className="utility-input-section">
                <GlassyContainer>
                    <h3 style={{ marginBottom: '15px', fontWeight: '700' }}>{getHeaderText()}</h3>
                    <p style={{ marginBottom: '15px', color: 'var(--color-text-secondary)' }}>
                        Search for compounds in PubChem using different identifiers.
                    </p>
                    
                    {/* Search Type Selection */}
                    <div style={{ marginBottom: '15px' }}>
                        <label style={{ 
                            display: 'block', 
                            marginBottom: '8px', 
                            color: 'var(--color-text-secondary)', 
                            fontWeight: '500',
                            fontSize: '0.9em'
                        }}>
                            Search Type:
                        </label>
                        <select
                            value={searchType}
                            onChange={(e) => {
                                setSearchType(e.target.value);
                                setIsValidInput(validateInput(searchValue));
                            }}
                            disabled={isLoading}
                            style={{
                                width: '100%',
                                padding: '12px 16px',
                                background: 'var(--glassy-color)',
                                border: '1px solid var(--c-light-border)',
                                borderRadius: '8px',
                                color: 'var(--color-text-primary)',
                                fontSize: '14px',
                                fontWeight: '500'
                            }}
                        >
                            <option value="name">Compound Name</option>
                            <option value="smiles">SMILES</option>
                            <option value="inchikey">InChI Key</option>
                        </select>
                    </div>

                    {/* Add inline styles for input placeholder */}
                    <style>
                        {`
                            .search-input::placeholder {
                                color: var(--color-text-secondary);
                                font-weight: 400;
                            }
                            .search-input:focus {
                                outline: none;
                                background-color: var(--glassy-color);
                            }
                        `}
                    </style>
                    
                    {/* Custom Search Input with Validation */}
                    <div style={{ display: 'flex', flexDirection: 'column', gap: '10px' }}>
                        <div style={{ display: 'flex', alignItems: 'center', gap: '10px', width: '100%' }}>
                            <input
                                type="text"
                                className="search-input"
                                value={searchValue}
                                onChange={(e) => handleInputChange(e.target.value)}
                                onKeyDown={(e) => {
                                    if (e.key === 'Enter' && isValidInput && !isLoading) {
                                        handleFormSubmit();
                                    }
                                }}
                                placeholder={getPlaceholderText()}
                                disabled={isLoading}
                                style={{
                                    flexGrow: 1,
                                    padding: '12px',
                                    height: 'auto',
                                    borderRadius: '15px',
                                    border: 'none',
                                    outline: 'none',
                                    backgroundColor: 'var(--glassy-color)',
                                    color: 'var(--color-text-primary)',
                                    fontSize: '16px',
                                    fontWeight: '600',
                                    minHeight: '30px',
                                    transition: 'background-color 0.2s ease'
                                }}
                                onFocus={(e) => {
                                    e.target.style.backgroundColor = 'var(--glassy-color)';
                                }}
                                onBlur={(e) => {
                                    e.target.style.backgroundColor = 'var(--glassy-color)';
                                }}
                            />
                            {/* Validation Icon */}
                            {isValidInput ? (
                                <FaCheckDouble style={{ fontSize: '2rem', color: 'var(--color-success)' }} />
                            ) : (
                                <IoWarningOutline style={{ 
                                    fontSize: '2rem', 
                                    color: searchValue.length > 0 ? 'var(--color-alert)' : 'var(--glassy-color)' 
                                }} />
                            )}
                            <button
                                onClick={handleFormSubmit}
                                disabled={!isValidInput || isLoading}
                                style={{
                                    padding: '12px 20px',
                                    borderRadius: '15px',
                                    border: 'none',
                                    backgroundColor: (!isValidInput || isLoading) ? 'var(--color-disabled, #6c757d)' : 'var(--color-success, #28a745)',
                                    color: '#fff',
                                    cursor: (!isValidInput || isLoading) ? 'not-allowed' : 'pointer',
                                    fontSize: '1em',
                                    fontWeight: '500',
                                    minHeight: '54px',
                                    whiteSpace: 'nowrap'
                                }}
                            >
                                {isLoading ? 'Searching...' : 'Search Compounds'}
                            </button>
                        </div>
                        {/* Validation message */}
                        {searchValue.length > 0 && !isValidInput && (
                            <p style={{ 
                                color: 'var(--color-alert, #ff6b6b)', 
                                fontSize: '0.9em', 
                                margin: '5px 0 0 0' 
                            }}>
                                Please enter a valid {searchType === 'inchikey' ? 'InChI Key format' : searchType}
                            </p>
                        )}
                    </div>
                    {(results || error) && (
                        <button 
                            onClick={handleReset}
                            style={{
                                marginTop: '10px',
                                padding: '8px 16px',
                                backgroundColor: 'var(--color-secondary)',
                                color: 'var(--color-text-primary)',
                                border: 'none',
                                borderRadius: '8px',
                                cursor: 'pointer',
                                transition: 'all 0.2s ease',
                                fontSize: '0.9em'
                            }}
                        >
                            Reset
                        </button>
                    )}
                </GlassyContainer>
            </div>

            {error && (
                <div className="utility-results-section" style={{ marginTop: '20px' }}>
                    <GlassyContainer>
                        <div style={{
                            padding: '16px',
                            backgroundColor: 'rgba(255, 71, 87, 0.1)',
                            borderRadius: '8px',
                            border: '1px solid rgba(255, 71, 87, 0.3)'
                        }}>
                            <h4 style={{ marginBottom: '8px', color: 'var(--color-error, #ff6b6b)', margin: '0 0 8px 0' }}>
                                ✗ Error
                            </h4>
                            <p style={{ color: 'var(--color-error, #ff6b6b)', margin: 0 }}>
                                {error}
                            </p>
                        </div>
                    </GlassyContainer>
                </div>
            )}

            {results && (
                <div className="utility-results-section" style={{ marginTop: '20px' }}>
                    <GlassyContainer>
                        <h4 style={{ marginBottom: '15px', color: 'var(--color-text-primary)', fontWeight: '600' }}>
                            Search Results
                        </h4>
                        <div className="search-info" style={{
                            marginBottom: '20px',
                            padding: '12px',
                            background: 'var(--color-bg-secondary)',
                            borderRadius: '8px',
                            border: '1px solid var(--c-light-border)',
                            fontSize: '0.9em'
                        }}>
                            <span style={{ color: 'var(--color-text-secondary)' }}>Search Type: </span>
                            <strong style={{ color: 'var(--color-accent)' }}>{results.search_type}</strong>
                            <span style={{ color: 'var(--color-text-secondary)' }}> | Query: </span>
                            <strong style={{ color: 'var(--color-text-primary)' }}>"{results.search_value}"</strong>
                            <span style={{ color: 'var(--color-text-secondary)' }}> | Found: </span>
                            <strong style={{ color: 'var(--color-success)' }}>{results.total_found} compounds</strong>
                        </div>
                        
                        {results.compound_details && results.compound_details.length > 0 ? (
                            <div className="results-grid" style={{
                                display: 'grid',
                                gridTemplateColumns: 'repeat(auto-fill, minmax(250px, 1fr))',
                                gap: '16px',
                                marginBottom: '20px'
                            }}>
                                {results.compound_details.map((compound, index) => (
                                    <div key={compound.cid || index} className="result-item" style={{
                                        background: 'var(--color-bg-secondary)',
                                        border: '1px solid var(--c-light-border)',
                                        borderRadius: '12px',
                                        padding: '16px',
                                        textAlign: 'center',
                                        transition: 'all 0.3s ease'
                                    }}
                                    onMouseEnter={(e) => {
                                        e.target.style.transform = 'translateY(-2px)';
                                        e.target.style.boxShadow = '0 8px 20px rgba(0, 0, 0, 0.1)';
                                        e.target.style.borderColor = 'var(--color-accent)';
                                    }}
                                    onMouseLeave={(e) => {
                                        e.target.style.transform = 'translateY(0)';
                                        e.target.style.boxShadow = 'none';
                                        e.target.style.borderColor = 'var(--c-light-border)';
                                    }}>
                                        <div className="cid-info" style={{ marginBottom: '12px' }}>
                                            <strong style={{ color: 'var(--color-accent)', fontSize: '15px', fontWeight: '700' }}>
                                                CID: {compound.cid}
                                            </strong>
                                        </div>
                                        {compound.compound_name && (
                                            <div style={{ 
                                                marginBottom: '10px', 
                                                fontSize: '14px', 
                                                color: 'var(--color-text-primary)',
                                                fontWeight: '500'
                                            }}>
                                                {compound.compound_name}
                                            </div>
                                        )}
                                        {compound.cid && (
                                            <div className="compound-image" style={{ marginTop: '10px' }}>
                                                <img 
                                                    src={`https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/${compound.cid}/PNG?record_type=2d&image_size=small`}
                                                    alt={`Structure for CID ${compound.cid}`}
                                                    style={{
                                                        maxWidth: '100%',
                                                        height: 'auto',
                                                        borderRadius: '8px',
                                                        background: 'white',
                                                        padding: '8px',
                                                        border: '1px solid var(--c-light-border)'
                                                    }}
                                                    onError={(e) => {
                                                        e.target.style.display = 'none';
                                                    }}
                                                />
                                            </div>
                                        )}
                                    </div>
                                ))}
                            </div>
                        ) : (
                            <div style={{
                                textAlign: 'center',
                                color: 'var(--color-text-secondary)',
                                fontStyle: 'italic',
                                padding: '40px 20px',
                                background: 'var(--color-bg-secondary)',
                                border: '2px dashed var(--c-light-border)',
                                borderRadius: '12px'
                            }}>
                                Found {results.total_found} CIDs but no detailed compound information available.
                            </div>
                        )}
                        
                        {results.total_found > 10 && (
                            <div className="results-note" style={{
                                color: 'var(--color-text-secondary)',
                                fontStyle: 'italic',
                                textAlign: 'center',
                                marginTop: '16px',
                                padding: '12px',
                                background: 'var(--color-bg-secondary)',
                                borderRadius: '8px',
                                fontWeight: '500'
                            }}>
                                Showing first 10 compounds out of {results.total_found} found.
                            </div>
                        )}
                    </GlassyContainer>
                </div>
            )}
        </motion.div>
    );
};

export default CompoundSearch;
