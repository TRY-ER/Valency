import React, { useState, useEffect } from 'react';
import { motion } from 'framer-motion';
import { fadeInUpVariantStatic } from '../../../../components/animations/framerAnim';
import GlassyContainer from '../../../../components/glassy_container/gc';
import { FaCheckDouble } from 'react-icons/fa';
import { IoWarningOutline } from 'react-icons/io5';
import { fastIdentitySearchByCid } from '../../../../services/api/mcpToolsService';

const IdentitySearch = ({ toolData = null }) => {
  const [cid, setCid] = useState('');
  const [identityData, setIdentityData] = useState(null);
  const [isLoading, setIsLoading] = useState(false);
  const [error, setError] = useState('');
  const [isValidCid, setIsValidCid] = useState(false);
  const [searchQuery, setSearchQuery] = useState('');

  // Effect to handle initial data passed as props (similar to CompoundProperties pattern)
  useEffect(() => {
    console.log('tool data received >>', toolData);
  }, [toolData]);

  useEffect(() => {
    console.log("identity data updated >>", identityData);
  }, [identityData]);

  // Effect to handle initial data passed as props
  useEffect(() => {
    if (toolData) {
      console.log('Processing tool data:', toolData);
      
      // Try to extract identity data from different possible structures
      let extractedIdentityData = null;
      let extractedCid = null;
      
      // Case 1: toolData already has the correct IdentifierList structure
      if (toolData.CID && (toolData.SMILES || toolData.InChI || toolData.InChIKey)) {
        extractedIdentityData = toolData;
        extractedCid = toolData.CID;
      }
      // Case 2: toolData has IdentifierList structure
      else if (toolData.IdentifierList && toolData.IdentifierList.CID) {
        extractedIdentityData = toolData.IdentifierList;
        extractedCid = toolData.IdentifierList.CID;
      }
      // Case 3: toolData has result.IdentifierList structure (similar to API response)
      else if (toolData.result && toolData.result.IdentifierList && toolData.result.IdentifierList.CID) {
        extractedIdentityData = toolData.result.IdentifierList;
        extractedCid = toolData.result.IdentifierList.CID;
      }
      // Case 4: toolData might have a nested structure like the API response
      else if (toolData.result && toolData.result["0"]) {
        try {
          const parsed = typeof toolData.result["0"] === 'string' ? JSON.parse(toolData.result["0"]) : toolData.result["0"];
          if (parsed.result && parsed.result.IdentifierList && parsed.result.IdentifierList.CID) {
            extractedIdentityData = parsed.result.IdentifierList;
            extractedCid = parsed.result.IdentifierList.CID;
          }
        } catch (e) {
          console.error('Error parsing nested toolData:', e);
        }
      }
      // Case 5: toolData has CID property that we can use to populate the input
      else if (toolData.CID) {
        extractedCid = toolData.CID;
      }
      
      // Set the extracted data
      if (extractedIdentityData) {
        setIdentityData(extractedIdentityData);
        console.log('Set identity data from toolData:', extractedIdentityData);
      }
      
      if (extractedCid) {
        // Handle both single CID and array of CIDs
        const cidToSet = Array.isArray(extractedCid) ? extractedCid[0] : extractedCid;
        setCid(String(cidToSet));
        setIsValidCid(true);
        console.log('Set CID from toolData:', cidToSet);
      }
      
      // Clear any existing errors when loading from toolData
      setError('');
    }
  }, [toolData]);

  const validateCid = (value) => {
    const cidValue = parseInt(value.trim());
    return !isNaN(cidValue) && cidValue > 0;
  };

  const handleCidChange = (value) => {
    setCid(value);
    setIsValidCid(validateCid(value));
    if (error) setError('');
  };

  const handleIdentitySearch = async () => {
    if (!cid.trim()) {
      setError('Please enter a Compound ID (CID).');
      return;
    }

    const cidValue = parseInt(cid.trim());
    if (isNaN(cidValue) || cidValue <= 0) {
      setError('Please enter a valid positive CID number');
      return;
    }

    setIsLoading(true);
    setError('');
    setIdentityData(null);

    try {
      const result = await fastIdentitySearchByCid({ cid: cid.trim() });
      console.log("result", result);
      if (result && result.status === "success") {
        let processible = result.result["0"];
        processible = JSON.parse(processible);
        console.log("processible", processible); 
        if (processible && processible.result.IdentifierList && processible.result.IdentifierList.CID) {
          setIdentityData(processible.result.IdentifierList);
        } else {
          setError('No identity data found for the given CID');
        }
      } else {
        setError(result.message || 'No identity data found for the given CID');
      }
    } catch (err) {
      console.error('Error fetching identity data:', err);
      setError('No identity data found for the given CID.');
    } finally {
      setIsLoading(false);
    }
  };

  const resetForm = () => {
    setCid('');
    setIdentityData(null);
    setError('');
    setIsValidCid(false);
    setSearchQuery('');
  };

  const handleKeyPress = (e) => {
    if (e.key === 'Enter' && isValidCid && !isLoading) {
      handleIdentitySearch();
    }
  };

  // Get all identifier entries for filtering
  const getIdentifierEntries = () => {
    if (!identityData) return [];
    
    const entries = [];
    Object.keys(identityData).forEach(key => {
      const value = identityData[key];
      if (value && key !== 'CID') {
        if (Array.isArray(value)) {
          value.forEach(item => {
            entries.push({ type: key, value: item.toString() });
          });
        } else {
          entries.push({ type: key, value: value.toString() });
        }
      }
    });
    
    // Add CID at the beginning - but only show the number, not "CID" prefix
    if (identityData.CID) {
      const cidValues = Array.isArray(identityData.CID) ? identityData.CID : [identityData.CID];
      cidValues.forEach(cidValue => {
        entries.unshift({ type: 'CID', value: cidValue.toString(), displayValue: cidValue.toString() });
      });
    }
    
    return entries;
  };

  // Filter identity data based on search query
  const filteredIdentifiers = getIdentifierEntries().filter(entry =>
    entry.value.toLowerCase().includes(searchQuery.toLowerCase()) ||
    entry.type.toLowerCase().includes(searchQuery.toLowerCase())
  );

  return (
    <motion.div 
      className="utility-component-container"
      variants={fadeInUpVariantStatic}
      initial="initial"
      animate="animate"
    >
      <div className="utility-input-section">
        <GlassyContainer>
          <h3 style={{ marginBottom: '15px', fontWeight: '700' }}>Identity Search</h3>
          <p style={{ marginBottom: '20px', color: 'var(--color-text-secondary)' }}>
            Search for compound identity information using PubChem CID from IdentifierList structure.
          </p>

          {/* Add inline styles for input placeholder */}
          <style>
            {`
              .cid-input::placeholder {
                color: var(--color-text-secondary);
                font-weight: 400;
              }
              .cid-input:focus {
                outline: none;
                background-color: var(--glassy-color);
              }
            `}
          </style>

          {/* Custom CID Input with Validation */}
          <div style={{ display: 'flex', flexDirection: 'column', gap: '10px' }}>
            <div style={{ display: 'flex', alignItems: 'center', gap: '20px', width: '100%' }}>
              <h4 style={{ fontWeight: '600', color: 'var(--color-text-primary)', margin: 0, minWidth: '180px' }}>
                Enter CID:
              </h4>
              <div style={{ display: 'flex', alignItems: 'center', gap: '10px', flex: 1 }}>
                <input
                  type="number"
                  className="cid-input"
                  value={cid}
                  onChange={(e) => handleCidChange(e.target.value)}
                  onKeyDown={handleKeyPress}
                  placeholder="Enter CID (e.g., 2244)"
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
                {isValidCid ? (
                  <FaCheckDouble style={{ fontSize: '2rem', color: 'var(--color-success)' }} />
                ) : (
                  <IoWarningOutline style={{
                    fontSize: '2rem',
                    color: cid.length > 0 ? 'var(--color-alert)' : 'var(--glassy-color)'
                  }} />
                )}
                <button
                  onClick={handleIdentitySearch}
                  disabled={!isValidCid || isLoading}
                  style={{
                    padding: '12px 20px',
                    borderRadius: '15px',
                    border: 'none',
                    backgroundColor: (!isValidCid || isLoading) ? 'var(--color-disabled, #6c757d)' : 'var(--color-success, #28a745)',
                    color: '#fff',
                    cursor: (!isValidCid || isLoading) ? 'not-allowed' : 'pointer',
                    fontSize: '1em',
                    fontWeight: '500',
                    minHeight: '54px',
                    whiteSpace: 'nowrap'
                  }}
                >
                  {isLoading ? 'Fetching...' : 'Get Identity'}
                </button>
              </div>
            </div>
            {/* Validation message for CID */}
            {cid.length > 0 && !isValidCid && (
              <p style={{
                color: 'var(--color-alert, #ff6b6b)',
                fontSize: '0.9em',
                margin: '5px 0 0 0'
              }}>
                Please enter a valid positive CID number
              </p>
            )}
          </div>
          
          {(identityData || error) && (
            <div style={{ display: 'flex', justifyContent: 'flex-end', marginTop: '20px' }}>
              <button
                onClick={resetForm}
                style={{
                  padding: '12px 24px',
                  backgroundColor: 'var(--color-error, #dc3545)',
                  color: '#fff',
                  border: 'none',
                  borderRadius: '12px',
                  cursor: 'pointer',
                  transition: 'all 0.2s ease',
                  fontSize: '0.95em',
                  fontWeight: '600',
                  boxShadow: '0 2px 6px rgba(0, 0, 0, 0.1)'
                }}
                onMouseEnter={(e) => {
                  e.target.style.backgroundColor = 'var(--color-error-dark, #c82333)';
                  e.target.style.transform = 'translateY(-1px)';
                  e.target.style.boxShadow = '0 4px 8px rgba(0, 0, 0, 0.15)';
                }}
                onMouseLeave={(e) => {
                  e.target.style.backgroundColor = 'var(--color-error, #dc3545)';
                  e.target.style.transform = 'translateY(0)';
                  e.target.style.boxShadow = '0 2px 6px rgba(0, 0, 0, 0.1)';
                }}
              >
                Reset All
              </button>
            </div>
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

      {identityData && (
        <div className="utility-results-section" style={{ marginTop: '20px' }}>
          <GlassyContainer>
            <div style={{ display: 'flex', justifyContent: 'space-between', alignItems: 'center', marginBottom: '20px' }}>
              <h4 style={{ color: 'var(--color-text-primary)', fontWeight: '600', margin: 0 }}>
                Identity Data for CID {cid} ({filteredIdentifiers.length} identifiers)
              </h4>
              {getIdentifierEntries().length > 5 && (
                <div style={{ display: 'flex', alignItems: 'center', gap: '10px', minWidth: '300px' }}>
                  <input
                    type="text"
                    placeholder="Search identifiers..."
                    value={searchQuery}
                    onChange={(e) => setSearchQuery(e.target.value)}
                    style={{
                      flex: 1,
                      padding: '8px 12px',
                      borderRadius: '10px',
                      border: 'none',
                      outline: 'none',
                      backgroundColor: 'var(--glassy-color)',
                      color: 'var(--color-text-primary)',
                      fontSize: '14px',
                      fontWeight: '500',
                      transition: 'background-color 0.2s ease'
                    }}
                    onFocus={(e) => {
                      e.target.style.backgroundColor = 'var(--color-bg-primary)';
                      e.target.style.boxShadow = '0 0 0 2px rgba(40, 167, 69, 0.3)';
                    }}
                    onBlur={(e) => {
                      e.target.style.backgroundColor = 'var(--glassy-color)';
                      e.target.style.boxShadow = 'none';
                    }}
                  />
                  {searchQuery && (
                    <button
                      onClick={() => setSearchQuery('')}
                      style={{
                        padding: '8px',
                        backgroundColor: 'var(--color-error, #dc3545)',
                        color: '#fff',
                        border: 'none',
                        borderRadius: '6px',
                        cursor: 'pointer',
                        fontSize: '12px',
                        fontWeight: '500',
                        transition: 'all 0.2s ease'
                      }}
                      onMouseEnter={(e) => {
                        e.target.style.backgroundColor = 'var(--color-error-dark, #c82333)';
                      }}
                      onMouseLeave={(e) => {
                        e.target.style.backgroundColor = 'var(--color-error, #dc3545)';
                      }}
                    >
                      Clear
                    </button>
                  )}
                </div>
              )}
            </div>
            
            {filteredIdentifiers.length > 0 ? (
              <>
                <div style={{
                  display: 'grid',
                  gridTemplateColumns: 'repeat(auto-fill, minmax(320px, 1fr))',
                  gap: '12px',
                  maxHeight: '500px',
                  overflowY: 'auto',
                  padding: '15px',
                  backgroundColor: 'var(--color-bg-secondary)',
                  borderRadius: '12px',
                  border: '1px solid var(--c-light-border)'
                }}>
                  {filteredIdentifiers.map((identifier, index) => (
                    <div 
                      key={index}
                      style={{
                        padding: '16px',
                        backgroundColor: 'var(--glassy-color)',
                        borderRadius: '12px',
                        border: '1px solid var(--c-light-border)',
                        cursor: 'pointer',
                        transition: 'all 0.3s ease',
                        fontSize: '14px',
                        lineHeight: '1.5',
                        wordBreak: 'break-word',
                        color: 'var(--color-text-primary)',
                        fontWeight: '500',
                        boxShadow: '0 2px 4px rgba(0, 0, 0, 0.05)'
                      }}
                    //   onMouseEnter={(e) => {
                    //     e.target.style.borderColor = 'var(--color-success, #28a745)';
                    //     e.target.style.transform = 'translateY(-3px)';
                    //     e.target.style.boxShadow = '0 6px 20px rgba(0, 0, 0, 0.15)';
                    //   }}
                    //   onMouseLeave={(e) => {
                    //     e.target.style.borderColor = 'var(--c-light-border)';
                    //     e.target.style.transform = 'translateY(0)';
                    //     e.target.style.boxShadow = '0 2px 4px rgba(0, 0, 0, 0.05)';
                    //   }}
                      onClick={() => {
                        navigator.clipboard.writeText(identifier.value);
                        // Optional: Add a toast notification here
                      }}
                      title="Click to copy to clipboard"
                    >
                      <div>{identifier.displayValue || identifier.value}</div>
                    </div>
                  ))}
                </div>
              </>
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
                No identifiers match your search "{searchQuery}"
              </div>
            )}
          </GlassyContainer>
        </div>
      )}

      {!isLoading && !identityData && cid && !error && (
        <div className="utility-results-section" style={{ marginTop: '20px' }}>
          <GlassyContainer>
            <div style={{
              textAlign: 'center',
              color: 'var(--color-text-secondary)',
              fontStyle: 'italic',
              padding: '40px 20px',
              background: 'var(--color-bg-secondary)',
              border: '2px dashed var(--c-light-border)',
              borderRadius: '12px'
            }}>
              No identity data found for the specified compound.
            </div>
          </GlassyContainer>
        </div>
      )}
    </motion.div>
  );
};

export default IdentitySearch;
