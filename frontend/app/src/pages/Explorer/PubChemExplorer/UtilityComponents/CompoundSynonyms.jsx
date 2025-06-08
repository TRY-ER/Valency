import React, { useState, useEffect } from 'react';
import { motion } from 'framer-motion';
import { fadeInUpVariantStatic } from '../../../../components/animations/framerAnim';
import GlassyContainer from '../../../../components/glassy_container/gc';
import { FaCheckDouble } from 'react-icons/fa';
import { IoWarningOutline } from 'react-icons/io5';
import { getCompoundSynonymsByCid } from '../../../../services/api/mcpToolsService';

const CompoundSynonyms = ({ toolData = null }) => {
  const [cid, setCid] = useState('');
  const [synonyms, setSynonyms] = useState([]);
  const [isLoading, setIsLoading] = useState(false);
  const [error, setError] = useState('');
  const [isValidCid, setIsValidCid] = useState(false);
  const [searchQuery, setSearchQuery] = useState('');

  // Effect to handle initial data passed as props (similar to CompoundProperties pattern)
  useEffect(() => {
    console.log('tool data received >>', toolData);
  }, [toolData]);

  useEffect(() => {
    console.log("synonyms updated >>", synonyms);
  }, [synonyms]);

  // Effect to handle initial data passed as props
  useEffect(() => {
    if (toolData) {
      console.log('Processing tool data:', toolData);
      
      // Try to extract synonyms from different possible structures
      let extractedSynonyms = null;
      let extractedCid = null;
      
      // Case 1: toolData already has the correct synonyms structure (array of strings)
      if (Array.isArray(toolData) && toolData.length > 0 && typeof toolData[0] === 'string') {
        extractedSynonyms = toolData;
      }
      // Case 2: toolData has InformationList structure
      else if (toolData.InformationList && toolData.InformationList.Information && 
               toolData.InformationList.Information[0] && toolData.InformationList.Information[0].Synonym) {
        extractedSynonyms = toolData.InformationList.Information[0].Synonym;
        extractedCid = toolData.InformationList.Information[0].CID || null;
      }
      // Case 3: toolData has result.InformationList structure (similar to API response)
      else if (toolData.result && toolData.result.InformationList && 
               toolData.result.InformationList.Information && 
               toolData.result.InformationList.Information[0] && 
               toolData.result.InformationList.Information[0].Synonym) {
        extractedSynonyms = toolData.result.InformationList.Information[0].Synonym;
        extractedCid = toolData.result.InformationList.Information[0].CID || null;
      }
      // Case 4: toolData might have a nested structure like the API response
      else if (toolData.result && toolData.result["0"]) {
        try {
          const parsed = typeof toolData.result["0"] === 'string' ? JSON.parse(toolData.result["0"]) : toolData.result["0"];
          if (parsed.result && parsed.result.InformationList && 
              parsed.result.InformationList.Information && 
              parsed.result.InformationList.Information[0] && 
              parsed.result.InformationList.Information[0].Synonym) {
            extractedSynonyms = parsed.result.InformationList.Information[0].Synonym;
            extractedCid = parsed.result.InformationList.Information[0].CID || null;
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
      if (extractedSynonyms && Array.isArray(extractedSynonyms)) {
        setSynonyms(extractedSynonyms);
        console.log('Set synonyms from toolData:', extractedSynonyms);
      }
      
      if (extractedCid) {
        setCid(String(extractedCid));
        setIsValidCid(true);
        console.log('Set CID from toolData:', extractedCid);
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

  const handleSynonymSearch = async () => {
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
    setSynonyms([]);

    try {
      const result = await getCompoundSynonymsByCid({ cid: cid.trim() });
      console.log("result >>", result)
      
      if (result && result.status === "success") {
        let processible = result.result["0"];
        processible = JSON.parse(processible);
        console.log("processible >>", processible)
        
        if (processible && processible.result.InformationList && processible.result.InformationList.Information && 
            processible.result.InformationList.Information[0] && processible.result.InformationList.Information[0].Synonym) {
          setSynonyms(processible.result.InformationList.Information[0].Synonym);
        } else {
          setSynonyms([]);
        }
      } else {
        setError(result.message || 'No synonyms found for the given CID');
        setSynonyms([]);
      }
    } catch (error) {
      console.error('Error fetching synonyms:', error);
      setError('No synonyms found for the given CID.');
      setSynonyms([]);
    } finally {
      setIsLoading(false);
    }
  };

  const resetForm = () => {
    setCid('');
    setSynonyms([]);
    setError('');
    setIsValidCid(false);
    setSearchQuery('');
  };

  const handleKeyPress = (e) => {
    if (e.key === 'Enter' && isValidCid && !isLoading) {
      handleSynonymSearch();
    }
  };

  // Filter synonyms based on search query
  const filteredSynonyms = synonyms.filter(synonym =>
    synonym.toLowerCase().includes(searchQuery.toLowerCase())
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
          <h3 style={{ marginBottom: '15px', fontWeight: '700' }}>Compound Synonyms</h3>
          <p style={{ marginBottom: '20px', color: 'var(--color-text-secondary)' }}>
            Retrieve all available synonyms for a compound by CID.
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
                  onClick={handleSynonymSearch}
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
                  {isLoading ? 'Fetching...' : 'Get Synonyms'}
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
          
          {(synonyms.length > 0 || error) && (
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

      {synonyms.length > 0 && (
        <div className="utility-results-section" style={{ marginTop: '20px' }}>
          <GlassyContainer>
            <div style={{ display: 'flex', justifyContent: 'space-between', alignItems: 'center', marginBottom: '20px' }}>
              <h4 style={{ color: 'var(--color-text-primary)', fontWeight: '600', margin: 0 }}>
                Synonyms for CID {cid} ({filteredSynonyms.length} of {synonyms.length} shown)
              </h4>
              {synonyms.length > 5 && (
                <div style={{ display: 'flex', alignItems: 'center', gap: '10px', minWidth: '300px' }}>
                  <input
                    type="text"
                    placeholder="Search synonyms..."
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
            
            {filteredSynonyms.length > 0 ? (
              <>
                <div style={{
                  display: 'grid',
                  gridTemplateColumns: 'repeat(auto-fill, minmax(280px, 1fr))',
                  gap: '12px',
                  maxHeight: '500px',
                  overflowY: 'auto',
                  padding: '15px',
                  backgroundColor: 'var(--color-bg-secondary)',
                  borderRadius: '12px',
                  border: '1px solid var(--c-light-border)'
                }}>
                  {filteredSynonyms.map((synonym, index) => (
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
                      onMouseEnter={(e) => {
                        e.target.style.backgroundColor = 'var(--color-bg-primary)';
                        e.target.style.borderColor = 'var(--color-success, #28a745)';
                        e.target.style.transform = 'translateY(-3px)';
                        e.target.style.boxShadow = '0 6px 20px rgba(0, 0, 0, 0.15)';
                      }}
                      onMouseLeave={(e) => {
                        e.target.style.backgroundColor = 'var(--glassy-color)';
                        e.target.style.borderColor = 'var(--c-light-border)';
                        e.target.style.transform = 'translateY(0)';
                        e.target.style.boxShadow = '0 2px 4px rgba(0, 0, 0, 0.05)';
                      }}
                      onClick={() => {
                        navigator.clipboard.writeText(synonym);
                        // Optional: Add a toast notification here
                      }}
                      title="Click to copy to clipboard"
                    >
                      {synonym}
                    </div>
                  ))}
                </div>
                <div style={{
                  marginTop: '15px',
                  fontSize: '13px',
                  color: 'var(--color-text-secondary)',
                  textAlign: 'center',
                  fontStyle: 'italic',
                  backgroundColor: 'rgba(40, 167, 69, 0.1)',
                  padding: '12px',
                  borderRadius: '8px',
                  border: '1px solid rgba(40, 167, 69, 0.2)'
                }}>
                  💡 Click on any synonym to copy it to clipboard
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
                No synonyms match your search "{searchQuery}"
              </div>
            )}
          </GlassyContainer>
        </div>
      )}

      {!isLoading && synonyms.length === 0 && cid && !error && (
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
              No synonyms found for the specified compound.
            </div>
          </GlassyContainer>
        </div>
      )}
    </motion.div>
  );
};

export default CompoundSynonyms;
