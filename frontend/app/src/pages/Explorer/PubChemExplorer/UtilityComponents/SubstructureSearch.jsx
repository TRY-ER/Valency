import React, { useState } from 'react';
import { motion } from 'framer-motion';
import { fadeInUpVariantStatic } from '../../../../components/animations/framerAnim';
import GlassyContainer from '../../../../components/glassy_container/gc';
import { FaCheckDouble } from 'react-icons/fa';
import { IoWarningOutline } from 'react-icons/io5';
import { call_endpoint_async } from '../../../../endpoints/caller';
import { endpoints } from '../../../../endpoints/endpoints';
import { fastSubstructureSearchBySmiles } from '../../../../services/api/mcpToolsService';

const SubstructureSearch = () => {
  const [smiles, setSmiles] = useState('');
  const [substructureResults, setSubstructureResults] = useState([]);
  const [isLoading, setIsLoading] = useState(false);
  const [error, setError] = useState('');
  const [isValidSmiles, setIsValidSmiles] = useState(null);
  const [isValidating, setIsValidating] = useState(false);
  const [searchQuery, setSearchQuery] = useState('');

  const validateSmiles = async (smilesValue) => {
    if (!smilesValue.trim()) {
      setIsValidSmiles(null);
      return;
    }

    setIsValidating(true);
    try {
      const payload = {
        type: 'MOL',
        value: smilesValue.trim()
      };
      const response = await call_endpoint_async(endpoints.validate, payload);
      if (response.data.status === "success") {
        setIsValidSmiles(response.data.valid);
      }
    } catch (error) {
      console.error('SMILES validation error:', error);
      setIsValidSmiles(false);
    } finally {
      setIsValidating(false);
    }
  };

  const handleSmilesChange = (e) => {
    const value = e.target.value;
    setSmiles(value);
    setError('');
    
    // Debounce validation
    clearTimeout(window.smilesValidationTimeout);
    window.smilesValidationTimeout = setTimeout(() => {
      validateSmiles(value);
    }, 500);
  };

  const handleSubstructureSearch = async () => {
    if (!smiles.trim()) {
      setError('Please enter a SMILES string.');
      return;
    }

    if (isValidSmiles === false) {
      setError('Please enter a valid SMILES string.');
      return;
    }

    setIsLoading(true);
    setError('');
    setSubstructureResults([]);

    try {
      const response = await fastSubstructureSearchBySmiles({ smiles: smiles.trim() });
      console.log('Substructure search response:', response); 

      if (response && response.status === "success") {
        let processedResults = response.result["0"];
        processedResults = JSON.parse(processedResults);
         
        // Extract CIDs from the response structure
        if (processedResults && processedResults.data && Array.isArray(processedResults.data)) {
          setSubstructureResults(processedResults.data.slice(0, 100)); // Limit to first 100 results
        } else {
          setSubstructureResults([]);
        }
      } else {
        setError(response.message || 'No compounds found with the given substructure.');
      }
    } catch (error) {
      console.error('Error during substructure search:', error);
      setError('No compounds found with the given substructure.');
    } finally {
      setIsLoading(false);
    }
  };

  const resetForm = () => {
    setSmiles('');
    setSubstructureResults([]);
    setError('');
    setIsValidSmiles(null);
    setIsValidating(false);
    setSearchQuery('');
    clearTimeout(window.smilesValidationTimeout);
  };

  const getValidationIcon = () => {
    if (isValidating) {
      return <div className="validation-spinner"></div>;
    }
    if (isValidSmiles === true) {
      return <FaCheckDouble style={{ fontSize: '2rem', color: 'var(--color-success)' }} />;
    }
    if (isValidSmiles === false) {
      return <IoWarningOutline style={{
        fontSize: '2rem',
        color: smiles.length > 0 ? 'var(--color-alert)' : 'var(--glassy-color)'
      }} />;
    }
    return null;
  };

  // Filter results based on search query
  const filteredResults = substructureResults.filter(cid =>
    cid.toString().toLowerCase().includes(searchQuery.toLowerCase())
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
          <h3 style={{ marginBottom: '15px', fontWeight: '700' }}>Substructure Search</h3>
          <p style={{ marginBottom: '20px', color: 'var(--color-text-secondary)' }}>
            Search for compounds containing a specific substructure using SMILES notation.
          </p>

          {/* Add inline styles for input placeholder */}
          <style>
            {`
              .smiles-input::placeholder {
                color: var(--color-text-secondary);
                font-weight: 400;
              }
              .smiles-input:focus {
                outline: none;
                background-color: var(--glassy-color);
              }
            `}
          </style>

          {/* Custom SMILES Input with Validation */}
          <div style={{ display: 'flex', flexDirection: 'column', gap: '10px' }}>
            <div style={{ display: 'flex', alignItems: 'center', gap: '20px', width: '100%' }}>
              <h4 style={{ fontWeight: '600', color: 'var(--color-text-primary)', margin: 0, minWidth: '180px' }}>
                Enter SMILES:
              </h4>
              <div style={{ display: 'flex', alignItems: 'center', gap: '10px', flex: 1 }}>
                <input
                  type="text"
                  className="smiles-input"
                  value={smiles}
                  onChange={handleSmilesChange}
                  placeholder="Enter SMILES string (e.g., CCO)"
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
                {isValidating ? (
                  <div className="validation-spinner"></div>
                ) : isValidSmiles === true ? (
                  <FaCheckDouble style={{ fontSize: '2rem', color: 'var(--color-success)' }} />
                ) : (
                  <IoWarningOutline style={{
                    fontSize: '2rem',
                    color: smiles.length > 0 ? 'var(--color-alert)' : 'var(--glassy-color)'
                  }} />
                )}
                <button
                  onClick={handleSubstructureSearch}
                  disabled={isLoading || !smiles.trim() || isValidSmiles === false}
                  style={{
                    padding: '12px 20px',
                    borderRadius: '15px',
                    border: 'none',
                    backgroundColor: (isLoading || !smiles.trim() || isValidSmiles === false) ? 'var(--color-disabled, #6c757d)' : 'var(--color-success, #28a745)',
                    color: '#fff',
                    cursor: (isLoading || !smiles.trim() || isValidSmiles === false) ? 'not-allowed' : 'pointer',
                    fontSize: '1em',
                    fontWeight: '500',
                    minHeight: '54px',
                    whiteSpace: 'nowrap'
                  }}
                >
                  {isLoading ? 'Searching...' : 'Search Substructures'}
                </button>
              </div>
            </div>
            {/* Validation message for SMILES */}
            {smiles.length > 0 && isValidSmiles === false && (
              <p style={{
                color: 'var(--color-alert, #ff6b6b)',
                fontSize: '0.9em',
                margin: '5px 0 0 0'
              }}>
                Please enter a valid SMILES string
              </p>
            )}
          </div>
          
          {(substructureResults.length > 0 || error) && (
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

      {substructureResults.length > 0 && (
        <div className="utility-results-section" style={{ marginTop: '20px' }}>
          <GlassyContainer>
            <div style={{ display: 'flex', justifyContent: 'space-between', alignItems: 'center', marginBottom: '20px' }}>
              <h4 style={{ color: 'var(--color-text-primary)', fontWeight: '600', margin: 0 }}>
                Substructure Search Results ({filteredResults.length} of {substructureResults.length} shown)
              </h4>
              {substructureResults.length > 5 && (
                <div style={{ display: 'flex', alignItems: 'center', gap: '10px', minWidth: '300px' }}>
                  <input
                    type="text"
                    placeholder="Search CIDs..."
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
            
            {filteredResults.length > 0 ? (
              <>
                <div style={{
                  display: 'grid',
                  gridTemplateColumns: 'repeat(auto-fill, minmax(120px, 1fr))',
                  gap: '12px',
                  maxHeight: '500px',
                  overflowY: 'auto',
                  padding: '15px',
                  backgroundColor: 'var(--color-bg-secondary)',
                  borderRadius: '12px',
                  border: '1px solid var(--c-light-border)'
                }}>
                  {filteredResults.map((cid, index) => (
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
                        textAlign: 'center',
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
                      onClick={() => window.open(`https://pubchem.ncbi.nlm.nih.gov/compound/${cid}`, '_blank')}
                      title="Click to open in PubChem"
                    >
                      CID: {cid}
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
                  💡 Click on any compound to open it in PubChem
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
                No compounds match your search "{searchQuery}"
              </div>
            )}
          </GlassyContainer>
        </div>
      )}

      {!isLoading && substructureResults.length === 0 && smiles && !error && (
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
              No compounds found with the specified substructure.
            </div>
          </GlassyContainer>
        </div>
      )}
    </motion.div>
  );
};

export default SubstructureSearch;
