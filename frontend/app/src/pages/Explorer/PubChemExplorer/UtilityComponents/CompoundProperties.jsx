import React, { useState } from 'react';
import { motion } from "framer-motion";
import { fadeInUpVariantStatic } from '../../../../components/animations/framerAnim';
import { getCompoundProperties } from '../../../../services/api/mcpToolsService';
import GlassyContainer from '../../../../components/glassy_container/gc';
import { FaCheckDouble } from "react-icons/fa";
import { IoWarningOutline } from "react-icons/io5";

const CompoundProperties = () => {
  const [cid, setCid] = useState('');
  const [isLoading, setIsLoading] = useState(false);
  const [properties, setProperties] = useState(null);
  const [error, setError] = useState('');
  const [isValidCid, setIsValidCid] = useState(false);

  // Validate CID in real-time
  const validateCid = (value) => {
    const cidValue = parseInt(value.trim());
    return !isNaN(cidValue) && cidValue > 0;
  };

  const handleCidChange = (value) => {
    setCid(value);
    setIsValidCid(validateCid(value));
    if (error) setError('');
  };

  const handleSearch = async () => {
    if (!cid.trim()) {
      setError('Please enter a CID');
      return;
    }

    const cidValue = parseInt(cid.trim());
    if (isNaN(cidValue) || cidValue <= 0) {
      setError('Please enter a valid positive CID number');
      return;
    }

    setIsLoading(true);
    setError('');
    setProperties(null);

    try {
      const response = await getCompoundProperties(cidValue);
      if (response && response.properties) {
        setProperties(response.properties);
      } else {
        setProperties(null);
      }
    } catch (err) {
      setError('Error fetching properties: ' + err.message);
    } finally {
      setIsLoading(false);
    }
  };

  const handleKeyPress = (e) => {
    if (e.key === 'Enter' && isValidCid && !isLoading) {
      handleSearch();
    }
  };

  const formatPropertyValue = (value) => {
    if (typeof value === 'number') {
      return value.toFixed(4);
    }
    return value;
  };

  const handleReset = () => {
    setProperties(null);
    setCid("");
    setError("");
    setIsValidCid(false);
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
          <h3 style={{ marginBottom: '15px', fontWeight: '700' }}>Compound Properties</h3>
          <p style={{ marginBottom: '15px', color: 'var(--color-text-secondary)' }}>
            Retrieve detailed physicochemical properties for a compound by CID.
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
            <div style={{ display: 'flex', alignItems: 'center', gap: '10px', width: '100%' }}>
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
                onClick={handleSearch}
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
                {isLoading ? 'Fetching...' : 'Get Properties'}
              </button>
            </div>
            {/* Validation message */}
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
          {(properties || error) && (
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

      {properties && (
        <div className="utility-results-section" style={{ marginTop: '20px' }}>
          <GlassyContainer>
            <h4 style={{ marginBottom: '15px', color: 'var(--color-text-primary)', fontWeight: '600' }}>
              Properties for CID {cid}
            </h4>
            <div className="compound-preview" style={{ marginBottom: '20px', textAlign: 'center' }}>
              <img 
                src={`https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/${cid}/PNG`}
                alt={`Compound ${cid}`}
                style={{
                  maxWidth: '300px',
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
            <div className="properties-grid" style={{
              display: 'grid',
              gridTemplateColumns: 'repeat(auto-fit, minmax(250px, 1fr))',
              gap: '16px'
            }}>
              {Object.entries(properties).map(([key, value]) => (
                <div key={key} className="property-item" style={{
                  background: 'var(--color-bg-secondary)',
                  border: '1px solid var(--c-light-border)',
                  borderRadius: '8px',
                  padding: '16px'
                }}>
                  <div className="property-label" style={{
                    fontWeight: '500',
                    color: 'var(--color-text-secondary)',
                    marginBottom: '8px',
                    fontSize: '0.9em'
                  }}>
                    {key.replace(/([A-Z])/g, ' $1').replace(/^./, str => str.toUpperCase())}:
                  </div>
                  <div className="property-value" style={{
                    fontWeight: '600',
                    color: 'var(--color-text-primary)',
                    fontSize: '1em',
                    fontFamily: typeof value === 'number' ? 'Monaco, Menlo, monospace' : 'inherit'
                  }}>
                    {formatPropertyValue(value)}
                  </div>
                </div>
              ))}
            </div>
          </GlassyContainer>
        </div>
      )}

      {!isLoading && !properties && cid && !error && (
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
              No properties found for the specified compound.
            </div>
          </GlassyContainer>
        </div>
      )}
    </motion.div>
  );
};

export default CompoundProperties;
