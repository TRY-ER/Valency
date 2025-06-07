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

  const validateSmiles = async (smilesValue) => {
    if (!smilesValue.trim()) {
      setIsValidSmiles(null);
      return;
    }

    setIsValidating(true);
    try {
      const payload = {
        smiles: smilesValue.trim(),
        type: 'smiles'
      };
      const response = await call_endpoint_async(endpoints.validate, payload);
      setIsValidSmiles(response.valid);
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
      const response = await fetch(`https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/substructure/smiles/${encodeURIComponent(smiles)}/cids/JSON`);
      if (!response.ok) {
        throw new Error('Substructure search failed');
      }
      const data = await response.json();
      if (data.IdentifierList && data.IdentifierList.CID) {
        setSubstructureResults(data.IdentifierList.CID.slice(0, 100)); // Limit to first 100 results
      } else {
        setSubstructureResults([]);
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
    clearTimeout(window.smilesValidationTimeout);
  };

  const getValidationIcon = () => {
    if (isValidating) {
      return <div className="validation-spinner"></div>;
    }
    if (isValidSmiles === true) {
      return <FaCheckDouble style={{ color: 'var(--success-color)', fontSize: '14px' }} />;
    }
    if (isValidSmiles === false) {
      return <IoWarningOutline style={{ color: 'var(--error-color)', fontSize: '16px' }} />;
    }
    return null;
  };

  return (
    <motion.div 
      className="utility-component-container"
      variants={fadeInUpVariantStatic}
      initial="hidden"
      animate="visible"
    >
      <GlassyContainer>
        <h3>Substructure Search</h3>
        <div style={{ position: 'relative', display: 'flex', alignItems: 'center', gap: '10px' }}>
          <input
            type="text"
            value={smiles}
            onChange={handleSmilesChange}
            placeholder="Enter SMILES string (e.g., CCO)"
            disabled={isLoading}
            style={{
              flex: 1,
              padding: '12px 45px 12px 15px',
              border: `2px solid ${isValidSmiles === false ? 'var(--error-color)' : 
                       isValidSmiles === true ? 'var(--success-color)' : 'var(--border-color)'}`,
              borderRadius: '8px',
              fontSize: '16px',
              backgroundColor: 'var(--input-bg)',
              color: 'var(--text-color)',
              transition: 'all 0.3s ease',
              outline: 'none'
            }}
            onFocus={(e) => {
              e.target.style.borderColor = 'var(--primary-color)';
              e.target.style.boxShadow = '0 0 0 3px rgba(var(--primary-color-rgb), 0.1)';
            }}
            onBlur={(e) => {
              e.target.style.borderColor = isValidSmiles === false ? 'var(--error-color)' : 
                                         isValidSmiles === true ? 'var(--success-color)' : 'var(--border-color)';
              e.target.style.boxShadow = 'none';
            }}
          />
          <div style={{ 
            position: 'absolute', 
            right: '15px', 
            top: '50%', 
            transform: 'translateY(-50%)',
            display: 'flex',
            alignItems: 'center',
            gap: '10px'
          }}>
            {getValidationIcon()}
          </div>
        </div>
        
        <div style={{ display: 'flex', gap: '10px', marginTop: '15px' }}>
          <button 
            onClick={handleSubstructureSearch} 
            disabled={isLoading || !smiles.trim() || isValidSmiles === false}
            style={{
              flex: 1,
              padding: '12px 20px',
              border: 'none',
              borderRadius: '8px',
              backgroundColor: 'var(--primary-color)',
              color: 'white',
              cursor: isLoading || !smiles.trim() || isValidSmiles === false ? 'not-allowed' : 'pointer',
              opacity: isLoading || !smiles.trim() || isValidSmiles === false ? 0.6 : 1,
              transition: 'all 0.3s ease',
              fontSize: '16px',
              fontWeight: '500'
            }}
          >
            {isLoading ? 'Searching...' : 'Search Substructures'}
          </button>
          
          <button 
            onClick={resetForm}
            style={{
              padding: '12px 20px',
              border: '2px solid var(--border-color)',
              borderRadius: '8px',
              backgroundColor: 'transparent',
              color: 'var(--text-color)',
              cursor: 'pointer',
              transition: 'all 0.3s ease',
              fontSize: '16px',
              fontWeight: '500'
            }}
            onMouseEnter={(e) => {
              e.target.style.backgroundColor = 'var(--hover-bg)';
              e.target.style.borderColor = 'var(--primary-color)';
            }}
            onMouseLeave={(e) => {
              e.target.style.backgroundColor = 'transparent';
              e.target.style.borderColor = 'var(--border-color)';
            }}
          >
            Reset
          </button>
        </div>
      </GlassyContainer>

      {error && (
        <div style={{
          marginTop: '20px',
          padding: '15px',
          backgroundColor: 'rgba(var(--error-color-rgb), 0.1)',
          border: '1px solid var(--error-color)',
          borderRadius: '8px',
          color: 'var(--error-color)',
          fontSize: '14px',
          display: 'flex',
          alignItems: 'center',
          gap: '8px'
        }}>
          <IoWarningOutline />
          {error}
        </div>
      )}

      {substructureResults.length > 0 && (
        <motion.div 
          className="substructure-results"
          initial={{ opacity: 0, y: 20 }}
          animate={{ opacity: 1, y: 0 }}
          transition={{ duration: 0.5 }}
          style={{ marginTop: '20px' }}
        >
          <h4>Substructure Search Results ({substructureResults.length} compounds found)</h4>
          <div style={{
            display: 'grid',
            gridTemplateColumns: 'repeat(auto-fill, minmax(80px, 1fr))',
            gap: '10px',
            marginTop: '15px',
            maxHeight: '400px',
            overflowY: 'auto',
            padding: '10px',
            backgroundColor: 'var(--card-bg)',
            borderRadius: '8px',
            border: '1px solid var(--border-color)'
          }}>
            {substructureResults.map((cid, index) => (
              <div 
                key={index}
                style={{
                  padding: '8px',
                  backgroundColor: 'var(--input-bg)',
                  borderRadius: '6px',
                  textAlign: 'center',
                  border: '1px solid var(--border-color)',
                  cursor: 'pointer',
                  transition: 'all 0.3s ease',
                  fontSize: '14px',
                  fontWeight: '500'
                }}
                onMouseEnter={(e) => {
                  e.target.style.backgroundColor = 'var(--hover-bg)';
                  e.target.style.borderColor = 'var(--primary-color)';
                  e.target.style.transform = 'translateY(-2px)';
                  e.target.style.boxShadow = '0 4px 12px rgba(0, 0, 0, 0.1)';
                }}
                onMouseLeave={(e) => {
                  e.target.style.backgroundColor = 'var(--input-bg)';
                  e.target.style.borderColor = 'var(--border-color)';
                  e.target.style.transform = 'translateY(0)';
                  e.target.style.boxShadow = 'none';
                }}
                onClick={() => window.open(`https://pubchem.ncbi.nlm.nih.gov/compound/${cid}`, '_blank')}
              >
                CID: {cid}
              </div>
            ))}
          </div>
        </motion.div>
      )}
    </motion.div>
  );
};

export default SubstructureSearch;
