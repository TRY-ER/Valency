import React, { useState } from 'react';
import { motion } from 'framer-motion';
import { fadeInUpVariantStatic } from '../../../../components/animations/framerAnim';
import GlassyContainer from '../../../../components/glassy_container/gc';
import { FaCheckDouble } from 'react-icons/fa';
import { IoWarningOutline } from 'react-icons/io5';
import { getCompoundSynonymsByCid } from '../../../../services/api/mcpToolsService';

const CompoundSynonyms = () => {
  const [cid, setCid] = useState('');
  const [synonyms, setSynonyms] = useState([]);
  const [isLoading, setIsLoading] = useState(false);
  const [error, setError] = useState('');
  const [isValidCid, setIsValidCid] = useState(null);

  const validateCid = (cidValue) => {
    if (!cidValue.trim()) {
      setIsValidCid(null);
      return;
    }

    const trimmedCid = cidValue.trim();
    // CID should be a positive integer
    const isValid = /^\d+$/.test(trimmedCid) && parseInt(trimmedCid) > 0;
    setIsValidCid(isValid);
  };

  const handleCidChange = (e) => {
    const value = e.target.value;
    setCid(value);
    setError('');
    validateCid(value);
  };

  const handleSynonymSearch = async () => {
    if (!cid.trim()) {
      setError('Please enter a Compound ID (CID).');
      return;
    }

    if (isValidCid === false) {
      setError('Please enter a valid CID (positive integer).');
      return;
    }

    setIsLoading(true);
    setError('');
    setSynonyms([]);

    try {
      const response = await fetch(`https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/${cid}/synonyms/JSON`);
      if (!response.ok) {
        throw new Error('Synonym search failed');
      }
      const data = await response.json();
      if (data.InformationList && data.InformationList.Information[0] && data.InformationList.Information[0].Synonym) {
        setSynonyms(data.InformationList.Information[0].Synonym);
      } else {
        setSynonyms([]);
      }
    } catch (error) {
      console.error('Error fetching synonyms:', error);
      setError('No synonyms found for the given CID.');
    } finally {
      setIsLoading(false);
    }
  };

  const resetForm = () => {
    setCid('');
    setSynonyms([]);
    setError('');
    setIsValidCid(null);
  };

  const getValidationIcon = () => {
    if (isValidCid === true) {
      return <FaCheckDouble style={{ color: 'var(--success-color)', fontSize: '14px' }} />;
    }
    if (isValidCid === false) {
      return <IoWarningOutline style={{ color: 'var(--error-color)', fontSize: '16px' }} />;
    }
    return null;
  };

  const handleKeyPress = (e) => {
    if (e.key === 'Enter') {
      handleSynonymSearch();
    }
  };

  return (
    <motion.div 
      className="utility-component-container"
      variants={fadeInUpVariantStatic}
      initial="hidden"
      animate="visible"
    >
      <GlassyContainer>
        <h3>Compound Synonyms</h3>
        <div style={{ position: 'relative', display: 'flex', alignItems: 'center', gap: '10px' }}>
          <input
            type="text"
            value={cid}
            onChange={handleCidChange}
            placeholder="Enter Compound ID (CID), e.g., 2244"
            disabled={isLoading}
            style={{
              flex: 1,
              padding: '12px 45px 12px 15px',
              border: `2px solid ${isValidCid === false ? 'var(--error-color)' : 
                       isValidCid === true ? 'var(--success-color)' : 'var(--border-color)'}`,
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
              e.target.style.borderColor = isValidCid === false ? 'var(--error-color)' : 
                                         isValidCid === true ? 'var(--success-color)' : 'var(--border-color)';
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
            onClick={handleSynonymSearch} 
            disabled={isLoading || !cid.trim() || isValidCid === false}
            style={{
              flex: 1,
              padding: '12px 20px',
              border: 'none',
              borderRadius: '8px',
              backgroundColor: 'var(--primary-color)',
              color: 'white',
              cursor: isLoading || !cid.trim() || isValidCid === false ? 'not-allowed' : 'pointer',
              opacity: isLoading || !cid.trim() || isValidCid === false ? 0.6 : 1,
              transition: 'all 0.3s ease',
              fontSize: '16px',
              fontWeight: '500'
            }}
          >
            {isLoading ? 'Fetching...' : 'Get Synonyms'}
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

      {synonyms.length > 0 && (
        <motion.div 
          className="synonyms-results"
          initial={{ opacity: 0, y: 20 }}
          animate={{ opacity: 1, y: 0 }}
          transition={{ duration: 0.5 }}
          style={{ marginTop: '20px' }}
        >
          <h4>Synonyms for CID {cid} ({synonyms.length} found)</h4>
          <div style={{
            display: 'grid',
            gridTemplateColumns: 'repeat(auto-fill, minmax(250px, 1fr))',
            gap: '12px',
            marginTop: '15px',
            maxHeight: '400px',
            overflowY: 'auto',
            padding: '15px',
            backgroundColor: 'var(--card-bg)',
            borderRadius: '8px',
            border: '1px solid var(--border-color)'
          }}>
            {synonyms.map((synonym, index) => (
              <div 
                key={index}
                style={{
                  padding: '12px 15px',
                  backgroundColor: 'var(--input-bg)',
                  borderRadius: '8px',
                  border: '1px solid var(--border-color)',
                  cursor: 'pointer',
                  transition: 'all 0.3s ease',
                  fontSize: '14px',
                  lineHeight: '1.4',
                  wordBreak: 'break-word'
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
            marginTop: '10px',
            fontSize: '12px',
            color: 'var(--text-secondary)',
            textAlign: 'center',
            fontStyle: 'italic'
          }}>
            Click on any synonym to copy it to clipboard
          </div>
        </motion.div>
      )}
    </motion.div>
  );
};

export default CompoundSynonyms;
