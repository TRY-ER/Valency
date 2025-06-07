import React, { useState } from 'react';
import { motion } from 'framer-motion';
import { fadeInUpVariantStatic } from '../../../../components/animations/framerAnim';
import GlassyContainer from '../../../../components/glassy_container/gc';
import { FaCheckDouble } from 'react-icons/fa';
import { IoWarningOutline } from 'react-icons/io5';
import { getCidsByXref } from '../../../../services/api/mcpToolsService';

const XrefSearch = () => {
  const [xrefType, setXrefType] = useState('RegistryID');
  const [xrefValue, setXrefValue] = useState('');
  const [isLoading, setIsLoading] = useState(false);
  const [results, setResults] = useState([]);
  const [error, setError] = useState('');
  const [isValidInput, setIsValidInput] = useState(null);

  const xrefTypes = [
    { value: 'RegistryID', label: 'Registry ID' },
    { value: 'RN', label: 'Registry Number' },
    { value: 'PubMedID', label: 'PubMed ID' },
    { value: 'MMDBID', label: 'MMDB ID' },
    { value: 'ProteinGI', label: 'Protein GI' },
    { value: 'NucleotideGI', label: 'Nucleotide GI' },
    { value: 'TaxonomyID', label: 'Taxonomy ID' },
    { value: 'MIMID', label: 'MIM ID' },
    { value: 'GeneID', label: 'Gene ID' },
    { value: 'ProbeID', label: 'Probe ID' },
    { value: 'PatentID', label: 'Patent ID' }
  ];

  const validateInput = (value, source) => {
    if (!value.trim()) {
      setIsValidInput(null);
      return;
    }

    // Basic validation based on source type
    let isValid = false;
    const trimmedValue = value.trim();

    switch (source) {
      case 'name':
        isValid = trimmedValue.length >= 2; // At least 2 characters for compound name
        break;
      case 'smiles':
        // Basic SMILES validation - contains valid characters
        isValid = /^[A-Za-z0-9@+\-\[\]()=#$.:\/\\%]+$/.test(trimmedValue);
        break;
      case 'inchi':
        isValid = trimmedValue.startsWith('InChI=');
        break;
      case 'inchikey':
        isValid = /^[A-Z]{14}-[A-Z]{10}-[A-Z]$/.test(trimmedValue);
        break;
      case 'formula':
        isValid = /^[A-Za-z0-9]+$/.test(trimmedValue);
        break;
      default:
        isValid = trimmedValue.length > 0;
    }

    setIsValidInput(isValid);
  };

  const handleInputChange = (e) => {
    const value = e.target.value;
    setXrefValue(value);
    setError('');
    validateInput(value, xrefType);
  };

  const handleSearch = async () => {
    if (!xrefValue.trim()) {
      setError('Please enter a cross-reference value');
      return;
    }

    setIsLoading(true);
    setError('');
    setResults([]);

    try {
      const response = await getCidsByXref(xrefType, xrefValue.trim());
      if (response && response.cids) {
        setResults(response.cids);
      } else {
        setResults([]);
      }
    } catch (err) {
      setError('Error performing cross-reference search: ' + err.message);
    } finally {
      setIsLoading(false);
    }
  };

  const handleKeyPress = (e) => {
    if (e.key === 'Enter') {
      handleSearch();
    }
  };

  const handleXrefSearch = async () => {
    if (!xrefValue.trim()) {
      setError('Please enter a search term.');
      return;
    }

    if (isValidInput === false) {
      setError('Please enter a valid search term for the selected source.');
      return;
    }

    setIsLoading(true);
    setError('');
    setResults([]);

    try {
      const response = await fetch(`https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/${xrefType}/${encodeURIComponent(xrefValue)}/cid/JSON`);
      if (!response.ok) {
        throw new Error('Cross-reference search failed');
      }
      const data = await response.json();
      
      if (data.IdentifierList && data.IdentifierList.CID) {
        setResults(data.IdentifierList.CID.map(cid => ({ type: 'CID', value: cid })));
      } else if (data.InformationList && data.InformationList.Information) {
        const names = data.InformationList.Information.flatMap(info => 
          info.Synonym ? info.Synonym.slice(0, 5) : []
        );
        setResults(names.map(name => ({ type: 'Name', value: name })));
      } else if (data.PropertyTable && data.PropertyTable.Properties) {
        setResults(data.PropertyTable.Properties.map(prop => ({ type: 'SMILES', value: prop.CanonicalSMILES })));
      } else {
        setResults([]);
      }
    } catch (error) {
      console.error('Error during cross-reference search:', error);
      setError('No results found for the given search term.');
    } finally {
      setIsLoading(false);
    }
  };

  const resetForm = () => {
    setXrefValue('');
    setXrefType('RegistryID');
    setResults([]);
    setError('');
    setIsValidInput(null);
  };

  const getValidationIcon = () => {
    if (isValidInput === true) {
      return <FaCheckDouble style={{ color: 'var(--success-color)', fontSize: '14px' }} />;
    }
    if (isValidInput === false) {
      return <IoWarningOutline style={{ color: 'var(--error-color)', fontSize: '16px' }} />;
    }
    return null;
  };

  const getInputPlaceholder = () => {
    switch (xrefType) {
      case 'name': return 'Enter compound name (e.g., aspirin)';
      case 'smiles': return 'Enter SMILES string (e.g., CC(=O)OC1=CC=CC=C1C(=O)O)';
      case 'inchi': return 'Enter InChI string';
      case 'inchikey': return 'Enter InChI Key';
      case 'formula': return 'Enter molecular formula (e.g., C9H8O4)';
      default: return 'Enter search term';
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
        <h3>Cross-Reference Search</h3>
        
        <div style={{ display: 'grid', gridTemplateColumns: '1fr 1fr', gap: '15px', marginBottom: '15px' }}>
          <div>
            <label style={{ display: 'block', marginBottom: '5px', fontWeight: '500', color: 'var(--text-color)' }}>
              Source Database:
            </label>
            <select
              value={xrefType}
              onChange={(e) => setXrefType(e.target.value)}
              style={{
                width: '100%',
                padding: '12px',
                border: '2px solid var(--border-color)',
                borderRadius: '8px',
                fontSize: '16px',
                backgroundColor: 'var(--input-bg)',
                color: 'var(--text-color)',
                cursor: 'pointer',
                outline: 'none',
                transition: 'all 0.3s ease'
              }}
              onFocus={(e) => {
                e.target.style.borderColor = 'var(--primary-color)';
                e.target.style.boxShadow = '0 0 0 3px rgba(var(--primary-color-rgb), 0.1)';
              }}
              onBlur={(e) => {
                e.target.style.borderColor = 'var(--border-color)';
                e.target.style.boxShadow = 'none';
              }}
            >
              {xrefTypes.map((type) => (
                <option key={type.value} value={type.value}>
                  {type.label}
                </option>
              ))}
            </select>
          </div>

          <div>
            <label style={{ display: 'block', marginBottom: '5px', fontWeight: '500', color: 'var(--text-color)' }}>
              Target Database:
            </label>
            <select
              value="cid"
              onChange={(e) => {}}
              style={{
                width: '100%',
                padding: '12px',
                border: '2px solid var(--border-color)',
                borderRadius: '8px',
                fontSize: '16px',
                backgroundColor: 'var(--input-bg)',
                color: 'var(--text-color)',
                cursor: 'pointer',
                outline: 'none',
                transition: 'all 0.3s ease'
              }}
              onFocus={(e) => {
                e.target.style.borderColor = 'var(--primary-color)';
                e.target.style.boxShadow = '0 0 0 3px rgba(var(--primary-color-rgb), 0.1)';
              }}
              onBlur={(e) => {
                e.target.style.borderColor = 'var(--border-color)';
                e.target.style.boxShadow = 'none';
              }}
            >
              <option value="cid">Compound ID (CID)</option>
              <option value="name">Compound Names</option>
              <option value="smiles">SMILES</option>
            </select>
          </div>
        </div>

        <div style={{ position: 'relative', display: 'flex', alignItems: 'center', gap: '10px' }}>
          <input
            type="text"
            value={xrefValue}
            onChange={handleInputChange}
            placeholder={getInputPlaceholder()}
            disabled={isLoading}
            style={{
              flex: 1,
              padding: '12px 45px 12px 15px',
              border: `2px solid ${isValidInput === false ? 'var(--error-color)' : 
                       isValidInput === true ? 'var(--success-color)' : 'var(--border-color)'}`,
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
              e.target.style.borderColor = isValidInput === false ? 'var(--error-color)' : 
                                         isValidInput === true ? 'var(--success-color)' : 'var(--border-color)';
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
            onClick={handleXrefSearch} 
            disabled={isLoading || !xrefValue.trim() || isValidInput === false}
            style={{
              flex: 1,
              padding: '12px 20px',
              border: 'none',
              borderRadius: '8px',
              backgroundColor: 'var(--primary-color)',
              color: 'white',
              cursor: isLoading || !xrefValue.trim() || isValidInput === false ? 'not-allowed' : 'pointer',
              opacity: isLoading || !xrefValue.trim() || isValidInput === false ? 0.6 : 1,
              transition: 'all 0.3s ease',
              fontSize: '16px',
              fontWeight: '500'
            }}
          >
            {isLoading ? 'Searching...' : 'Search Cross-References'}
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

      {results.length > 0 && (
        <motion.div 
          className="xref-results"
          initial={{ opacity: 0, y: 20 }}
          animate={{ opacity: 1, y: 0 }}
          transition={{ duration: 0.5 }}
          style={{ marginTop: '20px' }}
        >
          <h4>Cross-Reference Results ({results.length} found)</h4>
          <div style={{
            display: 'grid',
            gridTemplateColumns: 'repeat(auto-fill, minmax(300px, 1fr))',
            gap: '15px',
            marginTop: '15px',
            maxHeight: '400px',
            overflowY: 'auto',
            padding: '15px',
            backgroundColor: 'var(--card-bg)',
            borderRadius: '8px',
            border: '1px solid var(--border-color)'
          }}>
            {results.map((result, index) => (
              <div 
                key={index}
                style={{
                  padding: '15px',
                  backgroundColor: 'var(--input-bg)',
                  borderRadius: '8px',
                  border: '1px solid var(--border-color)',
                  cursor: 'pointer',
                  transition: 'all 0.3s ease'
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
                  if (result.type === 'CID') {
                    window.open(`https://pubchem.ncbi.nlm.nih.gov/compound/${result.value}`, '_blank');
                  }
                }}
              >
                <div style={{ fontWeight: '600', color: 'var(--primary-color)', marginBottom: '8px' }}>
                  {result.type}
                </div>
                <div style={{ 
                  fontSize: '14px', 
                  color: 'var(--text-color)',
                  wordBreak: 'break-all',
                  lineHeight: '1.4'
                }}>
                  {result.value}
                </div>
              </div>
            ))}
          </div>
        </motion.div>
      )}
    </motion.div>
  );
};

export default XrefSearch;
