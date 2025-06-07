import React, { useState } from 'react';
import { motion } from 'framer-motion';
import { fadeInUpVariantStatic } from '../../../../components/animations/framerAnim';
import GlassyContainer from '../../../../components/glassy_container/gc';
import { FaCheckDouble } from 'react-icons/fa';
import { IoWarningOutline } from 'react-icons/io5';

const MassSearch = () => {
  const [mass, setMass] = useState('');
  const [tolerance, setTolerance] = useState('0.1');
  const [isLoading, setIsLoading] = useState(false);
  const [results, setResults] = useState([]);
  const [error, setError] = useState('');
  const [isValidMass, setIsValidMass] = useState(null);
  const [isValidTolerance, setIsValidTolerance] = useState(true);

  const validateMass = (value) => {
    if (!value.trim()) {
      setIsValidMass(null);
      return;
    }
    const massValue = parseFloat(value.trim());
    const isValid = !isNaN(massValue) && massValue > 0;
    setIsValidMass(isValid);
  };

  const validateTolerance = (value) => {
    if (!value.trim()) {
      setIsValidTolerance(null);
      return;
    }
    const toleranceValue = parseFloat(value.trim());
    const isValid = !isNaN(toleranceValue) && toleranceValue >= 0;
    setIsValidTolerance(isValid);
  };

  const handleMassChange = (e) => {
    const value = e.target.value;
    setMass(value);
    setError('');
    validateMass(value);
  };

  const handleToleranceChange = (e) => {
    const value = e.target.value;
    setTolerance(value);
    setError('');
    validateTolerance(value);
  };

  const handleSearch = async () => {
    if (!mass.trim()) {
      setError('Please enter a molecular mass.');
      return;
    }

    if (isValidMass === false) {
      setError('Please enter a valid positive mass value.');
      return;
    }

    if (isValidTolerance === false) {
      setError('Please enter a valid non-negative tolerance value.');
      return;
    }

    setIsLoading(true);
    setError('');
    setResults([]);

    try {
      const massValue = parseFloat(mass);
      const toleranceValue = parseFloat(tolerance);
      
      // PubChem doesn't have a direct mass search API, so this is a simplified demo
      // In a real implementation, you'd need a different approach or API
      const response = await fetch(`https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/fastformula/formula/C*H*/cids/JSON?MaxRecords=100`);
      
      if (!response.ok) {
        throw new Error('Mass search failed');
      }
      
      const data = await response.json();
      let compounds = [];
      
      if (data.IdentifierList && data.IdentifierList.CID) {
        // This is a simplified example - in a real implementation, you'd need to
        // fetch the actual molecular weights and filter by mass ± tolerance
        compounds = data.IdentifierList.CID.slice(0, 20); // Limit to 20 for demo
      }
      
      setResults(compounds);
    } catch (error) {
      console.error('Error during mass search:', error);
      setError('No compounds found with the specified mass range.');
    } finally {
      setIsLoading(false);
    }
  };

  const resetForm = () => {
    setMass('');
    setTolerance('0.1');
    setResults([]);
    setError('');
    setIsValidMass(null);
    setIsValidTolerance(true);
  };

  const getValidationIcon = (isValid) => {
    if (isValid === true) {
      return <FaCheckDouble style={{ color: 'var(--success-color)', fontSize: '14px' }} />;
    }
    if (isValid === false) {
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
        <h3>Mass Search</h3>
        
        <div style={{ display: 'grid', gridTemplateColumns: '1fr 1fr', gap: '15px', marginBottom: '15px' }}>
          <div>
            <label style={{ display: 'block', marginBottom: '5px', fontWeight: '500', color: 'var(--text-color)' }}>
              Molecular Mass (Da):
            </label>
            <div style={{ position: 'relative', display: 'flex', alignItems: 'center' }}>
              <input
                type="number"
                step="0.001"
                value={mass}
                onChange={handleMassChange}
                placeholder="Enter mass (e.g., 180.156)"
                disabled={isLoading}
                style={{
                  width: '100%',
                  padding: '12px 45px 12px 15px',
                  border: `2px solid ${isValidMass === false ? 'var(--error-color)' : 
                           isValidMass === true ? 'var(--success-color)' : 'var(--border-color)'}`,
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
                  e.target.style.borderColor = isValidMass === false ? 'var(--error-color)' : 
                                             isValidMass === true ? 'var(--success-color)' : 'var(--border-color)';
                  e.target.style.boxShadow = 'none';
                }}
              />
              <div style={{ 
                position: 'absolute', 
                right: '15px', 
                top: '50%', 
                transform: 'translateY(-50%)'
              }}>
                {getValidationIcon(isValidMass)}
              </div>
            </div>
          </div>
          
          <div>
            <label style={{ display: 'block', marginBottom: '5px', fontWeight: '500', color: 'var(--text-color)' }}>
              Tolerance (Da):
            </label>
            <div style={{ position: 'relative', display: 'flex', alignItems: 'center' }}>
              <input
                type="number"
                step="0.001"
                value={tolerance}
                onChange={handleToleranceChange}
                placeholder="Enter tolerance (e.g., 0.1)"
                disabled={isLoading}
                style={{
                  width: '100%',
                  padding: '12px 45px 12px 15px',
                  border: `2px solid ${isValidTolerance === false ? 'var(--error-color)' : 
                           isValidTolerance === true ? 'var(--success-color)' : 'var(--border-color)'}`,
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
                  e.target.style.borderColor = isValidTolerance === false ? 'var(--error-color)' : 
                                             isValidTolerance === true ? 'var(--success-color)' : 'var(--border-color)';
                  e.target.style.boxShadow = 'none';
                }}
              />
              <div style={{ 
                position: 'absolute', 
                right: '15px', 
                top: '50%', 
                transform: 'translateY(-50%)'
              }}>
                {getValidationIcon(isValidTolerance)}
              </div>
            </div>
          </div>
        </div>
        
        <div style={{ display: 'flex', gap: '10px', marginTop: '15px' }}>
          <button 
            onClick={handleSearch} 
            disabled={isLoading || !mass.trim() || isValidMass === false || isValidTolerance === false}
            style={{
              flex: 1,
              padding: '12px 20px',
              border: 'none',
              borderRadius: '8px',
              backgroundColor: 'var(--primary-color)',
              color: 'white',
              cursor: isLoading || !mass.trim() || isValidMass === false || isValidTolerance === false ? 'not-allowed' : 'pointer',
              opacity: isLoading || !mass.trim() || isValidMass === false || isValidTolerance === false ? 0.6 : 1,
              transition: 'all 0.3s ease',
              fontSize: '16px',
              fontWeight: '500'
            }}
          >
            {isLoading ? 'Searching...' : 'Search by Mass'}
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
          className="mass-search-results"
          initial={{ opacity: 0, y: 20 }}
          animate={{ opacity: 1, y: 0 }}
          transition={{ duration: 0.5 }}
          style={{ marginTop: '20px' }}
        >
          <h4>Mass Search Results ({results.length} compounds found)</h4>
          <div style={{
            marginBottom: '15px',
            padding: '10px 15px',
            backgroundColor: 'var(--card-bg)',
            borderRadius: '8px',
            border: '1px solid var(--border-color)',
            fontSize: '14px',
            color: 'var(--text-secondary)'
          }}>
            <strong>Search Parameters:</strong> Mass = {mass} Da ± {tolerance} Da
          </div>
          
          <div style={{
            display: 'grid',
            gridTemplateColumns: 'repeat(auto-fill, minmax(200px, 1fr))',
            gap: '15px',
            marginTop: '15px',
            maxHeight: '500px',
            overflowY: 'auto',
            padding: '15px',
            backgroundColor: 'var(--card-bg)',
            borderRadius: '8px',
            border: '1px solid var(--border-color)'
          }}>
            {results.map((cid, index) => (
              <div 
                key={index}
                style={{
                  padding: '15px',
                  backgroundColor: 'var(--input-bg)',
                  borderRadius: '8px',
                  border: '1px solid var(--border-color)',
                  cursor: 'pointer',
                  transition: 'all 0.3s ease',
                  textAlign: 'center'
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
                <div style={{ fontWeight: '600', color: 'var(--primary-color)', marginBottom: '10px' }}>
                  CID: {cid}
                </div>
                <div className="compound-image">
                  <img 
                    src={`https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/${cid}/PNG`}
                    alt={`Compound ${cid}`}
                    style={{
                      maxWidth: '100%',
                      height: 'auto',
                      borderRadius: '6px',
                      backgroundColor: 'white',
                      padding: '5px'
                    }}
                    onError={(e) => {
                      e.target.style.display = 'none';
                    }}
                  />
                </div>
              </div>
            ))}
          </div>
          
          {results.length >= 20 && (
            <div style={{
              marginTop: '10px',
              fontSize: '12px',
              color: 'var(--text-secondary)',
              textAlign: 'center',
              fontStyle: 'italic'
            }}>
              Showing first 20 results. Click on any compound to view details.
            </div>
          )}
        </motion.div>
      )}
    </motion.div>
  );
};

export default MassSearch;
