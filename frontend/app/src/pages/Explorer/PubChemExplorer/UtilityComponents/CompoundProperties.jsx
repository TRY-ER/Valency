import React, { useState, useRef, useEffect } from 'react';
import { motion, AnimatePresence } from "framer-motion";
import { fadeInUpVariantStatic } from '../../../../components/animations/framerAnim';
import GlassyContainer from '../../../../components/glassy_container/gc';
import { getCompoundProperties } from '../../../../services/api/mcpToolsService';
import { FaCheckDouble, FaChevronDown, FaChevronUp } from "react-icons/fa";
import { IoWarningOutline } from "react-icons/io5";

const AVAILABLE_PROPERTIES = [
  { key: 'MolecularFormula', label: 'Molecular Formula', description: 'Molecular formula' },
  { key: 'MolecularWeight', label: 'Molecular Weight', description: 'Sum of all atomic weights (g/mol)' },
  { key: 'SMILES', label: 'SMILES', description: 'Simplified Molecular Input Line Entry System' },
  { key: 'InChI', label: 'InChI', description: 'Standard IUPAC International Chemical Identifier' },
  { key: 'InChIKey', label: 'InChI Key', description: 'Hashed version of InChI (27 characters)' },
  { key: 'IUPACName', label: 'IUPAC Name', description: 'Chemical name according to IUPAC nomenclatures' },
  { key: 'Title', label: 'Title', description: 'Compound summary page title' },
  { key: 'XLogP', label: 'XLogP', description: 'Octanol-water partition coefficient' },
  { key: 'ExactMass', label: 'Exact Mass', description: 'Mass of most likely isotopic composition' },
  { key: 'MonoisotopicMass', label: 'Monoisotopic Mass', description: 'Mass using most abundant isotopes' },
  { key: 'TPSA', label: 'TPSA', description: 'Topological polar surface area' },
  { key: 'Complexity', label: 'Complexity', description: 'Molecular complexity rating' },
  { key: 'Charge', label: 'Charge', description: 'Total charge of molecule' },
  { key: 'HBondDonorCount', label: 'H-Bond Donors', description: 'Number of hydrogen-bond donors' },
  { key: 'HBondAcceptorCount', label: 'H-Bond Acceptors', description: 'Number of hydrogen-bond acceptors' },
  { key: 'RotatableBondCount', label: 'Rotatable Bonds', description: 'Number of rotatable bonds' },
  { key: 'HeavyAtomCount', label: 'Heavy Atoms', description: 'Number of non-hydrogen atoms' },
  { key: 'IsotopeAtomCount', label: 'Isotope Atoms', description: 'Number of atoms with enriched isotopes' },
  { key: 'AtomStereoCount', label: 'Atom Stereo Count', description: 'Atoms with tetrahedral stereo' },
  { key: 'DefinedAtomStereoCount', label: 'Defined Atom Stereo', description: 'Atoms with defined tetrahedral stereo' },
  { key: 'UndefinedAtomStereoCount', label: 'Undefined Atom Stereo', description: 'Atoms with undefined tetrahedral stereo' },
  { key: 'BondStereoCount', label: 'Bond Stereo Count', description: 'Bonds with planar stereo' },
  { key: 'DefinedBondStereoCount', label: 'Defined Bond Stereo', description: 'Bonds with defined planar stereo' },
  { key: 'UndefinedBondStereoCount', label: 'Undefined Bond Stereo', description: 'Bonds with undefined planar stereo' },
  { key: 'CovalentUnitCount', label: 'Covalent Units', description: 'Number of covalently bound units' },
  { key: 'Volume3D', label: '3D Volume', description: 'Analytic volume of first conformer' },
  { key: 'FeatureCount3D', label: '3D Features', description: 'Total number of 3D features' },
  { key: 'FeatureAcceptorCount3D', label: '3D Acceptors', description: 'H-bond acceptors in 3D' },
  { key: 'FeatureDonorCount3D', label: '3D Donors', description: 'H-bond donors in 3D' },
  { key: 'ConformerCount3D', label: 'Conformers', description: 'Number of conformers in model' }
];

const CompoundProperties = ({ toolData = null }) => {
  const [cid, setCid] = useState('');
  const [isLoading, setIsLoading] = useState(false);
  const [properties, setProperties] = useState(null);
  const [error, setError] = useState('');
  const [isValidCid, setIsValidCid] = useState(false);
  const [selectedProperties, setSelectedProperties] = useState(['MolecularWeight', 'InChIKey']);
  const [isDropdownOpen, setIsDropdownOpen] = useState(false);
  const dropdownRef = useRef(null);

  // Effect to handle initial data passed as props (similar to PubChemGetter pattern)
  useEffect(() => {
    console.log('tool data received >>', toolData);
  }, [toolData]);

  useEffect(() => {
    console.log("properties updated >>", properties);
  }, [properties]);

  // Effect to handle initial data passed as props
  useEffect(() => {
    if (toolData) {
      console.log('Processing tool data:', toolData);
      
      // toolData should contain the properties directly (no type/content structure)
      // Set properties directly from toolData
      // setProperties(toolData);
      
      // Clear any existing errors
      setError('');
      
      // If toolData has CID information, extract and set it
      if (toolData.CID) {
        setCid(String(toolData.CID));
        setIsValidCid(true);
      }
    }
  }, [toolData]);

  // Close dropdown when clicking outside
  useEffect(() => {
    const handleClickOutside = (event) => {
      if (dropdownRef.current && !dropdownRef.current.contains(event.target)) {
        setIsDropdownOpen(false);
      }
    };

    document.addEventListener('mousedown', handleClickOutside);
    return () => {
      document.removeEventListener('mousedown', handleClickOutside);
    };
  }, []);

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

  const handlePropertyToggle = (propertyKey) => {
    setSelectedProperties(prev => {
      if (prev.includes(propertyKey)) {
        return prev.filter(p => p !== propertyKey);
      } else {
        return [...prev, propertyKey];
      }
    });
  };

  const handleSelectAll = () => {
    if (selectedProperties.length === AVAILABLE_PROPERTIES.length) {
      setSelectedProperties([]);
    } else {
      setSelectedProperties(AVAILABLE_PROPERTIES.map(p => p.key));
    }
  };

  const getSelectedPropertiesText = () => {
    if (selectedProperties.length === 0) return 'No properties selected';
    if (selectedProperties.length === AVAILABLE_PROPERTIES.length) return 'All properties selected';
    if (selectedProperties.length === 1) return '1 property selected';
    return `${selectedProperties.length} properties selected`;
  };

  const handleSearch = async () => {
    if (!cid.trim()) {
      setError('Please enter a CID');
      return;
    }

    if (selectedProperties.length === 0) {
      setError('Please select at least one property to fetch');
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
      console.log('Fetching properties for CID:', cidValue);
      console.log('Selected properties:', selectedProperties);
      const response = await getCompoundProperties({ 
        cids: [String(cidValue)], 
        properties_list: selectedProperties 
      });
      console.log('Response:', response);
      if (response && response.status === "success") {
        var processible = response.result["0"]
        processible = JSON.parse(processible);

        console.log('Processed Properties:', processible);
        if (processible && processible.status === "success") {
          setProperties(processible.result.PropertyTable.Properties[0])
        }
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
    if (typeof value === 'object' && value !== null) {
      // If the value is an object, convert it to a JSON string
      return JSON.stringify(value);
    }
    if (value === null || value === undefined) {
      return 'N/A';
    }
    return String(value);
  };

  const getPropertyLabel = (key) => {
    const property = AVAILABLE_PROPERTIES.find(p => p.key === key);
    return property ? property.label : key.replace(/([A-Z])/g, ' $1').replace(/^./, str => str.toUpperCase());
  };

  const handleReset = () => {
    setProperties(null);
    setCid("");
    setError("");
    setIsValidCid(false);
    setSelectedProperties(['MolecularWeight', 'InChIKey']);
    setIsDropdownOpen(false);
  };

  // Effect to handle initial data passed as props (similar to PubChemGetter pattern)
  useEffect(() => {
    console.log('tool data received >>', toolData);
  }, [toolData]);

  useEffect(() => {
    console.log("properties updated >>", properties);
  }, [properties]);

  // Effect to handle initial data passed as props
  useEffect(() => {
    if (toolData) {
      console.log('Processing tool data:', toolData);
      
      // Try to extract properties from different possible structures
      let extractedProperties = null;
      
      // Case 1: toolData already has the correct flat structure
      if (toolData.MolecularWeight || toolData.MolecularFormula || toolData.SMILES) {
        extractedProperties = toolData;
      }
      // Case 2: toolData has PropertyTable structure
      else if (toolData.PropertyTable && toolData.PropertyTable.Properties && toolData.PropertyTable.Properties[0]) {
        extractedProperties = toolData.PropertyTable.Properties[0];
      }
      // Case 3: toolData has result.PropertyTable structure (similar to API response)
      else if (toolData.result && toolData.result.PropertyTable && toolData.result.PropertyTable.Properties && toolData.result.PropertyTable.Properties[0]) {
        extractedProperties = toolData.result.PropertyTable.Properties[0];
      }
      // Case 4: toolData might have a nested structure like the API response
      else if (toolData.result && toolData.result["0"]) {
        try {
          const parsed = typeof toolData.result["0"] === 'string' ? JSON.parse(toolData.result["0"]) : toolData.result["0"];
          if (parsed.result && parsed.result.PropertyTable && parsed.result.PropertyTable.Properties && parsed.result.PropertyTable.Properties[0]) {
            extractedProperties = parsed.result.PropertyTable.Properties[0];
          }
        } catch (e) {
          console.error('Error parsing nested toolData:', e);
        }
      }
      
      if (extractedProperties) {
        setProperties(extractedProperties);
        // Clear any existing errors
        setError('');
        
        // If toolData has CID information, extract and set it
        if (toolData.CID) {
          setCid(String(toolData.CID));
          setIsValidCid(true);
        } else if (extractedProperties.CID) {
          setCid(String(extractedProperties.CID));
          setIsValidCid(true);
        }
      } else {
        console.warn('Could not extract properties from toolData:', toolData);
        // Set the raw toolData but it might not render properly
        setProperties(toolData);
      }
    }
  }, [toolData]);

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

          {/* Property Selector */}
          <div style={{ marginBottom: '20px', display: 'flex', alignItems: 'center', gap: '20px' }}>
            <h4 style={{ fontWeight: '600', color: 'var(--color-text-primary)', margin: 0, minWidth: '180px' }}>
              Select Properties:
            </h4>
            <div style={{ position: 'relative', maxWidth: '400px', flex: 1 }} ref={dropdownRef}>
              <button
                onClick={() => setIsDropdownOpen(!isDropdownOpen)}
                disabled={isLoading}
                style={{
                  width: '100%',
                  padding: '12px',
                  borderRadius: '15px',
                  border: 'none',
                  backgroundColor: 'var(--glassy-color)',
                  color: 'var(--color-text-primary)',
                  fontSize: '14px',
                  fontWeight: '600',
                  cursor: isLoading ? 'not-allowed' : 'pointer',
                  display: 'flex',
                  justifyContent: 'space-between',
                  alignItems: 'center',
                  minHeight: '48px'
                }}
              >
                <span style={{ fontSize: '14px' }}>{getSelectedPropertiesText()}</span>
                {isDropdownOpen ? <FaChevronUp /> : <FaChevronDown />}
              </button>
              
              <AnimatePresence>
                {isDropdownOpen && (
                <motion.div 
                  style={{
                    position: 'absolute',
                    top: '100%',
                    left: 0,
                    right: 0,
                    zIndex: 1000,
                    backgroundColor: 'var(--color-bg-primary)',
                    border: '1px solid var(--c-light-border)',
                    borderRadius: '15px',
                    marginTop: '5px',
                    maxHeight: '300px',
                    overflowY: 'auto',
                    boxShadow: '0 4px 12px rgba(0, 0, 0, 0.15)'
                  }}
                  initial={{ opacity: 0, y: -10, scale: 0.95 }}
                  animate={{ opacity: 1, y: 0, scale: 1 }}
                  exit={{ opacity: 0, y: -10, scale: 0.95 }}
                  transition={{ duration: 0.2, ease: "easeOut" }}
                >
                  {/* Select All Option */}
                  <motion.div
                    onClick={handleSelectAll}
                    style={{
                      padding: '12px 16px',
                      cursor: 'pointer',
                      borderBottom: '1px solid var(--c-light-border)',
                      backgroundColor: 'var(--color-bg-secondary)',
                      fontWeight: '600',
                      color: 'var(--color-text-primary)',
                      textAlign: 'left',
                      display: 'flex',
                      alignItems: 'center'
                    }}
                    initial={{ opacity: 0, x: -10 }}
                    animate={{ opacity: 1, x: 0 }}
                    transition={{ duration: 0.2, delay: 0.05 }}
                    whileHover={{ 
                      backgroundColor: 'var(--color-bg-secondary)',
                      scale: 1.02,
                      transition: { duration: 0.15 }
                    }}
                  >
                    <input
                      type="checkbox"
                      checked={selectedProperties.length === AVAILABLE_PROPERTIES.length}
                      onChange={() => {}}
                      style={{ 
                        marginRight: '10px',
                        accentColor: 'var(--color-success, #28a745)'
                      }}
                    />
                    {selectedProperties.length === AVAILABLE_PROPERTIES.length ? 'Deselect All' : 'Select All'}
                  </motion.div>
                  
                  {/* Individual Properties */}
                  {AVAILABLE_PROPERTIES.map((property, index) => (
                    <motion.div
                      key={property.key}
                      onClick={() => handlePropertyToggle(property.key)}
                      style={{
                        padding: '10px 16px',
                        cursor: 'pointer',
                        borderBottom: '1px solid var(--c-light-border)',
                        backgroundColor: selectedProperties.includes(property.key) 
                          ? 'rgba(40, 167, 69, 0.1)' 
                          : 'transparent',
                        textAlign: 'left'
                      }}
                      initial={{ opacity: 0, x: -10 }}
                      animate={{ opacity: 1, x: 0 }}
                      transition={{ 
                        duration: 0.2, 
                        delay: 0.1 + (index * 0.02),
                        ease: "easeOut"
                      }}
                      whileHover={{ 
                        backgroundColor: selectedProperties.includes(property.key) 
                          ? 'rgba(40, 167, 69, 0.15)' 
                          : 'rgba(255, 255, 255, 0.05)',
                        scale: 1.01,
                        transition: { duration: 0.15 }
                      }}
                    >
                      <div style={{ display: 'flex', alignItems: 'flex-start', gap: '10px' }}>
                        <input
                          type="checkbox"
                          checked={selectedProperties.includes(property.key)}
                          onChange={() => {}}
                          style={{ 
                            marginTop: '2px',
                            accentColor: 'var(--color-success, #28a745)'
                          }}
                        />
                        <div style={{ flex: 1, textAlign: 'left' }}>
                          <div style={{
                            fontWeight: '500',
                            color: 'var(--color-text-primary)',
                            marginBottom: '2px',
                            textAlign: 'left'
                          }}>
                            {property.label}
                          </div>
                          <div style={{
                            fontSize: '0.85em',
                            color: 'var(--color-text-secondary)',
                            lineHeight: '1.3',
                            textAlign: 'left'
                          }}>
                            {property.description}
                          </div>
                        </div>
                      </div>
                    </motion.div>
                  ))}
                </motion.div>
                )}
              </AnimatePresence>
            </div>
          </div>

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
                  <motion.div
                    initial={{ scale: 0, opacity: 0 }}
                    animate={{ scale: 1, opacity: 1 }}
                    transition={{ duration: 0.3, ease: "easeOut" }}
                  >
                    <FaCheckDouble style={{ fontSize: '2rem', color: 'var(--color-success)' }} />
                  </motion.div>
                ) : (
                  <motion.div
                    initial={{ scale: 0, opacity: 0 }}
                    animate={{ scale: 1, opacity: 1 }}
                    transition={{ duration: 0.3, ease: "easeOut" }}
                  >
                    <IoWarningOutline style={{
                      fontSize: '2rem',
                      color: cid.length > 0 ? 'var(--color-alert)' : 'var(--glassy-color)'
                    }} />
                  </motion.div>
                )}
                <motion.button
                  onClick={handleSearch}
                  disabled={!isValidCid || isLoading || selectedProperties.length === 0}
                  style={{
                    padding: '12px 20px',
                    borderRadius: '15px',
                    border: 'none',
                    backgroundColor: (!isValidCid || isLoading || selectedProperties.length === 0) ? 'var(--color-disabled, #6c757d)' : 'var(--color-success, #28a745)',
                    color: '#fff',
                    cursor: (!isValidCid || isLoading || selectedProperties.length === 0) ? 'not-allowed' : 'pointer',
                    fontSize: '1em',
                    fontWeight: '500',
                    minHeight: '54px',
                    whiteSpace: 'nowrap'
                  }}
                  whileHover={
                    (!isValidCid || isLoading || selectedProperties.length === 0) ? {} : {
                      scale: 1.05,
                      backgroundColor: '#1e7e34',
                      transition: { duration: 0.2 }
                    }
                  }
                  whileTap={
                    (!isValidCid || isLoading || selectedProperties.length === 0) ? {} : {
                      scale: 0.98,
                      transition: { duration: 0.1 }
                    }
                  }
                >
                  {isLoading ? 'Fetching...' : 'Get Properties'}
                </motion.button>
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
            {/* Validation message for properties */}
            {selectedProperties.length === 0 && (
              <p style={{
                color: 'var(--color-alert, #ff6b6b)',
                fontSize: '0.9em',
                margin: '5px 0 0 0'
              }}>
                Please select at least one property to fetch
              </p>
            )}
          </div>
          {(properties || error) && (
            <div style={{ display: 'flex', justifyContent: 'flex-end', marginTop: '20px' }}>
              <motion.button
                onClick={handleReset}
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
                whileHover={{
                  backgroundColor: '#c82333',
                  scale: 1.05,
                  boxShadow: '0 4px 8px rgba(0, 0, 0, 0.15)',
                  transition: { duration: 0.2 }
                }}
                whileTap={{
                  scale: 0.98,
                  transition: { duration: 0.1 }
                }}
                initial={{ opacity: 0, y: 10 }}
                animate={{ opacity: 1, y: 0 }}
                transition={{ duration: 0.3, delay: 0.2 }}
              >
                Reset All
              </motion.button>
            </div>
          )}
        </GlassyContainer>
      </div>

      {error && (
        <motion.div 
          className="utility-results-section" 
          style={{ marginTop: '20px' }}
          initial={{ opacity: 0, y: 30 }}
          animate={{ opacity: 1, y: 0 }}
          exit={{ opacity: 0, y: -30 }}
          transition={{ duration: 0.4, ease: "easeOut" }}
        >
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
        </motion.div>
      )}

      {properties && (
        <motion.div 
          className="utility-results-section" 
          style={{ marginTop: '20px' }}
          initial={{ opacity: 0, y: 50 }}
          animate={{ opacity: 1, y: 0 }}
          exit={{ opacity: 0, y: -30 }}
          transition={{ duration: 0.5, ease: "easeOut", delay: 0.1 }}
        >
          <GlassyContainer>
            <motion.h4 
              style={{ marginBottom: '20px', color: 'var(--color-text-primary)', fontWeight: '600' }}
              initial={{ opacity: 0, y: 20 }}
              animate={{ opacity: 1, y: 0 }}
              transition={{ duration: 0.4, delay: 0.2 }}
            >
              Properties for CID {cid}
            </motion.h4>
            <div className="properties-list" style={{
              display: 'flex',
              flexDirection: 'column',
              gap: '0'
            }}>
              {Object.entries(properties).map(([key, value], index) => (
                <motion.div 
                  key={key} 
                  className="property-item" 
                  style={{
                    display: 'flex',
                    justifyContent: 'space-between',
                    alignItems: 'flex-start',
                    padding: '16px 20px',
                    borderBottom: index < Object.entries(properties).length - 1 ? '1px solid var(--c-light-border)' : 'none',
                    backgroundColor: index % 2 === 0 ? 'transparent' : 'rgba(255, 255, 255, 0.02)',
                    transition: 'background-color 0.2s ease',
                    borderRadius: index === 0 ? '8px 8px 0 0' : index === Object.entries(properties).length - 1 ? '0 0 8px 8px' : '0'
                  }}
                  initial={{ opacity: 0, x: -20 }}
                  animate={{ opacity: 1, x: 0 }}
                  transition={{ 
                    duration: 0.4, 
                    delay: 0.3 + (index * 0.05),
                    ease: "easeOut"
                  }}
                  whileHover={{ 
                    backgroundColor: 'rgba(255, 255, 255, 0.05)',
                    transition: { duration: 0.2 }
                  }}
                >
                  <div className="property-label" style={{
                    fontWeight: '600',
                    color: 'var(--color-text-primary)',
                    fontSize: '0.95em',
                    minWidth: '200px',
                    maxWidth: '300px',
                    paddingRight: '20px'
                  }}>
                    {getPropertyLabel(key)}
                  </div>
                  <div className="property-value" style={{
                    fontWeight: '500',
                    color: 'var(--color-text-secondary)',
                    fontSize: '0.9em',
                    fontFamily: typeof value === 'number' ? 'Monaco, Menlo, monospace' : 'inherit',
                    textAlign: 'right',
                    flex: 1,
                    wordBreak: 'break-word',
                    lineHeight: '1.4'
                  }}>
                    {formatPropertyValue(value)}
                  </div>
                </motion.div>
              ))}
            </div>
          </GlassyContainer>
        </motion.div>
      )}

      {!isLoading && !properties && cid && !error && (
        <motion.div 
          className="utility-results-section" 
          style={{ marginTop: '20px' }}
          initial={{ opacity: 0, y: 30 }}
          animate={{ opacity: 1, y: 0 }}
          exit={{ opacity: 0, y: -30 }}
          transition={{ duration: 0.4, ease: "easeOut" }}
        >
          <GlassyContainer>
            <motion.div 
              style={{
                textAlign: 'center',
                color: 'var(--color-text-secondary)',
                fontStyle: 'italic',
                padding: '40px 20px',
                background: 'var(--color-bg-secondary)',
                border: '2px dashed var(--c-light-border)',
                borderRadius: '12px'
              }}
              initial={{ opacity: 0, scale: 0.95 }}
              animate={{ opacity: 1, scale: 1 }}
              transition={{ duration: 0.3, delay: 0.1 }}
            >
              No properties found for the specified compound.
            </motion.div>
          </GlassyContainer>
        </motion.div>
      )}
    </motion.div>
  );
};

export default CompoundProperties;
