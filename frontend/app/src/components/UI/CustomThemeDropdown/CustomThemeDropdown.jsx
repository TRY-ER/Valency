import React, { useState, useEffect, useRef } from 'react';
import { FaChevronDown, FaChevronUp } from 'react-icons/fa';
import './CustomThemeDropdown.css';
import { useTheme } from '../../../contexts/ThemeContext.js'; // Adjusted path

const CustomThemeDropdown = ({ options, selectedValue, onChange }) => {
  const [isOpen, setIsOpen] = useState(false);
  const dropdownRef = useRef(null);
  const { theme } = useTheme(); // To apply dark/light class

  const handleToggle = () => setIsOpen(!isOpen);

  const handleOptionClick = (value) => {
    onChange(value);
    setIsOpen(false);
  };

  useEffect(() => {
    const handleClickOutside = (event) => {
      if (dropdownRef.current && !dropdownRef.current.contains(event.target)) {
        setIsOpen(false);
      }
    };
    document.addEventListener('mousedown', handleClickOutside);
    return () => {
      document.removeEventListener('mousedown', handleClickOutside);
    };
  }, []);

  const selectedOption = options.find(option => option.value === selectedValue);

  return (
    <div className={`custom-theme-dropdown glassy-feel ${theme === 'dark' ? 'dark' : ''}`} ref={dropdownRef}>
      <div className="dropdown-selected" onClick={handleToggle}>
        <span>{selectedOption ? selectedOption.label : 'Select Theme'}</span>
        {isOpen ? <FaChevronUp /> : <FaChevronDown />}
      </div>
      {isOpen && (
        <ul className="dropdown-options">
          {options.map((option) => (
            <li
              key={option.value}
              onClick={() => handleOptionClick(option.value)}
              className={selectedValue === option.value ? 'selected' : ''}
            >
              {option.label}
            </li>
          ))}
        </ul>
      )}
    </div>
  );
};

export default CustomThemeDropdown;
