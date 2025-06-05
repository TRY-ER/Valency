import React, { useState, useRef, useEffect } from 'react';
import './CustomSelectDropdown.css';

const CustomSelectDropdown = ({ 
    value, 
    onChange, 
    options, 
    className = '', 
    placeholder = 'Select an option...',
    disabled = false 
}) => {
    const [isOpen, setIsOpen] = useState(false);
    const dropdownRef = useRef(null);

    // Close dropdown when clicking outside
    useEffect(() => {
        const handleClickOutside = (event) => {
            if (dropdownRef.current && !dropdownRef.current.contains(event.target)) {
                setIsOpen(false);
            }
        };

        document.addEventListener('mousedown', handleClickOutside);
        return () => document.removeEventListener('mousedown', handleClickOutside);
    }, []);

    // Close dropdown on escape key
    useEffect(() => {
        const handleKeyDown = (event) => {
            if (event.key === 'Escape') {
                setIsOpen(false);
            }
        };

        if (isOpen) {
            document.addEventListener('keydown', handleKeyDown);
            return () => document.removeEventListener('keydown', handleKeyDown);
        }
    }, [isOpen]);

    const handleOptionClick = (optionValue) => {
        onChange({ target: { value: optionValue } });
        setIsOpen(false);
    };

    const selectedOption = options.find(option => option.value === value);

    return (
        <div 
            className={`custom-select-dropdown ${className} ${disabled ? 'disabled' : ''}`}
            ref={dropdownRef}
        >
            <div 
                className={`select-trigger ${isOpen ? 'open' : ''}`}
                onClick={() => !disabled && setIsOpen(!isOpen)}
                role="button"
                tabIndex={disabled ? -1 : 0}
                onKeyDown={(e) => {
                    if (e.key === 'Enter' || e.key === ' ') {
                        e.preventDefault();
                        !disabled && setIsOpen(!isOpen);
                    }
                }}
            >
                <span className="selected-text">
                    {selectedOption ? selectedOption.label : placeholder}
                </span>
                <span className={`dropdown-arrow ${isOpen ? 'open' : ''}`}>
                    ▼
                </span>
            </div>
            
            {isOpen && !disabled && (
                <div className="select-options">
                    {options.map((option) => (
                        <div
                            key={option.value}
                            className={`select-option ${option.value === value ? 'selected' : ''}`}
                            onClick={() => handleOptionClick(option.value)}
                            role="option"
                            aria-selected={option.value === value}
                        >
                            {option.label}
                        </div>
                    ))}
                </div>
            )}
        </div>
    );
};

export default CustomSelectDropdown;
