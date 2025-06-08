import React from 'react';
import './SimpleInputBox.css';

const SimpleInputBox = ({
    value,
    onChange,
    onSubmit,
    header,
    placeholder,
    buttonText = "Submit",
    isLoading = false,
    disabled = false,
    error = null
}) => {
    const handleInputChange = (event) => {
        onChange(event.target.value);
    };

    const handleSubmit = () => {
        const safeValue = value || '';
        if (onSubmit && safeValue && safeValue.trim() !== "") {
            onSubmit(safeValue.trim());
        }
    };

    const handleKeyPress = (event) => {
        if (event.key === 'Enter') {
            handleSubmit();
        }
    };

    return (
        <div className="simple-input-box-container">
            {header && <h3 className="simple-input-box-header">{header}</h3>}
            <div className="simple-input-box-input-wrapper">
                <input
                    type="text"
                    value={value}
                    onChange={handleInputChange}
                    onKeyPress={handleKeyPress}
                    placeholder={placeholder}
                    className="simple-input-box-input"
                    disabled={isLoading || disabled}
                />
                <button 
                    onClick={handleSubmit} 
                    disabled={isLoading || disabled || !value || !(value && value.trim())}
                    className="simple-input-box-button"
                >
                    {isLoading ? 'Loading...' : buttonText}
                </button>
            </div>
            {error && <p className="simple-input-box-error">{error}</p>}
        </div>
    );
};

export default SimpleInputBox;
