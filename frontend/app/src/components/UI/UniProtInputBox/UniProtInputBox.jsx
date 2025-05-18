import React, { useState, useCallback } from 'react';
import axios from 'axios';
import './UniProtInputBox.css'; // We'll create this CSS file next

const UniProtInputBox = ({
    value,
    onChange,
    onValidationComplete,
    header,
    placeholder,
}) => {
    const [isLoading, setIsLoading] = useState(false);
    const [localError, setLocalError] = useState('');
    const [isCurrentValueValidatedSuccessfully, setIsCurrentValueValidatedSuccessfully] = useState(false);

    const handleInputChange = (event) => {
        onChange(event.target.value);
        setLocalError(''); // Clear local error on input change
        setIsCurrentValueValidatedSuccessfully(false); // Reset validation status on change
        // Optionally, clear parent validation status as well if desired
        // onValidationComplete(false, null, ''); 
    };

    const validateAndFetch = useCallback(async () => {
        if (!value || value.trim() === "") {
            const msg = "Please enter a UniProt Accession Key.";
            setLocalError(msg);
            onValidationComplete(false, null, msg);
            setIsCurrentValueValidatedSuccessfully(false);
            return;
        }

        setIsLoading(true);
        setLocalError('');

        try {
            const response = await axios.get(`https://alphafold.ebi.ac.uk/api/prediction/${value.trim()}`);
            if (response.data && response.data.length > 0) {
                onValidationComplete(true, response.data[0], null);
                setIsCurrentValueValidatedSuccessfully(true); // Mark as successfully validated
            } else {
                const errorMessage = "No data found for this accession key.";
                setLocalError(errorMessage);
                onValidationComplete(false, null, errorMessage);
                setIsCurrentValueValidatedSuccessfully(false);
            }
        } catch (err) {
            let errorMessage;
            if (err.response && err.response.status === 400) {
                errorMessage = "Invalid UniProt Accession Key. Please check the key and try again.";
            } else if (err.response) {
                errorMessage = `Error: ${err.response.status} - ${err.response.statusText}`;
            } else {
                errorMessage = `Failed to fetch data: ${err.message}`;
            }
            setLocalError(errorMessage);
            onValidationComplete(false, null, errorMessage);
            setIsCurrentValueValidatedSuccessfully(false);
        } finally {
            setIsLoading(false);
        }
    }, [value, onValidationComplete, onChange]);

    return (
        <div className="uniprot-input-box-container"> {/* Assuming glassy-feel for similar look */}
            {header && <h3 className="uniprot-input-box-header">{header}</h3>}
            <div className="uniprot-input-box-input-wrapper">
                <input
                    type="text"
                    value={value}
                    onChange={handleInputChange}
                    placeholder={placeholder}
                    className="uniprot-input-box-input"
                    disabled={isLoading}
                />
                <button 
                    onClick={validateAndFetch} 
                    disabled={isLoading || !value.trim() || isCurrentValueValidatedSuccessfully}
                    className="uniprot-input-box-button"
                >
                    {isLoading ? 'Validating...' : 'Validate Key'}
                </button>
            </div>
            {localError && <p className="uniprot-input-box-error">{localError}</p>}
        </div>
    );
};

export default UniProtInputBox;
