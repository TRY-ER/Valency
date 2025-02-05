import React, { useEffect } from "react";
import { useState } from "react";
import "./CustomDropdown.css";


function CustomDropdown({ options,
    configState,
    configNames,
    setConfig,
    isDark }) {

    const [isOpen, setIsOpen] = useState(false);

    const toggleDropdown = () => {
        setIsOpen(!isOpen);
    };

    const selectOption = (option) => {

        setConfig((prevConfig) => ({
            ...prevConfig,
            // add all config names to the new config object
            ...configNames.reduce((acc, name) => {
                acc[name] = option[name];
                return acc;
            }, {})
        }));

        setIsOpen(false);
    };

    useEffect(() => {
        setIsOpen(false);
    }, [configState])

    return (
        <div className={`custom-drop-tog-container ${isDark ? "dark" : ""}`}>
            <div className={`custom-dropdown ${isDark ? "dark" : ""}`}>
                <div className={`selected-option ${isDark ? "dark" : ""}`} onClick={toggleDropdown}>
                    {configState[configNames[0]] || "Select an option"}
                </div>
                {isOpen ? (
                    <ul className={`dropdown-options scrollbar-custom ${isDark ? "dark" : ""}`}>
                        {options.map((option, index) => (
                            <li key={index} onClick={() => selectOption(option)}>
                                {option[configNames[0]]}
                            </li>
                        ))}
                    </ul>
                ) : null}
            </div>
            {
                isOpen ?
                    <svg
                        fill="#000000"
                        height="10px"
                        width="10px"
                        viewBox="0 0 330 330"
                        onClick={toggleDropdown}
                    >
                        <path

                            id="XMLID_224_"

                            d="M325.606,229.393l-150.004-150C172.79,76.58,168.974,75,164.996,75c-3.979,0-7.794,1.581-10.607,4.394l-149.996,150c-5.858,5.858-5.858,15.355,0,21.213c5.857,5.857,15.355,5.858,21.213,0l139.39-139.393l139.397,139.393C307.322,253.536,311.161,255,315,255c3.839,0,7.678-1.464,10.607-4.394C331.464,244.748,331.464,235.251,325.606,229.393z" />
                    </svg>
                    :
                    <svg
                        fill="#000000"
                        height="10px"
                        width="10px"
                        viewBox="0 0 330 330"
                        onClick={toggleDropdown}
                    >
                        <path
                            d="M325.607,79.393c-5.857-5.857-15.355-5.858-21.213,0.001l-139.39,139.393L25.607,79.393c-5.857-5.857-15.355-5.858-21.213,0.001c-5.858,5.858-5.858,15.355,0,21.213l150.004,150c2.813,2.813,6.628,4.393,10.606,4.393s7.794-1.581,10.606-4.394l149.996-150C331.465,94.749,331.465,85.251,325.607,79.393z" />

                    </svg>
            }

        </div>
    );
}

export default CustomDropdown;