import React, { useState, useEffect } from "react"; // Removed useContext
import "./InputBox.css";
import { FaCheckDouble } from "react-icons/fa6";
import { IoWarningOutline } from "react-icons/io5";
import GlassyContainer from "../../glassy_container/gc";

import { call_endpoint_async } from "../../../endpoints/caller";
import { endpoints } from "../../../endpoints/endpoints";
// Removed ThemeContext import

const MolInputBox = ({
    activeMol,
    setActiveMol,
    isValidMol,
    setIsValidMol,
    inputType,
    header = "Enter Molecular SMILES Here",
    placeholder = "Try C for demo"
}) => {
    // Removed theme consumption

    useEffect(() => {
        // if (inputType !== "PROT") {
            if (activeMol.length > 0) {
                // Check if the molecule is valid
                const payload = {
                    type: inputType,
                    value: activeMol
                }
                call_endpoint_async(endpoints.validate, payload).then((response) => {
                    if (response.data.status === "success") {
                        console.log("prot val response >>", response.data)
                        setIsValidMol(response.data.valid);
                    }
                }).catch((error) => {
                    console.log('error>>', error);
                })
            }
            else {
                setIsValidMol(false);
            }
        // }
    }, [activeMol]);

    return <>
        <GlassyContainer> {/* Reverted className */}
            <div className="input-container">
                <p className="input-header">{header}</p>
                <div className="input-section">
                    <input className="input-box"
                        placeholder={placeholder}
                        value={activeMol}
                        onChange={(e) => setActiveMol(e.target.value)}
                    />
                    {isValidMol ?
                        <FaCheckDouble className="input-status success" />
                        :
                        <IoWarningOutline className={`input-status ${activeMol.length > 0 ? "fail" : "empty"}`} />
                    }
                </div>
            </div>
        </GlassyContainer>
    </>
}

export default MolInputBox;