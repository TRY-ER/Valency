import React, { useEffect, useState } from 'react';
import GlassyContainer from '../../glassy_container/gc';
import InterButtons from '../../buttons/InterButtons';
import { FaCheckDouble } from "react-icons/fa6";
import { IoWarningOutline } from "react-icons/io5";
import { endpoints } from '../../../endpoints/endpoints';
import { call_endpoint_async } from '../../../endpoints/caller';
import "./numValidSelector.css";
import { FaDownload } from "react-icons/fa";
import { FaRegCopy } from "react-icons/fa";


export default function NumValidSelector({ status,
    onStartProcess,
    setGenNum,
    genNum,
    numheader,
    numplaceholder,
    molheader,
    molplaceholder,
    activeMol,
    setActiveMol,
    inputType,
    handleReset,
    generatedData,
    showNumGen = true,
    labels = {
        "enabled": "Start Process",
        "disabled": "Select Valid Values to Begin",
        "cancel": "Reset Process"
    },
    actors = {
        "enabled": onStartProcess,
        "disabled": null,
        "cancel": handleReset
    }
}) {
    const [isValidMol, setIsValidMol] = useState(false)
    const [strGenData, setStrGenData] = useState("");

    useEffect(() => {
        if (generatedData && generatedData.length > 0) {
            let strData = "";
            generatedData.forEach((item) => {
                strData += `identifier: ${item.identifier}, score: ${item.score}\n`;
            });
            setStrGenData(strData);
        } 
    },[generatedData])

    const handleDownload = () => {
        const element = document.createElement("a");
        const file = new Blob([strGenData], { type: 'text/plain' });
        element.href = URL.createObjectURL(file);
        element.download = "similarity_data.txt";
        document.body.appendChild(element);
        element.click();
    }

    const handleCopy = () => {
        navigator.clipboard.writeText(strGenData);
    }

    useEffect(() => {
        if (activeMol.length > 0) {
            // Check if the molecule is valid
            const payload = {
                type: inputType,
                value: activeMol
            }
            call_endpoint_async(endpoints.validate, payload).then((response) => {
                if (response.data.status === "success") {
                    console.log("active mol", activeMol)
                    console.log("input type >>", inputType);
                    console.log("response >>", response.data);
                    setIsValidMol(response.data.valid);
                }
            }).catch((error) => {
                console.log('error>>', error);
            })
        }
        else {
            setIsValidMol(false);
        }
    }, [activeMol]);

    return (
        <div className="gen-num-section">
            <GlassyContainer>
                <div className="gen-num-container">
                    <div className="input-container">
                        <p className="input-header">{molheader}</p>
                        <div className="input-section">
                            <input className="input-box"
                                placeholder={molplaceholder}
                                key={"mol"}
                                value={activeMol}
                                onChange={(e) => setActiveMol(e.target.value)}
                                disabled={status !== "init"}
                            />
                            {isValidMol ?
                                <FaCheckDouble className="input-status success" />
                                :
                                <IoWarningOutline className={`input-status ${activeMol.length > 0 ? "fail" : "empty"}`} />
                            }
                        </div>
                        < br />
                        {
                            showNumGen &&
                            <>
                                <p className="input-header">{numheader}</p>
                                <div className="input-section">
                                    <input className="input-box"
                                        placeholder={numplaceholder}
                                        key={"number"}
                                        value={genNum}
                                        onChange={(e) => setGenNum(e.target.value)}
                                        inputMode='numeric'
                                        disabled={status !== "init"}
                                    />
                                </div>
                            </>
                        }
                    </div>
                    <br />
                    <div className="gen-num-action-btn-container">
                        <div className="action-btn-left-container">
                            {status === "init" ?
                                <InterButtons
                                    labels={labels}
                                    actors={actors}
                                    status={isValidMol && genNum > 0 ? 'enabled' : 'disabled'}
                                />
                                : ""}
                            {status === "completed" || status === "failed" || status === "no-data" ?
                                <InterButtons
                                    labels={labels}
                                    actors={actors}
                                    status="cancel"
                                />
                                : ""}
                        </div>
                        <div className="action-btn-right-container">
                            {
                                status === "completed" ?
                                    <>
                                        <FaDownload className='action-btn' onClick={handleDownload}/>
                                        <FaRegCopy className='action-btn' onClick={handleCopy}/>
                                    </>
                                    : ""
                            }

                        </div>
                    </div>
                </div>
            </GlassyContainer>
        </div>
    );
}