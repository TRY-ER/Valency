import React, { useEffect, useState } from 'react';
import GlassyContainer from '../../glassy_container/gc';
import InterButtons from '../../buttons/InterButtons';
import { FaCheckDouble } from "react-icons/fa6";
import { IoWarningOutline } from "react-icons/io5";
import { endpoints } from '../../../endpoints/endpoints';
import { call_endpoint_async } from '../../../endpoints/caller';

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
}) {
    const [isValidMol, setIsValidMol] = useState(false)
    const labels = {
        "enabled": "Start Process",
        "disabled": "Select a Positive Number to begin",
        "cancel": "Reset Process"
    }

    const actors = {
        "enabled": onStartProcess,
        "disabled": null,
        "cancel": null
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
                            />
                            {isValidMol ?
                                <FaCheckDouble className="input-status success" />
                                :
                                <IoWarningOutline className={`input-status ${activeMol.length > 0 ? "fail" : "empty"}`} />
                            }
                        </div>
                        < br />
                        <p className="input-header">{numheader}</p>
                        <div className="input-section">
                            <input className="input-box"
                                placeholder={numplaceholder}
                                key={"number"}
                                value={genNum}
                                onChange={(e) => setGenNum(e.target.value)}
                                inputMode='numeric'
                            />
                        </div>
                    </div>
                    <br />
                    {status === "init" ?
                        <InterButtons
                            labels={labels}
                            actors={actors}
                            status={genNum > 0 ? 'enabled' : 'disabled'}
                        />
                        : ""}
                </div>
            </GlassyContainer>
        </div>
    );
}