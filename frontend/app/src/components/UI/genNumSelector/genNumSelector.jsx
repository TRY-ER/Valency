import React, { useEffect, useState } from 'react';
import './genNumSelector.css';
import GlassyContainer from '../../glassy_container/gc';
import InterButtons from '../../buttons/InterButtons';

export default function GenNumSelector({ status, 
                                         onStartProcess,
                                         setGenNum,
                                         genNum,
                                         header,
                                         placeholder }) {
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

    return (
        <div className="gen-num-section">
            <GlassyContainer>
                <div className="gen-num-container">
                    <div className="input-container">
                        <p className="input-header">{header}</p>
                        <div className="input-section">
                            <input className="input-box"
                                placeholder={placeholder}
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