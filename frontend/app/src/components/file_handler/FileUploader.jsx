import React, { useEffect, useState } from 'react';
import './FileUploader.css';
import GlassyContainer from '../glassy_container/gc';
import InterButtons from '../buttons/InterButtons';

export default function FileUploader({ status, setRawFile, onStartProcess }) {
    const [lines, setLines] = useState([]);
    const [fileName, setFileName] = useState('');
    const [totalLines, setTotalLines] = useState(0);
    const [fileObj, setFileObj] = useState(null);

    const handleFileChange = (e) => {
        const file = e.target.files[0];
        if (!file) return;
        setFileObj(file);
        setRawFile(file);
        setFileName(file.name);

        const reader = new FileReader();
        reader.onload = (event) => {
            const textContent = event.target.result;
            const splittedLines = textContent.split('\n');
            setLines(splittedLines.slice(0, 5));
            setTotalLines(splittedLines.length);
        };
        reader.readAsText(file);
    };

    const labels = {
        "enabled": "Start Process",
        "disabled": "Select contents to start",
        "cancel": "Reset Process"
    }

    const actors = {
        "enabled": onStartProcess,
        "disabled": null,
        "cancel": null
    }

    return (
        <div className="file-upload-section">
            <GlassyContainer>
                <div className="file-uploader-container">
                    <div className="file-input">
                        <label htmlFor="fileInput">Upload Text File</label>
                        <input
                            id="fileInput"
                            type="file"
                            accept=".txt"
                            onChange={handleFileChange}
                            disabled={status !== "init"} 
                        />
                        {fileName && (
                            <div className="file-name">
                                {fileName} ({totalLines} Candidates)
                            </div>
                        )}
                    </div>
                    <div className="file-preview">
                        {lines.map((line, index) => (
                            <p key={index}>{index + 1}. {line}</p>
                        ))}
                        {lines.length > 0 ? <p>...</p> : <p>No content found !</p>}
                    </div>
                    {status === "init" ?
                        <InterButtons
                            labels={labels}
                            actors={actors}
                            status={fileObj ? 'enabled' : 'disabled'}
                        />
                        : ""}

                </div>
            </GlassyContainer>
        </div>
    );
}