import React, { useEffect, useState } from "react";
import ProgressSection from "../../../components/progress_section/ProgressSection";
import FileUploader from "../../../components/file_handler/FileUploader";
import FileDownloader from "../../../components/file_handler/FileDownloader";
import { call_endpoint_async, call_eventsource } from "../../../endpoints/caller";
import { generate_endpoints } from "../../../endpoints/endpoints";


const BRICSComponent = ({ inputType = "smiles" }) => {
    const [generatedData, setGeneratedData] = useState([]);
    const [total, setTotal] = useState(0);
    const [rawFile, setRawFile] = useState(null);
    const [status, setStatus] = useState("init");
    const [consoleData, setConsoleData] = useState([]);
    const [processSteps, setProcessSteps] = useState([
        {
            id: 1,
            title: "Upload",
            status: "init"
        },
        {
            id: 2,
            title: "Decompose",
            status: "init"
        },
        {
            id: 3,
            title: "Compose",
            status: "init"
        },
        {
            id: 4,
            title: "Output",
            status: "init"
        }]);

    const handleStartProcess = async () => {
        if (!rawFile) return;
        
        setStatus("running");
        // Update status
        setProcessSteps((prev) => {
            const newSteps = [...prev];
            newSteps[0].status = "completed"
            newSteps[1].status = "pending";
            newSteps[2].status = "pending";
            newSteps[3].status = "pending";
            return newSteps;
        });

        try {
            // Create form data
            const formData = new FormData();
            formData.append("file", rawFile);
            var generate_endpoint = ""

            if (inputType === "smiles"){
                generate_endpoint = generate_endpoints.brics.smiles
            }
            else if (inputType === "psmiles"){
                generate_endpoint = generate_endpoints.brics.psmiles
            }
            // Upload file first and get ID
            const uploadResponse = await call_endpoint_async(
                generate_endpoint,
                formData
            );

            if (uploadResponse.data.file_id) {
                // Update upload step
                const newSteps = [...processSteps];
                newSteps[0].status = "completed";
                setProcessSteps(newSteps);
                
                var gen_stream_endpoint= ""
                if (inputType === "smiles"){
                    gen_stream_endpoint = generate_endpoints.brics.smiles_stream
                }
                else if (inputType === "psmiles"){
                    gen_stream_endpoint = generate_endpoints.brics.psmiles_stream
                }
                // Start EventSource connection
                const eventSource = call_eventsource(
                    gen_stream_endpoint,
                    uploadResponse.data.file_id
                );

                eventSource.onmessage = (event) => {
                    const cleanJsonString = event.data.replace(/'/g, '"');
                    console.log("cleanJsonString >>", cleanJsonString)
                    const data = JSON.parse(cleanJsonString);

                    // Update steps based on event type
                    if (data.step === "decomposition") {
                        console.log("decompsition occuring")
                        setProcessSteps((prev) => {
                            const newSteps = [...prev];
                            newSteps[1].status = "completed";
                            return newSteps;
                        });
                        setConsoleData(prevData => [...prevData, <p className="console-normal">
                            Decomposition Done, {data.data} substructures were found
                        </p>]);
                    } else if (data.step === "composition") {
                        console.log("compsition occuring")
                        setProcessSteps((prev) => {
                            const newSteps = [...prev];
                            newSteps[2].status = "completed";
                            return newSteps;
                        });
                        setConsoleData(prevData => [...prevData, <p className="console-normal">
                            Composition Done, {data.data} candidates formed 
                        </p>]);
                    } else if (data.total !== undefined) {
                        console.log("got total count")
                        setTotal(data.total);
                        setConsoleData(prevData => [ ...prevData, <p className="console-normal">
                            Filtration Done, {data.total} candidates selected. 
                        </p>]);
                    } else if (data.hasOwnProperty('file')) {
                        console.log("got composed file")
                        // Decode and set generated data
                        const decodedData = atob(data.file);
                        const splittedData = decodedData.split("\n");
                        // make sure that empty strings are removed from the list
                        const filteredData = splittedData.filter(Boolean); 
                        setGeneratedData(filteredData);
                        setProcessSteps((prev) => {
                            const newSteps = [...prev];
                            newSteps[3].status = "completed";
                            return newSteps;
                        });
                        eventSource.close();
                        setStatus("completed");
                        
                    }
                    else if (data.type == "error") {
                        setProcessSteps((prev) => {
                            const newSteps = [...prev];
                            const pendingIndex = newSteps.findIndex(step => step.status === "pending");
                            if (pendingIndex !== -1) {
                                newSteps[pendingIndex].status = "failed";
                            }
                            return newSteps;
                        });
                        eventSource.close();
                        setConsoleData(prevData => [...prevData, <p className="console-alert">
                            Error in the Process: {data.data}
                        </p>]);
                        setStatus("failed");
                        return;
                    }
                };

                eventSource.onerror = () => {
                    eventSource.close();
                };
            }
        } catch (error) {
            console.error("Error in process:", error);
        }
    };

    const handleReset = () => {
        // Reset all states to initial values
        setRawFile(null);
        setStatus("init");
        setConsoleData([]);
        setGeneratedData(null);
        setTotal(0);
        setProcessSteps([
            {
                id: 1,
                title: "Upload",
                status: "init"
            },
            {
                id: 2,
                title: "Decompose",
                status: "init"
            },
            {
                id: 3,
                title: "Compose",
                status: "init"
            },
            {
                id: 4,
                title: "Output",
                status: "init"
            }
        ]);
    };    

    return (
        <>
            <div className="generator-container">
                <div className="generator-row-1">
                    <ProgressSection processSteps={processSteps} />
                </div>
                <div className="generator-row-2">
                    <FileUploader
                        status={status}
                        setRawFile={setRawFile}
                        onStartProcess={handleStartProcess}
                    />
                    <FileDownloader
                        status={status}
                        generatedData={generatedData}
                        consoleContent={consoleData}
                        handleReset={handleReset}
                    />
                </div>
            </div>
        </>
    );
};

export default BRICSComponent;