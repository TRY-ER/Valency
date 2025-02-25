import React, {  useRef, useState } from "react";
import FileDownloader from "../../../components/file_handler/FileDownloader";
import { call_endpoint_async, call_eventsource } from "../../../endpoints/caller";
import { generate_endpoints } from "../../../endpoints/endpoints";
import ProgressLoader from "../../../components/progress_loader/ProgressLoader";
import GenNumSelector from "../../../components/UI/genNumSelector/genNumSelector";

const LSTMComponent = ({ inputType = "psmiles" }) => {
    const [generatedData, setGeneratedData] = useState([]);
    const [status, setStatus] = useState("init");
    const [consoleData, setConsoleData] = useState([]);
    const [currentNum, setCurrentNum] = useState(0);
    const [maxNum, setMaxNum] = useState(0);
    const currentNumRef = useRef(currentNum);


    const incrementCurrentNum = () => {
        setCurrentNum((prev) => {
          const newNum = prev + 1;
          currentNumRef.current = newNum;
          return newNum;
        });
      };
      

    const handleStartProcess = async () => {
        try {
            const generate_endpoint = generate_endpoints.lstm.set
            const payload = {
                input_type: inputType,
                num_gen: maxNum
            }
            // Upload file first and get ID
            const setResponse = await call_endpoint_async(
                generate_endpoint,
                payload
            );
            console.log("setResponse >>", setResponse)
            setStatus("running");

            if (setResponse.data.id) {
                setConsoleData(prevData => [...prevData, <><p className="console-normal">
                    Generation initiated: Job ID {setResponse.data.id}
                </p><br /></>]);
                const gen_stream_endpoint = generate_endpoints.lstm.stream

                // Start EventSource connection
                const eventSource = call_eventsource(
                    gen_stream_endpoint,
                    setResponse.data.id
                );

                eventSource.onmessage = (event) => {
                    const cleanJsonString = event.data.replace(/'/g, '"');
                    console.log("cleanJsonString >>", cleanJsonString)
                    const data = JSON.parse(cleanJsonString);
                    console.log("data >>", data)
                    if (data.type === "generated") {
                        if (currentNumRef.current < 6) {
                            setConsoleData(prevData => [...prevData, <p className="console-normal">
                                {currentNumRef.current + 1}:  {data.data}
                            </p>]);
                        }
                        else if (currentNumRef.current === 6){
                            setConsoleData(prevData => [...prevData, <p className="console-normal">
                                ...
                            </p>]);
                        }
                        incrementCurrentNum();
                    }
                    else if (data.type === "completed"){
                        setGeneratedData(data.data);
                        eventSource.close();
                        setStatus("completed");
                    }
                    else if (data.type === "error"){
                        setConsoleData(prevData => [...prevData, <p className="console-error">
                            {data.data}
                        </p>]);
                        eventSource.close();
                        setStatus("failed");
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
        setStatus("init");
        setCurrentNum(0);
        setMaxNum(0);
        setGeneratedData([]);
        setConsoleData([]);
    };

    return (
        <>
            <div className="generator-container">
                <div className="generator-row-1">
                    <ProgressLoader progressNum={currentNum} maxNum={maxNum} status={status} />
                </div>
                <div className="generator-row-2">
                    <GenNumSelector
                        status={status}
                        onStartProcess={handleStartProcess}
                        header={"Number of Generations"}
                        setGenNum={setMaxNum}
                        genNum={maxNum}
                        placeholder={"Type Number of Generations"}
                    />
                    <FileDownloader
                        status={status}
                        generatedData={generatedData}
                        consoleContent={consoleData}
                        handleReset={handleReset}
                        showGenerated={false}
                        expandable={true}
                    />
                </div>
            </div>
        </>
    );
};

export default LSTMComponent;