import React, { useEffect, useRef, useState } from "react";
import ProgressSection from "../../../components/progress_section/ProgressSection";
import FileUploader from "../../../components/file_handler/FileUploader";
import FileDownloader from "../../../components/file_handler/FileDownloader";
import { call_endpoint_async, call_eventsource } from "../../../endpoints/caller";
import { generate_endpoints } from "../../../endpoints/endpoints";
import ProgressLoader from "../../../components/progress_loader/ProgressLoader";
import GenNumSelector from "../../../components/UI/genNumSelector/genNumSelector";
import NumValidSelector from "../../../components/UI/numValidSelector/numValidSelector";
import SimilarityCardContainer from "../../../components/UI/similarityCardContainer/SimilarityCardContainer";
import { discriminate_endpoints } from "../../../endpoints/endpoints";

const SSContent = ({ inputType = "MOL" }) => {
    const [generatedData, setGeneratedData] = useState([]);
    const [status, setStatus] = useState("init");
    const [maxNum, setMaxNum] = useState(12);
    const [activeMol, setActiveMol] = useState("");
    const [errorMsg, setErrorMsg] = useState("");

    const handleStartProcess = async () => {
        try{
            setStatus("loading");
            setGeneratedData([]);
            const payload = {
                input_type: inputType,
                data: activeMol,
                k: maxNum 
            };
            const response = await call_endpoint_async(discriminate_endpoints.ssearch, payload);
            if(response.data.status === "success"){
                console.log(response.data.results);
                setGeneratedData(response.data.results);
                setStatus("completed");
            }
            else if(response.data.status === "failed"){
                setStatus("failed");
                setErrorMsg(response.data.message);
            }
        } 
        catch(error){
            setStatus("failed");
        }
    };

    const handleReset = () => {
        setStatus("init");
        setGeneratedData([]);
        setActiveMol("");
        setMaxNum(12);
    };

    const headerMapper = {
        "MOL": "Enter Molecular SMILES Here",
        "POLY": "Enter Polymer SMILES Here",
        "PROT": "Enter Protein PDB Id Here",
    }

    const placeholderMapper  = {
        "MOL": "Try C for demo",
        "POLY": "Try [*]CC[*] for demo",
        "PROT": "Try 1MO8 for demo",
    }

    useEffect(() =>{
        console.log("generated data >>", generatedData);
    }, [generatedData])

    
    return (
        <>
            <div className="generator-container">
                <div className="generator-row-1">
                    <NumValidSelector
                        status={status}
                        onStartProcess={handleStartProcess}
                        numheader={"Number of Generations (1 - 30)"}
                        setGenNum={setMaxNum}
                        genNum={maxNum}
                        numplaceholder={"Type Number of Generations"}
                        molheader={headerMapper[inputType]}
                        molplaceholder={placeholderMapper[inputType]}
                        inputType={inputType}
                        activeMol={activeMol}
                        setActiveMol={setActiveMol}
                        handleReset={handleReset}
                        generatedData={generatedData}
                    />
                </div>
                <div className="generator-row-2">
                    <SimilarityCardContainer data={generatedData} status={status} />
                </div>
            </div>
        </>
    );
};

export default SSContent;