import React, { useEffect, useState } from "react";
import { call_endpoint_async } from "../../../endpoints/caller";
import NumValidSelector from "../../../components/UI/numValidSelector/numValidSelector";
import SimilarityCardContainer from "../../../components/UI/similarityCardContainer/SimilarityCardContainer";
import { discriminate_endpoints } from "../../../endpoints/endpoints";
import { motion } from "framer-motion";
import { fadeInUpVariantStatic } from "../../../components/animations/framerAnim";

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
            <motion.div 
                initial="hidden"
                animate="visible"
                variants={fadeInUpVariantStatic} 
            className="generator-container">
                <div className="generator-row-1">
                    <NumValidSelector
                        status={status}
                        onStartProcess={handleStartProcess}
                        numheader={"Number of Candidates Required (1 - 30)"}
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
            </motion.div>
        </>
    );
};

export default SSContent;