import React, {useState, useEffect} from "react";
import "./InfoBox.css";

import GlassyContainer from "../../glassy_container/gc";
import { call_endpoint_async } from "../../../endpoints/caller";
import { endpoints } from "../../../endpoints/endpoints";

const InfoBox = ({
    activeMol,
    isValidMol,
    infoType
}) => {
    const [infoVal, setInfoVal] = useState(null);
    useEffect(() =>{
        if(activeMol.length > 0 && isValidMol){
            const payload = {
                type: infoType, 
                value: activeMol
            }

            call_endpoint_async(endpoints.inform, payload).then((response) => {
                if(response.data.status === "success"){
                    setInfoVal(response.data.info);
                }
            }).catch((error) => {
                console.log('error in getting info val!>>', error);
            })
        }
        else{
            setInfoVal(null);
        }
    }, [activeMol, isValidMol])

    useEffect(() =>{
        console.log("info val >>", infoVal);
    }, [infoVal])

    const renderInfoContent = () => {
        if (!infoVal) return <p></p>;
        
        return Object.entries(infoVal).map(([key, value]) => (
            <div key={key} className="info-row">
                <span className="info-label">{key.replace(/_/g, ' ')}</span>
                <span className="info-value">{value}</span>
            </div>
        ));
    };

    return <>
        <div className="infobox-container">
            <GlassyContainer mode={activeMol.length > 0 && !isValidMol ? "loading" : "enabled"}>
                {renderInfoContent()}
            </GlassyContainer>
        </div>
    </>
}

export default InfoBox;