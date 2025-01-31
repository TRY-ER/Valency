import React, {useState, useEffect, isValidElement} from "react";
import "./TwoDViewer.css";

import GlassyContainer from "../../glassy_container/gc";
import { call_endpoint_async } from "../../../endpoints/caller";
import { endpoints } from "../../../endpoints/endpoints";

const TwoDViewer = ({
    activeMol,
    isValidMol,
    visType
}) => {
    const [imageData, setImageData] = useState(null);
    useEffect(() =>{
        if(activeMol.length > 0 && isValidMol){
            const payload = {
                type: visType,
                value: activeMol
            }

            console.log("payload >>", payload)

            call_endpoint_async(endpoints.visualize, payload).then((response) => {
                if(response.data.status === "success"){
                    setImageData(response.data.image_data);
                }
            }).catch((error) => {
                console.log('error in visualizing !>>', error);
            })
        } 
        else{
            setImageData(null);
        }
    }, [activeMol, isValidMol])  

    return <>
        <div className="twod-container">
            <GlassyContainer expandable={true} className="twod-viewer" mode={activeMol.length > 0 && !isValidMol ? "loading" : "enabled"}>
                {
                    imageData === null || activeMol.length == 0 ? 
                    <p></p> :
                    <img src={`data:image/png;base64,${imageData}`} alt="2D Molecule Image" />
                }
            </GlassyContainer>
        </div>
    </>
}

export default TwoDViewer; 