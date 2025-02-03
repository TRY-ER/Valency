import React, {useState, useEffect, isValidElement} from "react";
import "./ThreeDViewer.css";

import GlassyContainer from "../../glassy_container/gc";
import { call_endpoint_async } from "../../../endpoints/caller";
import { endpoints } from "../../../endpoints/endpoints";
import ProtViewer from "../ProtViewer/ProtViewer";

const ThreeDViewer = ({
    activeMol,
    isValidMol,
    setIsValidMol,
    visType
}) => {
    const [imageData, setImageData] = useState(null);
    
    return <>
        <div className="twod-container">
            <GlassyContainer expandable={true} className="twod-viewer" > 
                {
                    activeMol.length > 0 &&
                    <ProtViewer activeMol={activeMol} setIsValidMol={setIsValidMol} isValidMol={isValidMol} />
                }
            </GlassyContainer>
        </div>
    </>
}

export default ThreeDViewer; 