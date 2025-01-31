import React, {useState, useEffect} from "react";
import "./ProtExplorer.css";

import MolInputBox from "../../../components/UI/InputBox/InputBox";
import InfoBox from "../../../components/UI/InfoBox/InfoBox";
import TwoDViewer from "../../../components/UI/TwoDViewer/TwoDViewer";
import ThreeDViewer from "../../../components/UI/ThreeDViewer/ThreeDViewer";

const ProtEComponent = () => {
    const [activeMol, setActiveMol] = useState("");
    const [isValidMol, setIsValidMol] = useState(false);

    return <>
        <div className="explore-container">
            <div className="explorer-row-1">
                <MolInputBox 
                    activeMol={activeMol}
                    setActiveMol={setActiveMol} 
                    isValidMol={isValidMol}
                    setIsValidMol={setIsValidMol}
                    inputType={"PROT"}
                    placeholder="Try PDB ID: 1MO8"
                    header="Enter PDB Id Here"
                />
            </div>
            <div className="explorer-row-2">
                <InfoBox activeMol={activeMol} isValidMol={isValidMol} infoType={"PROT"}/>
                <ThreeDViewer activeMol={activeMol} isValidMol={isValidMol} setIsValidMol={setIsValidMol} visType={"PROT"}/>
            </div>
        </div>
    </>
}

export default ProtEComponent;