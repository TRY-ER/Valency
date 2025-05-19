import React, {useState} from "react"; // Removed useContext
import "./PolyExplorer.css";

import MolInputBox from "../../../components/UI/InputBox/InputBox";
import InfoBox from "../../../components/UI/InfoBox/InfoBox";
import TwoDViewer from "../../../components/UI/TwoDViewer/TwoDViewer";
import { motion } from "framer-motion";
import { fadeInUpVariantStatic } from "../../../components/animations/framerAnim";
// Removed ThemeContext import

const PolyEComponent = () => {
    const [activeMol, setActiveMol] = useState("");
    const [isValidMol, setIsValidMol] = useState(false);
    // Removed theme consumption

    return <>
        <motion.div 
            initial="hidden"
            animate="visible"
            variants={fadeInUpVariantStatic}
        className="explore-container"> {/* Reverted className */}
            <div className="explorer-row-1">
                <MolInputBox 
                    activeMol={activeMol}
                    setActiveMol={setActiveMol} 
                    isValidMol={isValidMol}
                    setIsValidMol={setIsValidMol}
                    inputType={"POLY"}
                    header = "Enter Polymer SMILES Here"
                    placeholder = "Try [*]CC[*] for demo"
                />
            </div>
            <div className="explorer-row-2">
                <InfoBox activeMol={activeMol} isValidMol={isValidMol} infoType={"POLY"}/>
                <TwoDViewer activeMol={activeMol} isValidMol={isValidMol} visType={"POLY"}/>
            </div>
        </motion.div>
    </> 
}

export default PolyEComponent;