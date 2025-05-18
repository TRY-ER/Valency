import React, {useState} from "react";
// import "./ProtExplorer.css"; // REMOVED
import "./UniProtViewer.css"; // ADDED - to get .explorer-row-2 and other relevant styles
import MolInputBox from "../../../components/UI/InputBox/InputBox";
import InfoBox from "../../../components/UI/InfoBox/InfoBox";
import ThreeDViewer from "../../../components/UI/ThreeDViewer/ThreeDViewer";
import { motion } from "framer-motion";
import { fadeInUpVariantStatic } from "../../../components/animations/framerAnim";

const BasicProtViewer = () => {
    const [activeMol, setActiveMol] = useState("");
    const [isValidMol, setIsValidMol] = useState(false);

    return <>
        <motion.div 
            initial="hidden"
            animate="visible"
            variants={fadeInUpVariantStatic} 
            className="explore-container">
            <div className="explorer-row-1">
                {/* MolInputBox is used directly as per original structure */}
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
                {/* InfoBox is wrapped in a div for flex layout; it contains its own GlassyContainer. */}
                <div style={{flex: 1, marginRight: '10px'}}>
                    <InfoBox activeMol={activeMol} isValidMol={isValidMol} infoType={"PROT"}/>
                </div>
                {/* ThreeDViewer is wrapped in a div for flex layout */}
                <div style={{flex: 2}}> 
                    <ThreeDViewer activeMol={activeMol} isValidMol={isValidMol} setIsValidMol={setIsValidMol} visType={"PROT"}/>
                </div>
            </div>
        </motion.div>
    </>
}

export default BasicProtViewer;

