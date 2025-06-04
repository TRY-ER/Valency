import React, { useState } from "react"; // Removed useContext
import "./MolExplorer.css";

import MolInputBox from "../../../components/UI/InputBox/InputBox";
import InfoBox from "../../../components/UI/InfoBox/InfoBox";
import TwoDViewer from "../../../components/UI/TwoDViewer/TwoDViewer";
import { motion } from "framer-motion";
import { fadeInUpVariantStatic } from "../../../components/animations/framerAnim";
// Removed ThemeContext import

const MolEComponent = () => {
    const [activeMol, setActiveMol] = useState("");
    const [isValidMol, setIsValidMol] = useState(false);
    // Removed theme consumption

    return <>
        <motion.div
            initial="hidden"
            animate="visible"
            variants={fadeInUpVariantStatic}
            className="explore-container">
            <div className="explorer-row-1">
                <MolInputBox
                    activeMol={activeMol}
                    setActiveMol={setActiveMol}
                    isValidMol={isValidMol}
                    setIsValidMol={setIsValidMol}
                    inputType={"MOL"}
                />
            </div>
            <div className="explorer-row-2">
                <InfoBox activeMol={activeMol} isValidMol={isValidMol} infoType={"MOL"} />
                <TwoDViewer activeMol={activeMol} isValidMol={isValidMol} visType={"MOL"} />
            </div>
        </motion.div>
    </>
}

export default MolEComponent;