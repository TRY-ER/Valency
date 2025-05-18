import React, { useState } from "react";
import "./MolStarViewer.css";
import { motion } from "framer-motion";
import { fadeInUpVariantStatic } from "../../../components/animations/framerAnim";
import GlassyContainer from "../../../components/glassy_container/gc";

const MolStarViewer = () => {
    const [activeMol, setActiveMol] = useState("7bv2"); // Example PDB ID
    const [isValidMol, setIsValidMol] = useState(false);

    return <>
        {/* <motion.div 
            initial="hidden"
            animate="visible"
            variants={fadeInUpVariantStatic} 
        className="explore-container">
            <div className="explorer-row-1">
                <Molstar />
            </div>
        </motion.div> */}
        <motion.div
            initial="hidden"
            animate="visible"
            variants={fadeInUpVariantStatic}
            className="viewer-container"
        >
            <GlassyContainer className="viewer-glassy-container">
                <div className="viewer-section">
                    <iframe src="https://molstar.org/me/viewer/?power-preference=high-performance" height="100%" width="100%" title="Mesoscale Explorer"></iframe>
                </div>
            </GlassyContainer>
        </motion.div>
    </>
}

export default MolStarViewer;