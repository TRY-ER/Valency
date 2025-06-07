import React, { useState } from "react";
import { motion } from "framer-motion";
import { fadeInUpVariantStatic, fadeInDownVariants } from "../../../components/animations/framerAnim";
import GlassyContainer from "../../../components/glassy_container/gc";
import Divider from "../../../components/divider";
import "./ChemBLUtilities.css";

// Import dummy subcomponents
import TargetByGeneName from "./UtilityComponents/TargetByGeneName";
import SmilesToCTAB from "./UtilityComponents/SmilesToCTAB";
import ComputeMolecularDescriptor from "./UtilityComponents/ComputeMolecularDescriptor";
import ComputeStructuralAlerts from "./UtilityComponents/ComputeStructuralAlerts";
import StandardizeMolecules from "./UtilityComponents/StandardizeMolecules";
import ParentMolecule from "./UtilityComponents/ParentMolecule";

// Utility tabs configuration
const utilityTabs = [
    {
        id: 1,
        title: 'Target by Gene Name',
        key: 'target-gene',
        component: <TargetByGeneName />,
        description: 'Find biological targets using gene names from ChEMBL database.'
    },
    {
        id: 2,
        title: 'SMILES to CTAB',
        key: 'smiles-ctab',
        component: <SmilesToCTAB />,
        description: 'Convert SMILES notation to CTAB format for molecular structures.'
    },
    {
        id: 3,
        title: 'Compute Molecular Descriptor',
        key: 'molecular-descriptor',
        component: <ComputeMolecularDescriptor />,
        description: 'Calculate molecular descriptors and properties from chemical structures.'
    },
    {
        id: 4,
        title: 'Compute Structural Alerts',
        key: 'structural-alerts',
        component: <ComputeStructuralAlerts />,
        description: 'Identify structural alerts and potential toxicity patterns in molecules.'
    },
    {
        id: 5,
        title: 'Standardize Molecules from SMILES',
        key: 'standardize-molecules',
        component: <StandardizeMolecules />,
        description: 'Standardize and normalize molecular structures from SMILES strings.'
    },
    {
        id: 6,
        title: 'Parent Molecule from SMILES',
        key: 'parent-molecule',
        component: <ParentMolecule />,
        description: 'Extract parent molecule structures from SMILES representations.'
    }
];

const ChemBLUtilities = () => {
    const [activeTabId, setActiveTabId] = useState(1);
    
    const activeTab = utilityTabs.find(tab => tab.id === activeTabId);

    const handleTabClick = (tabId) => {
        setActiveTabId(tabId);
    };

    return (
        <>
            <motion.div
                className="chembl-utilities-container"
                variants={fadeInUpVariantStatic}
                initial="initial"
                animate="animate"
            >
                {/* Tab Navigation */}
                <div className="utilities-tab-container">
                    {utilityTabs.map((tab, index) => (
                        <motion.div
                            key={tab.id}
                            initial="hidden"
                            animate="visible"
                            variants={fadeInDownVariants}
                            custom={index}
                            className={`utilities-tab ${activeTabId === tab.id ? 'active' : ''}`}
                            onClick={() => handleTabClick(tab.id)}
                        >
                            <div className="utilities-tab-tag glassy-feel">
                                <p className="utilities-tab-text">{tab.title}</p>
                            </div>
                        </motion.div>
                    ))}
                </div>

                <Divider />

                {/* Active Tab Description */}
                {/* {activeTab && (
                    <div className="utilities-tab-description">
                        <GlassyContainer>
                            <h3 className="utilities-tab-title">{activeTab.title}</h3>
                            <p className="utilities-tab-desc">{activeTab.description}</p>
                        </GlassyContainer>
                    </div>
                )} */}

                {/* Active Component */}
                <div className="utilities-content">
                    {activeTab && activeTab.component}
                </div>
            </motion.div>
        </>
    );
};

export default ChemBLUtilities;
