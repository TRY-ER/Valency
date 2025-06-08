import React, { useState } from "react";
import { motion } from "framer-motion";
import { fadeInUpVariantStatic, fadeInDownVariants } from "../../../components/animations/framerAnim";
import GlassyContainer from "../../../components/glassy_container/gc";
import Divider from "../../../components/divider";
import "./PubChemUtilities.css";

// Import PubChem utility subcomponents
import CompoundSearch from "./UtilityComponents/CompoundSearch";
import SubstructureSearch from "./UtilityComponents/SubstructureSearch";
import MassSearch from "./UtilityComponents/MassSearch";
import XrefSearch from "./UtilityComponents/XrefSearch";
import CompoundSynonyms from "./UtilityComponents/CompoundSynonyms";
import CompoundProperties from "./UtilityComponents/CompoundProperties";

// Utility tabs configuration
const utilityTabs = [
    {
        id: 1,
        title: 'Substructure Search',
        key: 'substructure-search',
        component: <SubstructureSearch />,
        description: 'Find compounds containing specific substructures using SMILES patterns.'
    },
    {
        id: 2,
        title: 'Compound Synonyms',
        key: 'compound-synonyms',
        component: <CompoundSynonyms />,
        description: 'Retrieve synonyms and alternative names for PubChem compounds.'
    },
    {
        id: 3,
        title: 'Compound Properties',
        key: 'compound-properties',
        component: <CompoundProperties />,
        description: 'Get computed properties and descriptors for compounds.'
    }
];

const PubChemUtilities = () => {
    const [activeTabId, setActiveTabId] = useState(1);
    
    const activeTab = utilityTabs.find(tab => tab.id === activeTabId);

    const handleTabClick = (tabId) => {
        setActiveTabId(tabId);
    };

    return (
        <>
            <motion.div
                className="pubchem-utilities-container"
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

export default PubChemUtilities;
