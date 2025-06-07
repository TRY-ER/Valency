import { React } from 'react';
import { motion } from "framer-motion";
import "./TabSection.css";
import { NavLink } from 'react-router-dom';
import { fadeInDownVariants } from "../animations/framerAnim";

const TabContainer = ({ tabDetails, basePath = "" }) => {
    
    return (
        <>
            <div className="tab-container">
                {
                    tabDetails.map((tab, index) => {
                        return (
                            <motion.div
                                key={tab.id}
                                initial="hidden"
                                animate="visible"
                                variants={fadeInDownVariants}
                                custom={index}
                            >
                                <NavLink 
                                    className={({ isActive }) => `tab-link ${isActive ? "active" : ""}`}
                                    to={basePath ? `${basePath}/${tab.path || tab.key || ''}` : (tab.path || tab.key || '')}
                                > 
                                    <div className={`tab-tag glassy-feel`}>
                                        <p className="tab-tag-text">{tab.title}</p>
                                    </div>
                                </NavLink>
                            </motion.div>
                        )
                    })
                }
            </div>
        </>
    )
}

export default TabContainer