import React, {useState, useEffect} from "react";
import TabContainer from "../../../components/tab_section/TabSection";
import Divider from "../../../components/divider";
import FunctionalSection from "../../../components/functional_section/Functional";
import { SimilaritySearchTabContents } from "../../../contents/tag_content/NestDiscriminatorTags";
import { Routes, Route, useLocation, useNavigate, Outlet, NavLink } from 'react-router-dom';
import { fadeInDownVariants } from "../../../components/animations/framerAnim";
import { motion } from "framer-motion";

const SSComponent = ({
    tabContent, basePath = "discriminators"
}) => {
    const location = useLocation();

    return (
        <>
            <div className="tab-container">
                {
                    tabContent.map((tab, index) => {
                        return (
                           <NavLink
                                key={tab.id}
                                to={ tab.link == "" ? `/${basePath}` : `/${basePath}/${tab.link}`}
                                // className={({ isActive }) => `tab-link ${isActive ? "active" : ""}`}
                                className={() =>{
                                    if (location.pathname === `/${basePath}` && tab.link === ""){
                                        return `tab-link active`
                                    }
                                    else if(location.pathname === `/${basePath}/${tab.link}`){
                                        return `tab-link active`
                                    }
                                    else{
                                        return `tab-link`
                                    }
                                }} 
                            >
                                <motion.div 
                                    initial="hidden"
                                    animate="visible"
                                    variants={fadeInDownVariants} 
                                    custom={index}
                                className={`tab-tag glassy-feel`}>
                                    <p className="tab-tag-text">{tab.title}</p>
                                </motion.div>
                            </NavLink>
                        )
                    })
                }
            </div>
            <Divider />
            <Outlet /> 
        </>
    ) 
}

export default SSComponent;