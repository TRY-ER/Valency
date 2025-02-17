import React, { useEffect, useState } from 'react';
import TabContainer from '../../../components/tab_section/TabSection';
import FunctionalSection from '../../../components/functional_section/Functional';
import Divider from '../../../components/divider';
import { BRICSGeneratorTabContents } from '../../../contents/tag_content/NestGeneratorTags';
import { useLocation, useNavigate, NavLink, Outlet } from 'react-router-dom';
import { motion } from 'framer-motion';
import { fadeInDownVariants } from '../../../components/animations/framerAnim';
import "./Generator.css";

export default function BRICSGenerator({
    tabContent, basePath = "generators"
}) {
    const location = useLocation();
    // const location = useLocation();
    // const navigate = useNavigate();

    // useEffect(() => {
    //     const queryParams = new URLSearchParams(location.search);
    //     const tabParam = parseInt(queryParams.get('sub-tab'), 10);

    //     if (tabParam >= 1 && tabParam <= BRICSGeneratorTabContents.length) {
    //         setActiveTabId(tabParam);
    //     }
    // }, [location.search]);

    // useEffect(() =>{
    //     const queryParams = new URLSearchParams(location.search);
    //     queryParams.set('sub-tab', activeTabId);
    //     navigate({
    //         pathname: location.pathname,
    //         search: queryParams.toString(),
    //     });
    // }, [activeTabId])

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
            {/* <FunctionalSection
                docElem={BRICSGeneratorTabContents[activeTabId - 1].docs}
                funcElem={BRICSGeneratorTabContents[activeTabId - 1].component}
                customClassName={"functional-container-nested"}
            /> */}
            {/* {functionalContent}  */}
        </>
    )
}