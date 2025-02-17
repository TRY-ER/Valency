import React, { useEffect, useState } from 'react';
import './Discriminator.css';
import { Outlet, useLocation, useNavigate } from 'react-router-dom';
import GlassyContainer from '../../components/glassy_container/gc';
import TabContainer from '../../components/tab_section/TabSection';
import DiscriminatorTabContents from '../../contents/tag_content/DiscriminatorTags';
import FunctionalSection from '../../components/functional_section/Functional';
import Divider from '../../components/divider';
import MoEDocs from '../../contents/doc_content/explorer_content/MoEDocs';
import ExploreTabContents from '../../contents/tag_content/ExploreTags';
import { Routes, Route, NavLink } from 'react-router-dom';
import { motion } from 'framer-motion';
import { fadeInDownVariants } from '../../components/animations/framerAnim';

export default function Discriminator({
    tabContent,
    basePath = "discriminators"
}) {
    const location = useLocation();

    const constructNestPath = (tab) => {
        let allPaths = [];
        if (tab.link === "") {
            allPaths.push(`/${basePath}`);
            if (tab.subElements) {
                for (let subTab of tab.subElements) {
                    if (subTab.link === "") {
                        allPaths.push(`/${basePath}`);
                    }
                    else {
                        allPaths.push(`/${basePath}/${subTab.link}`);
                    }
                }
            }
        } else {
            allPaths.push(`/${basePath}/${tab.link}`);
            if (tab.subElements) {
                for (let subTab of tab.subElements) {
                    if (subTab.link === "") {
                        allPaths.push(`/${basePath}/${tab.link}`);
                    }
                    else {
                        allPaths.push(`/${basePath}/${tab.link}/${subTab.link}`);
                    }
                }
            }
        }
        // console.log("disc nest paths >>", allPaths);
        return allPaths;
    }


    return (
        <div className="base-page-container">
            <div className="tab-container">
                {
                    tabContent.map((tab, index) => {
                        const allNestPaths = constructNestPath(tab);
                        return (
                            <NavLink
                                key={tab.id}
                                to={tab.link === "" ? `/${basePath}` : `/${basePath}/${tab.link}`}
                                className={() =>
                                    allNestPaths.includes(location.pathname)
                                        ? "tab-link active"
                                        : "tab-link"
                                }
                            >
                                <motion.div 
                                variants={fadeInDownVariants} 
                                initial="hidden"
                                animate="visible"
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
        </div>
    )
}