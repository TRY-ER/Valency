// setup react component

import React from 'react';
import './Explorer.css';
import { Outlet, NavLink } from 'react-router-dom';
import Divider from '../../components/divider';
import { motion } from 'framer-motion';
import { fadeInDownVariants } from '../../components/animations/framerAnim'; 

export default function Explorer({
    tabContent,basePath=""
}) {
    // const [activeTabId, setActiveTabId] = useState(1);

    return (
        <div className="base-page-container">
            <div className="tab-container">
                {
                    tabContent.map((tab, index) => {
                        return (
                            <NavLink
                                key={tab.id}
                                to={tab.link}
                                end={tab.link === ''} // Ensures exact match for the index tab
                                className={({ isActive }) => `tab-link ${isActive ? "active" : ""}`}
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
            {/* <Routes>
                <Route index element={
                    <FunctionalSection
                        docElem={ExploreTabContent[0].docs}
                        funcElem={ExploreTabContent[0].component}
                    />
                } />
                {
                    ExploreTabContent.map((tab, index) => {
                        return tab.link && (
                            <Route path={tab.link} element={<FunctionalSection
                                docElem={tab.docs}
                                funcElem={tab.component}
                            />} />
                        )
                    })
                }
            </Routes> */}
            {/* <FunctionalSection
                docElem={ExploreTabContent[activeTabId - 1].docs}
                funcElem={ExploreTabContent[activeTabId - 1].component}
            /> */}
        </div>
    )
}