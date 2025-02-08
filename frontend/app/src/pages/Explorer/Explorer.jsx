// setup react component

import React, { useEffect, useState } from 'react';
import './Explorer.css';
import { Routes, useLocation, useNavigate, Route, Outlet, NavLink } from 'react-router-dom';
import GlassyContainer from '../../components/glassy_container/gc';
import TabContainer from '../../components/tab_section/TabSection';
import ExploreTabContent from '../../contents/tag_content/ExploreTags';
import FunctionalSection from '../../components/functional_section/Functional';
import Divider from '../../components/divider';
import MoEDocs from '../../contents/doc_content/explorer_content/MoEDocs';

export default function Explorer({
    tabContent,basePath=""
}) {
    // const [activeTabId, setActiveTabId] = useState(1);

    return (
        <div className="base-page-container">
            <div className="tab-container">
                {
                    tabContent.map((tab) => {
                        return (
                            <NavLink
                                key={tab.id}
                                to={`${basePath}/${tab.link}`}
                                className={({ isActive }) => `tab-link ${isActive ? "active" : ""}`}
                            >
                                <div className={`tab-tag glassy-feel`}>
                                    <p className="tab-tag-text">{tab.title}</p>
                                </div>
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