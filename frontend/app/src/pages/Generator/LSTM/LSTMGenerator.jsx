import React, { useEffect, useState } from 'react';
import TabContainer from '../../../components/tab_section/TabSection';
import FunctionalSection from '../../../components/functional_section/Functional';
import Divider from '../../../components/divider';
import { LSTMGeneratorTabContents } from '../../../contents/tag_content/NestGeneratorTags';
import "./LSTMGenerator.css";
import { useLocation, useNavigate, Outlet, NavLink } from 'react-router-dom';

export default function LSTMGenerator({
    tabContent, basePath = "generators/lstm"
}) {
    const location = useLocation();

    return (
        <>
            <div className="tab-container">
                {
                    tabContent.map((tab) => {
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
       </>
    )
}