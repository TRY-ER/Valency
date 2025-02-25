import { React } from 'react';
import "./TabSection.css";
import { NavLink} from 'react-router-dom';

const TabContainer = ({ tabDetails, basePath = "" }) => {
    
    return (
        <>
            <div className="tab-container">
                {
                    tabDetails.map((tab) => {
                        return (
                            <NavLink 
                             key={tab.id} 
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
        </>
    )
}

export default TabContainer