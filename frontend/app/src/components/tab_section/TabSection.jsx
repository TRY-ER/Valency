import { React } from 'react';
import "./TabSection.css";

const TabContainer = ({ tabDetails, activeTabId, setActiveTabId }) => {
    return (
        <>
            <div className="tab-container">
                {
                    tabDetails.map((tab) => {
                        return (
                            <div className={`tab-tag glassy-feel ${tab.id === activeTabId ? "active" : ""}`}
                                onClick={() => setActiveTabId(tab.id)}>
                                <p className="tab-tag-text">{tab.title}</p>
                            </div>
                        )
                    })
                }
            </div>
        </>
    )
}

export default TabContainer