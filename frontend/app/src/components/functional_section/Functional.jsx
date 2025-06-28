import React, { useState } from 'react';
import "./Functional.css";
import DocsContainer from '../docs_section/Docs';
import FunctionContainer from '../function_container/FunctionContainer';

const FunctionalSection = ({ docElem, funcElem, customClassName = null }) => {
    const [isDocsCollapsed, setIsDocsCollapsed] = useState(false);

    const handleToggleDocsCollapse = () => {
        setIsDocsCollapsed(!isDocsCollapsed);
    };

    const handleDocsExpand = () => {
        // When docs expand to full screen, we might want to handle it differently
        // For now, we'll keep it simple
    };

    return (
        <div className={`${customClassName ? customClassName : "functional-container"} ${isDocsCollapsed ? "docs-collapsed" : ""}`}>
            <div className="functional-item">
                <FunctionContainer functionalComponents={funcElem} />
            </div>
            {
                docElem && <div className={`functional-doc ${isDocsCollapsed ? "collapsed" : ""}`}>
                    <DocsContainer 
                        docElem={docElem} 
                        isCollapsed={isDocsCollapsed}
                        onToggleCollapse={handleToggleDocsCollapse}
                        onExpand={handleDocsExpand}
                    />
                </div>
            }
        </div>
    )
}

export default FunctionalSection;