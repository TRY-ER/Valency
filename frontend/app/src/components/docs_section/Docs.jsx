import React, { useState } from 'react';
import GlassyContainer from '../glassy_container/gc';
import { FaCompressArrowsAlt, FaExpandArrowsAlt } from "react-icons/fa";
import { IoDocumentText } from "react-icons/io5";
import "./Docs.css";

const DocsContainer = ({ docElem, isCollapsed, onToggleCollapse, onExpand }) => {
    const [isFullScreen, setIsFullScreen] = useState(false);

    const handleExpand = () => {
        setIsFullScreen(true);
        if (onExpand) onExpand();
    };

    const handleCloseFullScreen = () => {
        setIsFullScreen(false);
    };

    // If collapsed, show only the documentation button
    if (isCollapsed) {
        return (
            <div className="docs-collapsed-container">
                <button 
                    className="docs-button"
                    onClick={onToggleCollapse}
                    title="Show Documentation"
                >
                    <IoDocumentText />
                    <span>Docs</span>
                </button>
            </div>
        );
    }

    return (
        <>
            <div className="docs-container-wrapper">
                <div className="docs-header">
                    <h3>Documentation</h3>
                    <div className="docs-controls">
                        <button 
                            className="docs-control-btn expand-btn"
                            onClick={handleExpand}
                            title="Expand to full screen"
                        >
                            <FaExpandArrowsAlt />
                        </button>
                        <button 
                            className="docs-control-btn collapse-btn"
                            onClick={onToggleCollapse}
                            title="Collapse documentation panel"
                        >
                            <FaCompressArrowsAlt />
                        </button>
                    </div>
                </div>
                <div className="docs-content-scroll">
                    <GlassyContainer expandable={false}>
                        <div className="doc-container">
                            {docElem}
                        </div>
                    </GlassyContainer>
                </div>
            </div>

            {/* Full screen modal */}
            {isFullScreen && (
                <div className="docs-modal-overlay" onClick={handleCloseFullScreen}>
                    <div className="docs-modal-content" onClick={(e) => e.stopPropagation()}>
                        <div className="docs-modal-header">
                            <h3>Documentation</h3>
                            <button 
                                className="docs-close-btn"
                                onClick={handleCloseFullScreen}
                                title="Close full screen"
                            >
                                ×
                            </button>
                        </div>
                        <div className="docs-modal-body">
                            {docElem}
                        </div>
                    </div>
                </div>
            )}
        </>
    );
};

export default DocsContainer;