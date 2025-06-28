import React, { useState, useEffect, useRef } from 'react';
import { motion, AnimatePresence } from 'framer-motion';
import { FiX, FiMinimize2, FiMaximize2, FiBook } from 'react-icons/fi';
import { IoDocumentText } from "react-icons/io5";
import DocRenderer from '../../../contents/doc_content/DocRenderer';
import './FloatingDocPanel.css';

const FloatingDocPanel = ({ toolName, isOpen, onClose }) => {
    const [isCollapsed, setIsCollapsed] = useState(false);
    const [isFullScreen, setIsFullScreen] = useState(false);
    const [docPath, setDocPath] = useState(null);
    const [dragConstraints, setDragConstraints] = useState({});

    // Calculate drag constraints based on window size and panel dimensions
    useEffect(() => {
        const updateConstraints = () => {
            const padding = 20;
            const panelWidth = isCollapsed ? 280 : 500;
            const panelHeight = isCollapsed ? 200 : window.innerHeight - 40;
            const windowWidth = window.innerWidth;
            const windowHeight = window.innerHeight;
            
            // Calculate how far the panel can move from its starting position
            let startTop, startRight;
            
            if (isCollapsed) {
                // Collapsed panel starts at 50% from top (centered)
                startTop = windowHeight / 2 - panelHeight / 2; // Actual top position when centered
                startRight = padding; // Distance from right edge
            } else {
                // Expanded panel starts at top: 20px, right: 20px
                startTop = padding; // Distance from top edge
                startRight = padding; // Distance from right edge
            }
            
            setDragConstraints({
                left: -(windowWidth - panelWidth - padding - startRight), // Can move left until padding from left edge
                right: 0, // Cannot move further right than starting position
                top: -startTop + padding, // Can move up until padding from top
                bottom: windowHeight - panelHeight - padding - startTop // Can move down until padding from bottom
            });
        };

        if (isOpen && !isFullScreen) {
            updateConstraints();
            window.addEventListener('resize', updateConstraints);
            return () => window.removeEventListener('resize', updateConstraints);
        }
    }, [isOpen, isCollapsed, isFullScreen]);

    // Import the documentation mapper dynamically
    useEffect(() => {
        const loadDocPath = async () => {
            try {
                const { getToolDocumentationPath } = await import('../../../utils/toolDocumentationMapper');
                const path = getToolDocumentationPath(toolName);
                setDocPath(path);
            } catch (error) {
                console.error('Error loading documentation mapper:', error);
                setDocPath(null);
            }
        };

        if (toolName && isOpen) {
            loadDocPath();
        }
    }, [toolName, isOpen]);

    // Handle escape key to close panel
    useEffect(() => {
        const handleEscape = (event) => {
            if (event.key === 'Escape' && isOpen) {
                if (isFullScreen) {
                    setIsFullScreen(false);
                } else {
                    onClose();
                }
            }
        };

        if (isOpen) {
            document.addEventListener('keydown', handleEscape);
            return () => document.removeEventListener('keydown', handleEscape);
        }
    }, [isOpen, isFullScreen, onClose]);

    // Prevent body scroll when panel is open
    useEffect(() => {
        if (isOpen && !isCollapsed) {
            document.body.style.overflow = 'hidden';
        } else {
            document.body.style.overflow = 'unset';
        }

        return () => {
            document.body.style.overflow = 'unset';
        };
    }, [isOpen, isCollapsed]);

    const handleToggleCollapse = () => {
        setIsCollapsed(!isCollapsed);
        if (isFullScreen) {
            setIsFullScreen(false);
        }
    };

    const handleToggleFullScreen = () => {
        setIsFullScreen(!isFullScreen);
        if (isCollapsed) {
            setIsCollapsed(false);
        }
    };

    const handleClosePanel = () => {
        setIsCollapsed(false);
        setIsFullScreen(false);
        onClose();
    };

    // Don't render if not open or no documentation path
    if (!isOpen || !docPath) {
        return null;
    }

    return (
        <AnimatePresence>
            {isOpen && (
                <>
                    {/* Overlay for fullscreen mode */}
                    {isFullScreen && (
                        <motion.div
                            className="floating-doc-overlay"
                            initial={{ opacity: 0 }}
                            animate={{ opacity: 1 }}
                            exit={{ opacity: 0 }}
                            transition={{ duration: 0.3 }}
                            onClick={() => setIsFullScreen(false)}
                        />
                    )}

                    <motion.div
                        className={`floating-doc-panel ${isCollapsed ? 'collapsed' : ''} ${isFullScreen ? 'fullscreen' : ''}`}
                        initial={{ 
                            x: isFullScreen ? 0 : '100%',
                            y: 0,
                            opacity: 0 
                        }}
                        animate={{ 
                            x: isFullScreen ? 0 : 0,
                            y: isFullScreen ? 0 : 0,
                            opacity: 1 
                        }}
                        exit={{ 
                            x: '100%',
                            opacity: 0 
                        }}
                        transition={{ 
                            type: 'spring',
                            damping: 25,
                            stiffness: 200
                        }}
                        drag={!isFullScreen}
                        dragMomentum={false}
                        dragElastic={0.1}
                        dragConstraints={!isFullScreen ? dragConstraints : false}
                        whileDrag={{ scale: 1.02, cursor: 'grabbing' }}
                        style={{
                            x: isFullScreen ? 0 : undefined,
                            y: isFullScreen ? 0 : undefined,
                            pointerEvents: 'auto'
                        }}
                        onClick={(e) => e.stopPropagation()}
                    >
                        {/* Panel Header */}
                        <div className="floating-doc-header" style={{ cursor: !isFullScreen ? 'grab' : 'default' }}>
                            <div className="floating-doc-title">
                                <IoDocumentText className="floating-doc-icon" />
                                {!isCollapsed && <h3>{toolName}</h3>}
                                {!isFullScreen && (
                                    <div className="floating-doc-drag-indicator">
                                        <div className="drag-dot"></div>
                                        <div className="drag-dot"></div>
                                        <div className="drag-dot"></div>
                                        <div className="drag-dot"></div>
                                        <div className="drag-dot"></div>
                                        <div className="drag-dot"></div>
                                    </div>
                                )}
                            </div>
                            <div className="floating-doc-controls">
                                {!isCollapsed && (
                                    <button
                                        className="floating-doc-control-btn"
                                        onClick={handleToggleFullScreen}
                                        title={isFullScreen ? "Exit fullscreen" : "Enter fullscreen"}
                                    >
                                        {isFullScreen ? <FiMinimize2 /> : <FiMaximize2 />}
                                    </button>
                                )}
                                {!isCollapsed && (
                                    <button
                                        className="floating-doc-control-btn"
                                        onClick={handleToggleCollapse}
                                        title="Collapse documentation"
                                    >
                                        <FiMinimize2 />
                                    </button>
                                )}
                                {isCollapsed && (
                                    <button
                                        className="floating-doc-control-btn"
                                        onClick={handleToggleFullScreen}
                                        title={isFullScreen ? "Exit fullscreen" : "Enter fullscreen"}
                                    >
                                        {isFullScreen ? <FiMinimize2 /> : <FiMaximize2 />}
                                    </button>
                                )}
                                <button
                                    className="floating-doc-control-btn close-btn"
                                    onClick={handleClosePanel}
                                    title="Close documentation"
                                >
                                    <FiX />
                                </button>
                            </div>
                        </div>

                        {/* Panel Content */}
                        <AnimatePresence>
                            {!isCollapsed && (
                                <motion.div
                                    className="floating-doc-content"
                                    initial={{ height: 0, opacity: 0 }}
                                    animate={{ height: 'auto', opacity: 1 }}
                                    exit={{ height: 0, opacity: 0 }}
                                    transition={{ duration: 0.3 }}
                                >
                                    <div className="floating-doc-scroll">
                                        <DocRenderer filePath={docPath} />
                                    </div>
                                </motion.div>
                            )}
                        </AnimatePresence>

                        {/* Collapsed state indicator */}
                        {isCollapsed && (
                            <motion.div
                                className="floating-doc-collapsed-indicator"
                                initial={{ opacity: 0 }}
                                animate={{ opacity: 1 }}
                                exit={{ opacity: 0 }}
                                transition={{ delay: 0.2 }}
                                onClick={handleToggleCollapse}
                            >
                                <h4 className="floating-doc-collapsed-title">{toolName}</h4>
                                <FiBook />
                                <span>Click to expand docs</span>
                            </motion.div>
                        )}
                    </motion.div>
                </>
            )}
        </AnimatePresence>
    );
};

export default FloatingDocPanel;
