import React, { useState, useEffect } from "react";
import "./common.css"
import CustomMarkdownRenderer from "../CustomMarkdown";
import ReactMarkdown from "react-markdown";
import { motion, AnimatePresence } from "framer-motion";
import remarkGfm from "remark-gfm";
import { FiTool, FiUser, FiAlertCircle, FiChevronDown, FiChevronUp, FiHelpCircle, FiMaximize, FiX } from "react-icons/fi"; // Added maximize and close icons
import { BsCheck2All } from 'react-icons/bs';
import { getToolResponse } from "../../../services/api/agentService"; // Import the getToolResponse function
import ToolResponseHandler from "../../ToolResponseHandler/ToolResponseHandler";
import toolTypeMapper from "../../../utils/toolTypeMapper";
import RedirectionButton from "./RedirectionButton";
import DocumentationButton from "./DocumentationButton";
import FloatingDocPanel from "../FloatingDocPanel/FloatingDocPanel";

const FunctionResponse = ({ data }) => {
    const [expanded, setExpanded] = useState(false);
    const [toolResponse, setToolResponse] = useState(null);
    const [loading, setLoading] = useState(false);
    const [error, setError] = useState(null);
    const [toolTransferData, setToolTransferData] = useState(null);
    const [isModalOpen, setIsModalOpen] = useState(false);
    const [isDocPanelOpen, setIsDocPanelOpen] = useState(false);

    // Fetch the tool response when the component mounts or when expanded is toggled to true
    useEffect(() => {
        if (expanded && !toolResponse && data.function_response_id) {
            setLoading(true);
            setError(null);

            getToolResponse(data.function_response_id)
                .then(response => {
                    setToolResponse(response);
                })
                .catch(err => {
                    console.error("Error fetching tool response:", err);
                    setError("Failed to load the complete tool response");
                })
                .finally(() => {
                    setLoading(false);
                });
        }
    }, [expanded, data.id, toolResponse]);

    useEffect(() => {
        console.log("tool response >>", toolResponse)
        if (toolResponse) {
            if (toolResponse.response_data.content){
                let content = toolResponse.response_data.content[0];
                if (content.text){
                    const main_return = JSON.parse(content.text);
                    let rawData;
                    if (main_return.data) {
                        rawData = main_return.data;
                    }
                    else {
                        rawData = main_return;
                    }
                    
                    // Process the tool data using the type mapper
                    const processedData = toolTypeMapper.processToolData(data.name, rawData);
                    setToolTransferData(processedData);
                }
            }
        } 
    }, [toolResponse])

    const toggleExpand = () => setExpanded(!expanded);
    const openModal = () => setIsModalOpen(true);
    const closeModal = () => setIsModalOpen(false);
    const handleShowDocs = (toolName) => {
        setIsDocPanelOpen(true);
    };
    const handleCloseDocs = () => setIsDocPanelOpen(false);

    // Handle Escape key to close modal
    useEffect(() => {
        const handleEscape = (event) => {
            if (event.key === 'Escape' && isModalOpen) {
                closeModal();
            }
        };

        if (isModalOpen) {
            document.addEventListener('keydown', handleEscape);
            // Prevent body scroll when modal is open
            document.body.style.overflow = 'hidden';
        }

        return () => {
            document.removeEventListener('keydown', handleEscape);
            document.body.style.overflow = 'unset';
        };
    }, [isModalOpen]);


    return (
        <motion.div
            initial={{ opacity: 0, y: 20 }}
            animate={{ opacity: 1, y: 0 }}
            transition={{ duration: 0.5 }}
            className="func-tag-container"
        >
            <div className="func-tag-holder" style={{ flexDirection: "column", width: "100%" }}>
                <div style={{ display: "flex", alignItems: "center", width: "100%" }}>
                    <div className="func-tag-inner-container" style={{ marginLeft: "10px" }}>
                        <p>{data.name}</p>
                    </div>
                    <div className="func-tag-inner-container" style={{ marginLeft: "10px" }}>
                        <p>ID: {data?.function_response_id}</p>
                    </div>
                    <div className="tool-icon-container">
                        <FiTool className="chat-tool-icon" />
                        <BsCheck2All className="chat-tool-icon" style={{ marginLeft: "5px", color: "#2ecc71" }} />
                    </div>
                    <RedirectionButton 
                        toolName={data.name} 
                        label={`Go to ${data.name} tool`}
                    />
                    <DocumentationButton 
                        toolName={data.name}
                        onShowDocs={handleShowDocs}
                    />
                    {toolResponse && toolTransferData && (
                        <motion.div
                            className="tool-icon-container"
                            onClick={openModal}
                            style={{ cursor: "pointer", marginLeft: "10px" }}
                            whileHover={{ scale: 1.1 }}
                            whileTap={{ scale: 0.95 }}
                            title="Open in fullscreen"
                        >
                            <FiMaximize className="chat-tool-icon" style={{ color: "var(--color-text-secondary)" }} />
                        </motion.div>
                    )}
                    <motion.div
                        className="tool-icon-container"
                        onClick={toggleExpand}
                        style={{ cursor: "pointer", marginLeft: "auto" }}
                        whileHover={{ scale: 1.1 }}
                        whileTap={{ scale: 0.95 }}
                    >
                        {expanded ? <FiChevronUp className="chat-tool-icon" /> : <FiChevronDown className="chat-tool-icon" />}
                    </motion.div>
                </div>

                <AnimatePresence>
                    {expanded && (
                        <motion.div
                            key="content"
                            initial={{ opacity: 0, height: 0 }}
                            animate={{ opacity: 1, height: "auto" }}
                            exit={{ opacity: 0, height: 0 }}
                            transition={{
                                duration: 0.3,
                                ease: "easeInOut"
                            }}
                            style={{ 
                                marginTop: "10px", 
                                width: "100%", 
                                maxHeight: "70vh",
                                overflowY: "auto",
                                overflowX: "hidden"
                            }}
                        >
                            {loading && (
                                <div className="func-tag-inner-container" style={{ width: "100%" }}>
                                    <p>Loading tool response...</p>
                                </div>
                            )}

                            {error && (
                                <div className="func-tag-inner-container" style={{ width: "100%", color: "red" }}>
                                    <p>{error}</p>
                                </div>
                            )}
                            {toolResponse && 
                                <ToolResponseHandler ToolName={data.name} ToolData={toolTransferData} />
                            }
                        </motion.div>
                    )}
                </AnimatePresence>
            </div>

            {/* Fullscreen Modal */}
            <AnimatePresence>
                {isModalOpen && (
                    <motion.div
                        className="fullscreen-modal-overlay"
                        initial={{ opacity: 0 }}
                        animate={{ opacity: 1 }}
                        exit={{ opacity: 0 }}
                        transition={{ duration: 0.3 }}
                        onClick={closeModal}
                    >
                        <motion.div
                            className="fullscreen-modal-content"
                            initial={{ scale: 0.9, opacity: 0 }}
                            animate={{ scale: 1, opacity: 1 }}
                            exit={{ scale: 0.9, opacity: 0 }}
                            transition={{ duration: 0.3 }}
                            onClick={(e) => e.stopPropagation()}
                        >
                            <div className="fullscreen-modal-header">
                                <h2>Tool Response: {data.name}</h2>
                                <button 
                                    className="fullscreen-modal-close"
                                    onClick={closeModal}
                                    aria-label="Close modal"
                                >
                                    <FiX />
                                </button>
                            </div>
                            <div className="fullscreen-modal-body">
                                {toolResponse && toolTransferData && (
                                    <ToolResponseHandler ToolName={data.name} ToolData={toolTransferData} />
                                )}
                            </div>
                        </motion.div>
                    </motion.div>
                )}
            </AnimatePresence>

            {/* Floating Documentation Panel */}
            <FloatingDocPanel 
                toolName={data.name}
                isOpen={isDocPanelOpen}
                onClose={handleCloseDocs}
            />
        </motion.div>
    )

}

const FunctionCall = ({ data }) => {
    const [expanded, setExpanded] = useState(false);
    const toggleExpand = () => setExpanded(!expanded);

    return (
        <motion.div
            initial={{ opacity: 0, y: 20 }}
            animate={{ opacity: 1, y: 0 }}
            transition={{ duration: 0.5 }}
            className="func-tag-container"
        >
            <div className="func-tag-holder" style={{ flexDirection: "column", width: "100%" }}>
                <div style={{ display: "flex", alignItems: "center", width: "100%" }}>
                    <div className="func-tag-inner-container" style={{ marginLeft: "10px" }}>
                        <p>{data.name}</p>
                    </div>
                    <div className="func-tag-inner-container" style={{ marginLeft: "10px" }}>
                        <p>{data.args ? `${Object.keys(data.args).length} arguments` : 'No arguments'}</p>
                    </div>
                    <div className="tool-icon-container">
                        <FiTool className="chat-tool-icon" />
                        <FiHelpCircle className="chat-tool-icon" style={{ marginLeft: "5px", color: "#2ecc71" }} />
                    </div>
                    <motion.div
                        className="tool-icon-container"
                        onClick={toggleExpand}
                        style={{ cursor: "pointer", marginLeft: "auto" }}
                        whileHover={{ scale: 1.1 }}
                        whileTap={{ scale: 0.95 }}
                    >
                        {expanded ? <FiChevronUp className="chat-tool-icon" /> : <FiChevronDown className="chat-tool-icon" />}
                    </motion.div>
                </div>

                <AnimatePresence>
                    {expanded && data.args && (
                        <motion.div
                            key="content"
                            initial={{ opacity: 0, height: 0 }}
                            animate={{ opacity: 1, height: "auto" }}
                            exit={{ opacity: 0, height: 0 }}
                            transition={{
                                duration: 0.3,
                                ease: "easeInOut"
                            }}
                            style={{ marginTop: "10px", width: "100%", overflow: "hidden" }}
                        >
                            <div className="func-tag-inner-container" style={{ width: "100%" }}>
                                <pre style={{
                                    fontFamily: 'monospace',
                                    padding: '1em',
                                    borderRadius: '3px',
                                    overflowX: 'auto',
                                    margin: 0,
                                }}>
                                    {JSON.stringify(data.args, null, 2)}
                                </pre>
                            </div>
                        </motion.div>
                    )}
                </AnimatePresence>
            </div>
        </motion.div>
    )

}

const TextResponse = ({ content }) => {
    const markdownInput = Array.isArray(content)
        ? content.join('')
        : (typeof content === 'string' ? content : '');

    return (
        <motion.div
            initial={{ opacity: 0, y: 20 }}
            animate={{ opacity: 1, y: 0 }}
            transition={{ duration: 0.5 }}
            className="func-tag-container"
        >
            <div className="func-tag-inner-container" style={{
                wordWrap: 'break-word',
                wordBreak: 'break-all',
                overflowWrap: 'break-word',
                maxWidth: '100%',
                overflow: 'hidden'
            }}>
                <ReactMarkdown
                    remarkPlugins={[remarkGfm]}
                    components={{
                        p: props => <p style={{ 
                            marginBlockStart: '1em', 
                            marginBlockEnd: '1em',
                            wordWrap: 'break-word',
                            wordBreak: 'break-word',
                            overflowWrap: 'break-word'
                        }} {...props} />,
                        strong: props => <strong style={{ fontWeight: 'bold' }} {...props} />,
                        em: props => <em style={{ fontStyle: 'italic' }} {...props} />,
                        code: ({ inline, className, children, ...props }) => {
                            const style = {
                                fontFamily: 'monospace',
                                backgroundColor: 'rgba(211, 211, 211, 0.2)', // Light gray background
                                padding: '0.2em 0.4em',
                                borderRadius: '3px',
                                color: 'rgb(69, 199, 37)', // A common color for code text
                                wordWrap: 'break-word',
                                wordBreak: 'break-all',
                                overflowWrap: 'break-word',
                                whiteSpace: inline ? 'pre-wrap' : 'pre-wrap'
                            };
                            if (inline) {
                                return <code style={style} className={className} {...props}>{children}</code>;
                            }
                            // For block code, ReactMarkdown with remark-gfm usually wraps <code> in <pre>
                            // We apply similar style to the code part, pre will handle block display
                            return <code style={{ ...style, padding: '0' /* Pre will have padding */ }} className={className} {...props}>{children}</code>;
                        },
                        pre: props => <pre style={{
                            fontFamily: 'monospace',
                            backgroundColor: 'rgba(211, 211, 211, 0.4)',
                            padding: '1em',
                            borderRadius: '3px',
                            overflowX: 'auto',
                            color: 'rgb(199, 37, 78)',
                            wordWrap: 'break-word',
                            wordBreak: 'break-all',
                            overflowWrap: 'break-word',
                            whiteSpace: 'pre-wrap',
                            maxWidth: '100%'
                        }} {...props} />,
                        ul: props => <ul style={{ 
                            listStyleType: 'disc', 
                            paddingInlineStart: '40px', 
                            marginBlockStart: '1em', 
                            marginBlockEnd: '1em',
                            wordWrap: 'break-word',
                            wordBreak: 'break-word',
                            overflowWrap: 'break-word'
                        }} {...props} />,
                        ol: props => <ol style={{ 
                            listStyleType: 'decimal', 
                            paddingInlineStart: '40px', 
                            marginBlockStart: '1em', 
                            marginBlockEnd: '1em',
                            wordWrap: 'break-word',
                            wordBreak: 'break-word',
                            overflowWrap: 'break-word'
                        }} {...props} />,
                        li: props => <li style={{ 
                            display: 'list-item',
                            wordWrap: 'break-word',
                            wordBreak: 'break-word',
                            overflowWrap: 'break-word'
                        }} {...props} />,
                        // You can add more elements here like h1, etc. if needed
                        // h1: props => <h1 style={{ fontSize: '2em', fontWeight: 'bold' }} {...props} />,
                    }}
                >
                    {markdownInput}
                </ReactMarkdown>
            </div>
        </motion.div>
    );
}

const AgentFunctionTransfer = ({ data }) => {
    return (
        <motion.div
            initial={{ opacity: 0, y: 20 }}
            animate={{ opacity: 1, y: 0 }}
            transition={{ duration: 0.5 }}
            className="func-tag-container"
        >
            <div className="func-tag-holder">
                <div className="func-tag-inner-container">
                    <p>Transfering to Agent</p>
                </div>
                <div className="tool-icon-container">
                    <FiUser className="chat-tool-icon" /> {/* Added FiUser icon */}
                </div>
            </div>
        </motion.div>
    );
}

const AgentFunctionRecieve = ({ data }) => {
    return (
        <motion.div
            initial={{ opacity: 0, y: 20 }}
            animate={{ opacity: 1, y: 0 }}
            transition={{ duration: 0.5 }}
            className="func-tag-container"
        >
            <div className="func-tag-holder">
                <div className="func-tag-inner-container">
                    <p>{data?.agent_name}</p>
                </div>
                <div className="tool-icon-container">
                    <FiUser className="chat-tool-icon" /> {/* Added FiUser icon */}
                </div>
            </div>
        </motion.div>
    );
}

const ErrorDisplay = ({ message, timestamp }) => {
    const options = {
        day: '2-digit',
        month: '2-digit',
        year: 'numeric',
        hour: '2-digit',
        minute: '2-digit',
        timeZoneName: 'short'
    };

    let timeString = '';
    if (timestamp) {
        try {
            const formattedDate = new Date(timestamp).toLocaleString(undefined, options);
            timeString = ` at: ${formattedDate}`;
        } catch (e) {
            console.error("Error formatting timestamp:", e);
            timeString = ' (invalid time)';
        }
    }

    return (
        <motion.div
            initial={{ opacity: 0, y: 20 }}
            animate={{ opacity: 1, y: 0 }}
            transition={{ duration: 0.5 }}
            className="func-tag-container"
        >
            <div className="func-tag-inner-container" style={{ display: 'flex', alignItems: 'center', color: 'red' }}>
                <FiAlertCircle style={{ marginRight: '8px', fontSize: '1.2em' }} />
                <span>Error: {message}{timeString}</span>
            </div>
        </motion.div>
    );
};

const CompleteDisplay = ({ timestamp }) => {
    const options = {
        day: '2-digit',
        month: '2-digit',
        year: 'numeric',
        hour: '2-digit',
        minute: '2-digit',
        timeZoneName: 'short'
    };
    const formattedDate = timestamp ? new Date(timestamp).toLocaleString(undefined, options) : 'N/A';

    return (
        <motion.div
            initial={{ opacity: 0, y: 20 }}
            animate={{ opacity: 1, y: 0 }}
            transition={{ duration: 0.5 }}
            className="func-tag-container"
        >
            <div className="func-tag-inner-container" style={{ display: 'flex', alignItems: 'center', color: 'green' }}>
                <BsCheck2All style={{ marginRight: '8px', fontSize: '1.2em' }} />
                <span>Stream completed{formattedDate !== 'N/A' ? ` at: ${formattedDate}` : ''}</span>
            </div>
        </motion.div>
    );
};

export {
    FunctionResponse,
    FunctionCall,
    TextResponse,
    AgentFunctionTransfer,
    AgentFunctionRecieve,
    ErrorDisplay,
    CompleteDisplay
};