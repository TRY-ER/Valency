import React from "react";
import "./common.css"
import CustomMarkdownRenderer from "../CustomMarkdown";
import ReactMarkdown from "react-markdown";
import { motion } from "framer-motion";
import remarkGfm from "remark-gfm";
import { FiTool, FiUser, FiAlertCircle } from "react-icons/fi"; // Import FiAlertCircle
import { BsCheck2All } from 'react-icons/bs'; // Import BsCheck2All

const FunctionResponse = ({ data }) => {
    return (
        <motion.div
            initial={{ opacity: 0, y: 20 }}
            animate={{ opacity: 1, y: 0 }}
            transition={{ duration: 0.5 }}
            className="func-tag-container"
        >
            <div className="func-tag-holder">
                <div className="func-tag-inner-container">
                    <p>{data.name}</p>
                </div>
                <div className="tool-icon-container">
                    <FiTool className="tool-icon" />
                </div>
            </div>
        </motion.div>
    )

}

const FunctionCall = ({ data }) => {
    return (
        <motion.div
            initial={{ opacity: 0, y: 20 }}
            animate={{ opacity: 1, y: 0 }}
            transition={{ duration: 0.5 }}
            className="func-tag-container"
        >
            <div className="func-tag-holder">
                <div className="func-tag-inner-container">
                    <p>{data.name}</p>
                </div>
                <div className="tool-icon-container">
                    <FiTool className="tool-icon" />
                </div>
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
            <div className="func-tag-inner-container">
                <ReactMarkdown
                    remarkPlugins={[remarkGfm]}
                    components={{
                        p: props => <p style={{ marginBlockStart: '1em', marginBlockEnd: '1em' }} {...props} />,
                        strong: props => <strong style={{ fontWeight: 'bold' }} {...props} />,
                        em: props => <em style={{ fontStyle: 'italic' }} {...props} />,
                        code: ({ inline, className, children, ...props }) => {
                            const style = {
                                fontFamily: 'monospace',
                                backgroundColor: 'rgba(211, 211, 211, 0.4)', // Light gray background
                                padding: '0.2em 0.4em',
                                borderRadius: '3px',
                                color: 'rgb(199, 37, 78)' // A common color for code text
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
                            color: 'rgb(199, 37, 78)'
                        }} {...props} />,
                        ul: props => <ul style={{ listStyleType: 'disc', paddingInlineStart: '40px', marginBlockStart: '1em', marginBlockEnd: '1em' }} {...props} />,
                        ol: props => <ol style={{ listStyleType: 'decimal', paddingInlineStart: '40px', marginBlockStart: '1em', marginBlockEnd: '1em' }} {...props} />,
                        li: props => <li style={{ display: 'list-item' }} {...props} />,
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
                    <FiUser className="tool-icon" /> {/* Added FiUser icon */}
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
                    <FiUser className="tool-icon" /> {/* Added FiUser icon */}
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