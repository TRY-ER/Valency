import React from 'react';
import { motion } from 'framer-motion';
import { FiExternalLink } from 'react-icons/fi';
import { getRedirectionPath, hasRedirection } from '../../../utils/redirectionMapper.js';

const RedirectionButton = ({ toolName, label, style = {} }) => {
    
    // Check if this tool has a redirection path
    const redirectPath = getRedirectionPath(toolName);
    
    // Don't render anything if no redirection is available
    if (!hasRedirection(toolName)) {
        return null;
    }
    
    const handleRedirection = () => {
        if (redirectPath) {
            // Open in new tab
            window.open(redirectPath, '_blank', 'noopener,noreferrer');
        }
    };
    
    return (
        <motion.div
            className="tool-icon-container"
            onClick={handleRedirection}
            style={{ 
                cursor: "pointer", 
                marginLeft: "10px",
                ...style 
            }}
            whileHover={{ scale: 1.1 }}
            whileTap={{ scale: 0.95 }}
            title={label || `Go to ${toolName} tool`}
        >
            <FiExternalLink 
                className="chat-tool-icon" 
                style={{ color: "var(--color-success)" }} 
            />
        </motion.div>
    );
};

export default RedirectionButton;
