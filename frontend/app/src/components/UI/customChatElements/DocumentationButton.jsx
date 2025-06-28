import React from 'react';
import { motion } from 'framer-motion';
import { IoDocumentText } from 'react-icons/io5';
import { hasDocumentation } from '../../../utils/toolDocumentationMapper';

const DocumentationButton = ({ toolName, onShowDocs, style = {} }) => {
    
    // Check if this tool has documentation available
    const hasDocsAvailable = hasDocumentation(toolName);
    
    // Don't render anything if no documentation is available
    if (!hasDocsAvailable) {
        return null;
    }
    
    const handleShowDocs = () => {
        if (onShowDocs) {
            onShowDocs(toolName);
        }
    };
    
    return (
        <motion.div
            className="tool-icon-container"
            onClick={handleShowDocs}
            style={{ 
                cursor: "pointer", 
                marginLeft: "10px",
                ...style 
            }}
            whileHover={{ scale: 1.1 }}
            whileTap={{ scale: 0.95 }}
            title={`View ${toolName} documentation`}
        >
            <IoDocumentText 
                className="chat-tool-icon" 
                style={{ color: "var(--color-info)" }} 
            />
        </motion.div>
    );
};

export default DocumentationButton;
