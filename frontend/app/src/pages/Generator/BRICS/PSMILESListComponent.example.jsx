// Example usage of PSMILESListComponent with toolData

import React from 'react';
import PSMILESListComponent from './PSMILESListComponent';

// Example toolData that matches the API response schema
const exampleToolData = {
    result: {
        "0": JSON.stringify({
            candidates: [
                "CCO",
                "CC(C)O", 
                "CCCO",
                "CC(C)(C)O",
                "CCCCO"
            ]
        })
    }
};

// Usage Examples:

// 1. Normal usage without pre-loaded data
const NormalUsage = () => {
    return <PSMILESListComponent />;
};

// 2. Usage with pre-loaded data (toolData)
const PreLoadedUsage = () => {
    return <PSMILESListComponent toolData={exampleToolData} />;
};

// 3. Usage with conditional toolData (e.g., from a parent component)
const ConditionalUsage = ({ someToolData }) => {
    return <PSMILESListComponent toolData={someToolData} />;
};

export { NormalUsage, PreLoadedUsage, ConditionalUsage };
