import React from "react";
import { getToolComponent, isToolAvailable } from "./ToolMapper";

const ToolResponseHandler = ({ ToolName, ToolData }) => {
    console.log("ToolName >>", ToolName);
    console.log("ToolData >>", ToolData);
    
    let returnable = null;
    
    // Get the appropriate component from ToolMapper
    const ToolComponent = getToolComponent(ToolName, ToolData);
    console.log("tool component >>", ToolComponent);
    
    if (ToolComponent && isToolAvailable(ToolName)) {
        // Render the mapped component with tool data as props
        returnable = (
            <div className="tool-interface-container" style={{ width: "100%" }}>
                {/* <div className="tool-header" style={{ 
                    marginBottom: '10px',
                    padding: '10px',
                    backgroundColor: '#f8f9fa',
                    borderRadius: '5px',
                    borderLeft: '4px solid #007bff'
                }}>
                    <h3 style={{ margin: 0, color: '#333' }}>
                        {ToolName.replace(/([A-Z])/g, ' $1').replace(/^./, str => str.toUpperCase())}
                    </h3>
                </div> */}
                <div className="tool-content">
                     <ToolComponent toolData={ToolData} />
                </div>
            </div>
        );
    } else {
        // Fallback to JSON display for unknown tools
        returnable = (
            <div className="func-tag-inner-container" style={{ width: "100%" }}>
                <div className="tool-header" style={{ 
                    marginBottom: '10px',
                    padding: '10px',
                    backgroundColor: '#fff3cd',
                    borderRadius: '5px',
                    borderLeft: '4px solid #ffc107'
                }}>
                    <h3 style={{ margin: 0, color: '#856404' }}>
                        Unknown Tool: {ToolName}
                    </h3>
                    <p style={{ margin: '5px 0 0 0', fontSize: '14px', color: '#856404' }}>
                        Displaying raw data. Component mapping not found.
                    </p>
                </div>
                <pre style={{
                    fontFamily: 'monospace',
                    padding: '1em',
                    borderRadius: '3px',
                    overflowX: 'auto',
                    margin: 0,
                    border: '1px solid #dee2e6'
                }}>
                    {JSON.stringify(ToolData, null, 2)}
                </pre>
            </div>
        );
    }

    return (
        <>
            {returnable}
        </>
    );
}

export default ToolResponseHandler;