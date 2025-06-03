import React from 'react';
import "./FlexRenderer.css"; // Assuming you have some styles defined in this file
import {
    FunctionResponse,
    FunctionCall,
    TextResponse,
    AgentFunctionTransfer,
    AgentFunctionRecieve,
    ErrorDisplay,      // Added ErrorDisplay import
    CompleteDisplay    // Added CompleteDisplay import
} from './customChatElements/FunctionResponse';

const FlexRenderer = ({ items }) => {
  if (!items || !Array.isArray(items)) {
    return null; // Or some fallback UI
  }

  return (
    <div style={{ display: 'flex', flexDirection: 'column' }} className='flex-renderer-container'>
      {items.map((item, index) => {
        // Placeholder for type-based rendering
        // You can expand this switch statement to handle different item types
        switch (item.type) {
          case 'function_response':
            if (item.data.name === "transfer_to_agent"){
                return <AgentFunctionTransfer key={index} data={item.data} />; 
            }
            return <FunctionResponse key={index} data={item.data} />;
          case 'function_call':
            if (item.data.name === "transfer_to_agent"){
                return <AgentFunctionRecieve key={index} data={item.data.args} />;
            }
            return <FunctionCall key={index} data={item.data} />;
          case 'text':
            console.log("text content recieved >>", item.content)
            return <TextResponse key={index} content={item.content} />;
          case 'error':
            // Use ErrorDisplay for error messages
            return <ErrorDisplay key={index} message={item.data} timestamp={item.timestamp} />;
          case 'complete':
            // Use CompleteDisplay for completion messages
            return <CompleteDisplay key={index} timestamp={item.data} />;
          // Add more cases as needed for different types of components
          default:
            // Handle unknown type or render nothing
            // It might be useful to log a warning for unhandled types
            console.warn(`Unhandled item type: ${item.type}`, item);
            // Fallback for string content if type is not specified or unhandled
            if (typeof item === 'string') {
                return <p key={index}>{item}</p>;
            }
            return null; 
        }
      })}
    </div>
  );
};

export default FlexRenderer;
