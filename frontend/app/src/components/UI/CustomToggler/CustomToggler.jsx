import React, { useState } from 'react';
import './CustomToggler.css'; // Import your custom CSS for styling

function CustomToggler({
    feature_name,
    isEnabled,
    setIsEnabled,
    isDark
}) {
  const handleToggle = () => {
    setIsEnabled(!isEnabled);
  };

  return (
    <div className={`custom-toggler ${isDark ? "dark" : ""}`}>
        <div
          className={`toggler-slider ${isDark ? "dark" : ""} ${isEnabled ? 'enabled' : 'disabled'} `}
          onClick={handleToggle}
        >
          <div className={`toggler-thumb ${isDark ? "dark" : ""} ${isEnabled ? 'enabled' : 'disabled'}`} />
        </div>
      <div className="toggle-instruction-text">{isEnabled ? `Press Ctrl+Shift+z while typing to disable ${feature_name}` : 
      `Press Ctrl+Shift+z while typing to enable ${feature_name}`}</div>
    </div>
  );
}

export default CustomToggler;