import React from 'react';
import './SpinLoader.css';


const SpinningLoader = ({isDark}) => {
  console.log("spin is dark>>",isDark)
  return (
  <div className={`loader-container ${isDark ? "dark" : ""}`}>
      <div className="loader">
        <div className="loader-dot"></div>
        <div className="loader-dot"></div>
        <div className="loader-dot"></div>
        <div className="loader-dot"></div>
      </div>
    </div>
  );
};

export default SpinningLoader;