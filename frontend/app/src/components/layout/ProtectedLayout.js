import React from 'react';
import { Outlet } from 'react-router-dom';
import Sidebar from '../sidebar/sidebar';
import '../../App.css'; // Assuming App.css contains the layout styles for .App, .sidebar-content, .side-content

const ProtectedLayout = () => {
  return (
    <div className="App"> {/* This div handles the overall flex layout */}
      <div className="sidebar-content">
        <Sidebar />
      </div>
      <div className="side-content">
        <Outlet /> {/* Protected page content will be rendered here */}
      </div>
    </div>
  );
};

export default ProtectedLayout;