import logo from './logo.svg';
import './App.css';
import React, { useEffect, useState } from 'react';
import { BrowserRouter as Router } from 'react-router-dom';
import { useLocation, Outlet } from 'react-router-dom';
import ReRoutes from './routes';
import Sidebar from './components/sidebar/sidebar';
import menuContent from './contents/menuContent';

function App() {
  return (
    <div className="App">
      <div className="sidebar-content">
        <Sidebar />
      </div>
      <div className="side-content">
        <ReRoutes />
      </div>
    </div>
  );
}

export default App;
