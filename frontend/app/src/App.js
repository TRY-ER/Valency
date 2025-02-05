import logo from './logo.svg';
import './App.css';
import React, { useEffect, useState } from 'react';
import { BrowserRouter as Router } from 'react-router-dom';
import { useLocation } from 'react-router-dom';
import ReRoutes from './routes';
import Sidebar from './components/sidebar/sidebar';
import menuContent from './contents/menuContent';

function App() {
  const location = useLocation();
  const [activeMenu, setActiveMenu] = useState(() => {
    const currentPath = location.pathname;
    const menuItem = menuContent.find(item =>
      currentPath.startsWith(item.link) && item.link !== '/'
    ) || menuContent.find(item => item.link === '/');
    return menuItem ? menuItem.id : 1;
  });

  useEffect(() => {
    const currentPath = location.pathname;
    const menuItem = menuContent.find(item =>
      currentPath.startsWith(item.link) && item.link !== '/'
    ) || menuContent.find(item => item.link === '/');
    if (menuItem) {
      setActiveMenu(menuItem.id);
    }
  }, [location]);


  return (
    <div className="App">
      <div className="sidebar-content">
        <Sidebar activeMenu={activeMenu} setActiveMenu={setActiveMenu} />
      </div>
      <div className="side-content">
        <ReRoutes />
      </div>
    </div>
  );
}

export default App;
