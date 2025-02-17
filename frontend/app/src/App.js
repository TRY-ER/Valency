import logo from './logo.svg';
import './App.css';
import React, { createContext } from 'react';
import ReRoutes from './routes';
import Sidebar from './components/sidebar/sidebar';
import { MasterToolProvider } from './contexts/MasterToolContexts';

function App() {
  return (
    <MasterToolProvider>
      <div className="App">
        <div className="sidebar-content">
          <Sidebar />
        </div>
        <div className="side-content">
          <ReRoutes />
        </div>
      </div>
    </MasterToolProvider>
  );
}

export default App;
