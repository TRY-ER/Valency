import './App.css';
import React from 'react';
import ReRoutes from './routes';
import { MasterToolProvider } from './contexts/MasterToolContexts';
import { AuthProvider } from './contexts/AuthContext.tsx';
import { ThemeProvider } from './contexts/ThemeContext.js';

function App() {
  return (
    <AuthProvider>
      <ThemeProvider>
        <MasterToolProvider>
          <ReRoutes />
        </MasterToolProvider>
      </ThemeProvider>
    </AuthProvider>
  );
}

export default App;
