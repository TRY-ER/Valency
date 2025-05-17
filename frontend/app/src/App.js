import './App.css';
import React from 'react';
import ReRoutes from './routes';
import { MasterToolProvider } from './contexts/MasterToolContexts';
import { AuthProvider } from './contexts/AuthContext.tsx';

function App() {
  return (
    <AuthProvider>
      <MasterToolProvider>
        <ReRoutes />
      </MasterToolProvider>
    </AuthProvider>
  );
}

export default App;
