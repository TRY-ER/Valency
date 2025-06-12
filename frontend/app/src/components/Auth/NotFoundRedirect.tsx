import React from 'react';
import { Navigate } from 'react-router-dom';
import { useAuthContext } from '../../contexts/AuthContext.tsx';

const NotFoundRedirect: React.FC = () => {
  const { isAuthenticated, isVerified } = useAuthContext();

  if (isAuthenticated) {
    if (!isVerified) {
      return <Navigate to="/verify-email" replace />;
    }
    return <Navigate to="/" replace />;
  } else {
    return <Navigate to="/login" replace />;
  }
};

export default NotFoundRedirect;
