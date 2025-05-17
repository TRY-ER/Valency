import React from 'react';
import { Navigate, Outlet } from 'react-router-dom';
import { useAuthContext } from '../../contexts/AuthContext.tsx';

const AuthRedirectRoute: React.FC = () => {
  const { isAuthenticated, isLoading } = useAuthContext();

  if (isLoading) {
    // Show a loading indicator while auth state is being determined
    return <div>Loading authentication status...</div>; 
  }

  if (isAuthenticated) {
    // If user is authenticated, redirect them to the main page (e.g., root)
    return <Navigate to="/" replace />;
  }

  // If user is not authenticated, render the requested auth page (Login, Signup, etc.)
  return <Outlet />;
};

export default AuthRedirectRoute;
