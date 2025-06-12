import React from 'react';
import { Navigate, Outlet } from 'react-router-dom';
import { useAuthContext } from '../../contexts/AuthContext.tsx';

const AuthRedirectRoute: React.FC = () => {
  const { isAuthenticated, isVerified, isLoading } = useAuthContext();

  if (isLoading) {
    // Show a loading indicator while auth state is being determined
    return <div>Loading authentication status...</div>; 
  }

  if (isAuthenticated) {
    // If user is authenticated but not verified, redirect to verify email page
    if (!isVerified) {
      return <Navigate to="/verify-email" replace />;
    }
    
    // If user is authenticated and verified, redirect them to the main page
    return <Navigate to="/" replace />;
  }

  // If user is not authenticated, render the requested auth page (Login, Signup, etc.)
  return <Outlet />;
};

export default AuthRedirectRoute;
