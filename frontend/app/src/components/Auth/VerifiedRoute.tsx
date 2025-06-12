import React from 'react';
import { Navigate, Outlet, useLocation } from 'react-router-dom';
import { useAuthContext } from '../../contexts/AuthContext.tsx';

interface VerifiedRouteProps {
  // You can add any specific props for VerifiedRoute here if needed
}

const VerifiedRoute: React.FC<VerifiedRouteProps> = () => {
  const { isAuthenticated, isVerified, isLoading } = useAuthContext();
  const location = useLocation();

  if (isLoading) {
    // Show a loading indicator while auth state is being determined
    return <div>Loading...</div>; 
  }

  if (!isAuthenticated) {
    // User not authenticated, redirect to login page
    // Pass the current location to redirect back after successful login
    return <Navigate to="/login" state={{ from: location }} replace />;
  }

  if (!isVerified) {
    // User is authenticated but email is not verified, redirect to verify email page
    return <Navigate to="/verify-email" replace />;
  }

  // User is authenticated and verified, render the child components
  return <Outlet />;
};

export default VerifiedRoute;
