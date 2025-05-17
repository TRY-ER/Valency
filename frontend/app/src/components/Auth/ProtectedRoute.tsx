import React from 'react';
import { Navigate, Outlet, useLocation } from 'react-router-dom';
import { useAuthContext } from '../../contexts/AuthContext.tsx';

interface ProtectedRouteProps {
  // You can add any specific props for ProtectedRoute here if needed
}

const ProtectedRoute: React.FC<ProtectedRouteProps> = () => {
  const { isAuthenticated, isLoading } = useAuthContext();
  const location = useLocation();

  if (isLoading) {
    // Optional: Show a loading indicator while auth state is being determined
    return <div>Loading...</div>; 
  }

  if (!isAuthenticated) {
    // User not authenticated, redirect to login page
    // Pass the current location to redirect back after successful login
    return <Navigate to="/login" state={{ from: location }} replace />;
  }

  // User is authenticated, render the child components
  return <Outlet />;
};

export default ProtectedRoute;
