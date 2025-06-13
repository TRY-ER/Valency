import React from 'react';
import { Navigate, Outlet } from 'react-router-dom';
import { useAuthContext } from '../../contexts/AuthContext.tsx';

const AuthRedirectRoute: React.FC = () => {
  const { isAuthenticated, isVerified, isLoading } = useAuthContext();

  if (isLoading) {
    // Show a loading indicator while auth state is being determined
    return <div>Loading authentication status...</div>; 
  }

  // If user is authenticated AND verified, they should not be on auth pages like /login, /signup.
  // So, redirect them to the main application page (e.g., '/').
  if (isAuthenticated && isVerified) {
    return <Navigate to="/" replace />;
  }

  // In all other cases:
  // 1. User is not authenticated (isAuthenticated is false) -> allow access to auth pages.
  // 2. User is authenticated but not verified (isAuthenticated is true, isVerified is false) -> allow access to auth pages.
  // This allows an unverified user to go back to /login or /signup if needed.
  return <Outlet />;
};

export default AuthRedirectRoute;
