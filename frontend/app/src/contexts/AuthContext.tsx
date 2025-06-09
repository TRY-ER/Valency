import React, { createContext, useContext, useState, useEffect } from 'react';
import type { ReactNode } from 'react'; // Import ReactNode as a type
import AuthService from '../services/api/AuthService.ts'; // Import AuthService

interface AuthContextType {
  isAuthenticated: boolean;
  isVerified: boolean; // Assuming verification status is part of auth state
  user: any; // Replace 'any' with your User type
  accessToken: string | null;
  isLoading: boolean;
  login: (accessToken: string, refreshToken: string, userData: any, isVerified: boolean) => void;
  logout: () => void;
  setIsVerifiedState: (verified: boolean) => void; // To update verification status post-login if needed
}

const AuthContext = createContext<AuthContextType | undefined>(undefined);

export const useAuthContext = () => {
  const context = useContext(AuthContext);
  if (context === undefined) {
    throw new Error('useAuthContext must be used within an AuthProvider');
  }
  return context;
};

export const AuthProvider: React.FC<{ children: ReactNode }> = ({ children }) => {
  const [user, setUser] = useState<any>(null); // Replace 'any' with your User type
  const [accessToken, setAccessToken] = useState<string | null>(localStorage.getItem('access_token'));
  const [isAuthenticated, setIsAuthenticated] = useState<boolean>(!!localStorage.getItem('access_token'));
  const [isVerified, setIsVerified] = useState<boolean>(localStorage.getItem('isUserVerified') === 'true');
  const [isLoading, setIsLoading] = useState<boolean>(true);

  // Listen for storage changes to handle logout from other tabs or services
  useEffect(() => {
    const handleStorageChange = (e: StorageEvent) => {
      if (e.key === 'access_token' && e.newValue === null) {
        // Token was removed, likely due to logout
        console.log('Auth token removed, updating auth state...');
        setAccessToken(null);
        setUser(null);
        setIsAuthenticated(false);
        setIsVerified(false);
      }
    };

    const handleLogoutEvent = (e: CustomEvent) => {
      console.log('Logout event received:', e.detail);
      // Force logout when triggered by external services
      logout();
    };

    window.addEventListener('storage', handleStorageChange);
    window.addEventListener('auth:logout', handleLogoutEvent as EventListener);
    
    return () => {
      window.removeEventListener('storage', handleStorageChange);
      window.removeEventListener('auth:logout', handleLogoutEvent as EventListener);
    };
  }, []);

  useEffect(() => {
    const initializeAuth = async () => {
      setIsLoading(true);
      let currentAccessToken = localStorage.getItem('access_token');

      if (currentAccessToken) {
        // Attempt to validate the current access token
        let validationResponse = await AuthService.getCurrentUser();

        if (!validationResponse.success) {
          // Access token is invalid or expired, try to refresh it
          console.log("Access token invalid, attempting refresh...");
          const refreshSuccess = await AuthService.refreshToken();

          if (refreshSuccess) {
            console.log("Token refresh successful.");
            // Token refreshed, get the new access token and retry getCurrentUser
            currentAccessToken = AuthService.getAccessToken();
            if (currentAccessToken) {
              // AuthService.storeAuthToken() already updated the auth header
              validationResponse = await AuthService.getCurrentUser();
            } else {
              // Should not happen if refreshToken() was successful and stored the token
              console.error("Failed to get new access token after refresh.");
              logout(); // Logout if new token is not available
              setIsLoading(false);
              return;
            }
          } else {
            // Refresh token failed (e.g., expired or invalid) - logout immediately
            console.log("Refresh token is invalid or expired. Logging out immediately.");
            logout();
            setIsLoading(false);
            return;
          }
        }

        // After validation or successful refresh, check the response
        if (validationResponse.success && validationResponse.data) {
          setUser(validationResponse.data);
          setAccessToken(currentAccessToken); // Set React state with the current (potentially new) token
          setIsAuthenticated(true);
          setIsVerified(validationResponse.data.email_verified);
          // AuthService.storeAuthToken would have updated localStorage if tokens were refreshed
          // If only validation, ensure user data in localStorage is up-to-date
          localStorage.setItem('user', JSON.stringify(validationResponse.data));
          localStorage.setItem('isUserVerified', String(validationResponse.data.email_verified));
        } else {
          // Validation failed even after potential refresh, or getCurrentUser failed for other reasons
          console.error("Failed to validate user after potential refresh:", validationResponse.error);
          logout();
        }
      } else {
        // No initial access token found
        setIsAuthenticated(false);
        setIsVerified(false);
        setUser(null);
        // Ensure all auth-related items are cleared if no token
        AuthService.clearAuthData(); // Clears localStorage and ApiClient header
      }
      setIsLoading(false);
    };

    initializeAuth();
  }, []);

  const login = (newAccessToken: string, newRefreshToken: string, userData: any, verifiedStatus: boolean) => {
    // When logging in, AuthService.storeAuthToken is usually called by the login method itself.
    // Here, we ensure the AuthContext state is also updated.
    AuthService.storeAuthToken({ 
      access_token: newAccessToken, 
      refresh_token: newRefreshToken, 
      user: userData, 
      token_type: 'bearer', // Assuming bearer token type
      mfa_required: userData.mfa_enabled // Assuming mfa_enabled is part of userData
    });
    setAccessToken(newAccessToken);
    setUser(userData);
    setIsAuthenticated(true);
    setIsVerified(verifiedStatus);
  };

  const logout = () => {
    // Clear all authentication data immediately
    AuthService.clearAuthData(); // Clears localStorage and ApiClient auth header
    setAccessToken(null);
    setUser(null);
    setIsAuthenticated(false);
    setIsVerified(false);
  };

  const setIsVerifiedState = (verified: boolean) => {
    localStorage.setItem('isUserVerified', String(verified));
    setIsVerified(verified);
  }

  if (isLoading) {
    return <div>Loading authentication...</div>; 
  }

  return (
    <AuthContext.Provider value={{ isAuthenticated, isVerified, user, accessToken, isLoading, login, logout, setIsVerifiedState }}>
      {children}
    </AuthContext.Provider>
  );
};
