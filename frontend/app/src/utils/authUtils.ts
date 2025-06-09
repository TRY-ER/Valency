/**
 * Utility functions for authentication management across the application
 */

import AuthService from '../services/api/AuthService.ts';

/**
 * Triggers a global logout event that all services can listen to
 * @param reason - The reason for logout (e.g., 'token_expired', 'user_action', 'security')
 * @param source - The source that triggered the logout (e.g., 'mcp_service', 'user_click')
 */
export const triggerGlobalLogout = (reason: string = 'unknown', source: string = 'unknown'): void => {
  console.log(`Global logout triggered - Reason: ${reason}, Source: ${source}`);
  
  // Clear auth data first
  AuthService.clearAuthData();
  
  // Dispatch custom event to notify all components
  if (typeof window !== 'undefined') {
    window.dispatchEvent(new CustomEvent('auth:logout', { 
      detail: { reason, source, timestamp: new Date().toISOString() }
    }));
  }
};

/**
 * Checks if the current access token is still valid without making a network request
 * This is a basic check based on localStorage presence
 */
export const isTokenPresent = (): boolean => {
  return !!AuthService.getAccessToken();
};

/**
 * Validates the current token by making a request to the server
 * @returns Promise<boolean> - true if token is valid, false otherwise
 */
export const validateCurrentToken = async (): Promise<boolean> => {
  try {
    const response = await AuthService.getCurrentUser();
    return response.success;
  } catch (error) {
    console.error('Token validation failed:', error);
    return false;
  }
};

/**
 * Attempts to refresh the token if it's expired
 * @returns Promise<boolean> - true if refresh was successful, false otherwise
 */
export const attemptTokenRefresh = async (): Promise<boolean> => {
  try {
    return await AuthService.refreshToken();
  } catch (error) {
    console.error('Token refresh failed:', error);
    return false;
  }
};

/**
 * Creates a wrapper function that automatically handles token refresh for API calls
 * @param apiCall - The API function to wrap
 * @returns A wrapped function that handles token refresh automatically
 */
export const withTokenRefresh = <T extends (...args: any[]) => Promise<any>>(
  apiCall: T
): T => {
  return (async (...args: any[]) => {
    try {
      // First attempt
      return await apiCall(...args);
    } catch (error: any) {
      // If it's a 401 error, try refreshing the token
      if (error.response?.status === 401 || error.status === 401) {
        console.log('API call failed with 401, attempting token refresh...');
        
        const refreshSuccess = await attemptTokenRefresh();
        if (refreshSuccess) {
          console.log('Token refresh successful, retrying API call...');
          return await apiCall(...args);
        } else {
          console.log('Token refresh failed, triggering logout...');
          triggerGlobalLogout('token_refresh_failed', 'api_wrapper');
          throw error;
        }
      }
      throw error;
    }
  }) as T;
};

/**
 * Custom hook utility for components that need to handle auth state changes
 * This can be used in React components to listen for logout events
 */
export const createAuthEventListener = (
  onLogout?: (detail: any) => void,
  onTokenRefresh?: () => void
) => {
  const handleLogout = (event: CustomEvent) => {
    console.log('Auth logout event received:', event.detail);
    if (onLogout) {
      onLogout(event.detail);
    }
  };

  const handleTokenRefresh = () => {
    console.log('Token refresh event received');
    if (onTokenRefresh) {
      onTokenRefresh();
    }
  };

  // Add event listeners
  if (typeof window !== 'undefined') {
    window.addEventListener('auth:logout', handleLogout as EventListener);
    window.addEventListener('auth:token_refreshed', handleTokenRefresh as EventListener);
  }

  // Return cleanup function
  return () => {
    if (typeof window !== 'undefined') {
      window.removeEventListener('auth:logout', handleLogout as EventListener);
      window.removeEventListener('auth:token_refreshed', handleTokenRefresh as EventListener);
    }
  };
};

/**
 * Utility to check if a user is authenticated and verified
 */
export const checkAuthStatus = () => {
  const user = AuthService.getUser();
  const isLoggedIn = AuthService.isLoggedIn();
  
  return {
    isAuthenticated: isLoggedIn,
    isVerified: user?.email_verified || false,
    isMfaRequired: AuthService.isMfaRequired(),
    user: user
  };
};
