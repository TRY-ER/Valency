import { useEffect } from 'react';
import { useAuthContext } from '../contexts/AuthContext.tsx';
import { createAuthEventListener, checkAuthStatus } from '../utils/authUtils.ts';

/**
 * Custom hook for handling authentication in components
 * Provides auth state and handles automatic logout on token expiration
 */
export const useAuth = () => {
  const authContext = useAuthContext();

  useEffect(() => {
    // Create auth event listener
    const cleanup = createAuthEventListener(
      (detail) => {
        console.log('useAuth: Logout event received:', detail);
        // The AuthContext already handles logout via storage events
        // But we can add additional cleanup here if needed
      },
      () => {
        console.log('useAuth: Token refresh event received');
        // Token was refreshed successfully
        // The AuthContext will handle updating the state
      }
    );

    return cleanup;
  }, []);

  return {
    ...authContext,
    // Add utility methods
    checkStatus: checkAuthStatus,
  };
};

/**
 * Hook for components that need to handle auth events without accessing full auth context
 */
export const useAuthEvents = (options: {
  onLogout?: (detail: any) => void;
  onTokenRefresh?: () => void;
}) => {
  useEffect(() => {
    const cleanup = createAuthEventListener(
      options.onLogout,
      options.onTokenRefresh
    );

    return cleanup;
  }, [options.onLogout, options.onTokenRefresh]);
};
