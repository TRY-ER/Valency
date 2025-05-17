// Export all API services for easy imports

import ApiClient from './ApiClient';
import AuthService from './AuthService';
import EndpointManager from './EndpointManager';

// Re-export types for convenience
export type { ApiResponse } from './ApiClient';
export type {
  LoginRequest,
  OtpLoginRequest,
  SignupRequest,
  ForgotPasswordRequest,
  ResetPasswordRequest,
  UserData,
  AuthResponse
} from './AuthService';

// Export service instances
export {
  ApiClient,
  AuthService,
  EndpointManager
};

// Default export for the main manager
export default EndpointManager;