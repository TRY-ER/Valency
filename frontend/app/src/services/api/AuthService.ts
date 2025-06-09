import type { ApiResponse } from './ApiClient.ts';
import ApiClient from './ApiClient.ts';

// Define types for auth-related requests and responses
export interface LoginRequest {
  email: string;
  password: string;
}

export interface OtpLoginRequest {
  email: string;
  otp: string;
}

export interface SignupRequest {
  email: string;
  password: string;
  username?: string;
}

export interface ForgotPasswordRequest {
  email: string;
}

export interface ResetPasswordRequest {
  token: string;
  password: string;
}

export interface UserData {
  id: string;
  email: string;
  username: string; // Ensure username is part of UserData
  mfa_enabled: boolean;
  email_verified: boolean;
  created_at: string; // Assuming ISO string format
  updated_at: string; // Assuming ISO string format
}

export interface AuthResponse {
  access_token: string;
  refresh_token: string;
  token_type: string;
  mfa_required: boolean; // This field might not be consistently present in all auth responses (e.g. refresh)
  user: UserData;
}

// Add RefreshTokenRequestData interface
interface RefreshTokenRequestData {
  refresh_token: string;
}

export interface SignUpResponse {
  id: string;
  email_verified: boolean;
  mfa_enabled: boolean;
  email: string;
  username: string;
}

class AuthService {
  private client: ApiClient;
  private static instance: AuthService;
  private isRefreshingToken = false;
  private tokenRefreshSubscribers: Array<(success: boolean) => void> = [];

  private constructor() {
    // Use a separate base URL if auth service is hosted separately
    const AUTH_API_URL = 'http://localhost:8001';
    this.client = new ApiClient(AUTH_API_URL);
    // Attempt to load and set the token from localStorage when AuthService is initialized
    const token = this.getAccessToken();
    if (token) {
      this.client.setAuthHeader(`Bearer ${token}`);
    }
  }

  // Singleton pattern to ensure we only have one instance
  public static getInstance(): AuthService {
    if (!AuthService.instance) {
      AuthService.instance = new AuthService();
    }
    return AuthService.instance;
  }

  // Login with email and password
  async login(credentials: LoginRequest): Promise<ApiResponse<AuthResponse>> {
    // The backend expects 'username' and 'password' for OAuth2PasswordRequestForm
    // We'll send email as username.
    const formData = new URLSearchParams();
    formData.append('username', credentials.email);
    formData.append('password', credentials.password);

    return this.client.post<AuthResponse>('/auth/login', formData, {
      headers: {
        'Content-Type': 'application/x-www-form-urlencoded',
      },
    });
  }

  // Register new user
  async signup(userData: SignupRequest): Promise<ApiResponse<UserData>> {
    return this.client.post<UserData>('/auth/signup', userData);
  }

  // Request OTP for email login
  async requestOtp(email: string): Promise<ApiResponse<{ message: string }>> {
    return this.client.post<{ message: string }>('/auth/request-otp', { email });
  }

  // Verify OTP and login
  async verifyOtp(data: OtpLoginRequest): Promise<ApiResponse<AuthResponse>> {
    return this.client.post<AuthResponse>('/auth/verify-otp', data);
  }

  // Request password reset
  async forgotPassword(data: ForgotPasswordRequest): Promise<ApiResponse<{ message: string }>> {
    return this.client.post<{ message: string }>('/auth/password/forgot', data);
  }

  // Reset password with token
  async resetPassword(data: ResetPasswordRequest): Promise<ApiResponse<{ message: string }>> {
    return this.client.post<{ message: string }>('/auth/password/reset', data);
  }

  // Method to refresh the access token
  public async refreshToken(): Promise<boolean> {
    if (this.isRefreshingToken) {
      return new Promise((resolve) => {
        this.tokenRefreshSubscribers.push(resolve);
      });
    }
    this.isRefreshingToken = true;

    const currentRefreshToken = this.getRefreshToken();
    if (!currentRefreshToken) {
      this.isRefreshingToken = false;
      this.notifyTokenRefreshSubscribers(false);
      this.clearAuthData();
      return false;
    }

    try {
      const API_BASE_URL = 'http://localhost:8001';
      const refreshEndpoint = `${API_BASE_URL}/auth/refresh`;

      const response = await fetch(refreshEndpoint, {
        method: 'POST',
        headers: { 'Content-Type': 'application/json' },
        body: JSON.stringify({ refresh_token: currentRefreshToken } as RefreshTokenRequestData),
      });

      if (!response.ok) {
        const errorData = await response.json().catch(() => ({ detail: 'Refresh token failed, and error response was not valid JSON.' }));
        console.error("Refresh token failed:", response.status, errorData);
        this.clearAuthData();
        this.isRefreshingToken = false;
        this.notifyTokenRefreshSubscribers(false);
        return false;
      }

      const responseData: AuthResponse = await response.json();

      if (responseData && responseData.access_token) {
        this.storeAuthToken(responseData); // This already calls setAuthHeader
        this.isRefreshingToken = false;
        this.notifyTokenRefreshSubscribers(true);
        
        // Dispatch token refresh success event
        if (typeof window !== 'undefined') {
          window.dispatchEvent(new CustomEvent('auth:token_refreshed', { 
            detail: { timestamp: new Date().toISOString() }
          }));
        }
        
        return true;
      } else {
        console.error("Refresh token response was OK but not in expected AuthResponse format:", responseData);
        this.clearAuthData();
        this.isRefreshingToken = false;
        this.notifyTokenRefreshSubscribers(false);
        return false;
      }

    } catch (error) {
      console.error('Error during token refresh request execution:', error);
      this.clearAuthData();
      this.isRefreshingToken = false;
      this.notifyTokenRefreshSubscribers(false);
      return false;
    }
  }

  private notifyTokenRefreshSubscribers(success: boolean): void {
    this.tokenRefreshSubscribers.forEach(callback => callback(success));
    this.tokenRefreshSubscribers = [];
  }

  // Verify email address with token
  async verifyEmail(token: string): Promise<ApiResponse<{ message: string }>> {
    return this.client.get<{ message: string }>(`/auth/email/verify/${token}`);
  }

  // Request a new email verification token
  async requestEmailVerification(): Promise<ApiResponse<{ message: string }>> {
    return this.client.post<{ message: string }>('/auth/email/request-verification');
  }

  // Get current user data
  async getCurrentUser(): Promise<ApiResponse<UserData>> {
    return this.client.get<UserData>('/auth/users/me');
  }

  // Logout (invalidate token)
  // async logout(): Promise<ApiResponse<{ message: string }>> {
  //   return this.client.post<{ message: string }>('/auth/logout');
  // }

  // Helper method to handle login response
  storeAuthToken(response: AuthResponse): void {
    if (response.access_token && response.user && typeof response.access_token === 'string') {
      localStorage.setItem('access_token', response.access_token);
      localStorage.setItem('refresh_token', response.refresh_token);
      localStorage.setItem('user', JSON.stringify(response.user));
      localStorage.setItem('isAuthenticated', 'true');
      localStorage.setItem('isUserVerified', String(response.user.email_verified));
      // Consistently use user.mfa_enabled for mfa_required in localStorage,
      // as refresh token response might not have a top-level mfa_required field.
      localStorage.setItem('mfa_required', String(response.user.mfa_enabled));
      // Update ApiClient to use the new token
      this.client.setAuthHeader(`Bearer ${response.access_token}`);
    } else {
      console.error("AuthService.storeAuthToken: access_token or user is missing or invalid in the response.", response);
      this.clearAuthData(); 
    }
  }

  // Helper method to clear auth data on logout
  clearAuthData(): void {
    localStorage.removeItem('access_token'); // Corrected key to match getAccessToken
    localStorage.removeItem('refresh_token');
    localStorage.removeItem('user');
    localStorage.removeItem('isAuthenticated');
    localStorage.removeItem('isUserVerified'); // Fixed to match the actual key being used
    localStorage.removeItem('mfa_required');
    // Clear auth header in ApiClient
    this.client.clearAuthHeader();
  }

  // Check if user is logged in
  isLoggedIn(): boolean {
    return localStorage.getItem('access_token') !== null; // Corrected key to match getAccessToken
  }

  // Get access token
  getAccessToken(): string | null {
    return localStorage.getItem('access_token'); // Corrected key to 'access_token' from 'accessToken'
  }

  // Get refresh token
  getRefreshToken(): string | null {
    return localStorage.getItem('refresh_token');
  }

  // Get MFA status
  isMfaRequired(): boolean {
    return localStorage.getItem('mfa_required') === 'true';
  }

  // Get current user from local storage
  getUser(): UserData | null {
    const userJson = localStorage.getItem('user');
    if (!userJson) return null;
    
    try {
      return JSON.parse(userJson) as UserData;
    } catch (e) {
      console.error('Error parsing user data from localStorage', e);
      return null;
    }
  }
}

export default AuthService.getInstance();