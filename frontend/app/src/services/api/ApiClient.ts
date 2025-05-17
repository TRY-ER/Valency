/**
 * Base API client for making HTTP requests
 * This handles common configuration, error handling, and request/response processing
 */

const API_BASE_URL = 'http://localhost:8001';

export interface ApiResponse<T = any> {
  data: T;
  success: boolean;
  message?: string;
  error?: string;
}

class ApiClient {
  private baseUrl: string;
  private authToken: string | null = null; // Store token internally

  constructor(baseUrl = API_BASE_URL) {
    this.baseUrl = baseUrl;
    try {
      const tokenFromStorage = localStorage.getItem('access_token'); // Corrected key
      if (tokenFromStorage) {
        this.authToken = `Bearer ${tokenFromStorage}`; // Set in Bearer format
      } else {
        this.authToken = null;
      }
    } catch (e) {
      console.error("[ApiClient] Constructor: Error accessing localStorage or processing token:", e);
      this.authToken = null; // Ensure authToken is null in case of error
    }
  }

  // Method to set the authorization header dynamically
  public setAuthHeader(token: string | null): void {
    this.authToken = token;
  }

  // Method to clear the authorization header
  public clearAuthHeader(): void {
    this.authToken = null;
  }

  // Helper method to get authentication headers
  private getHeaders(customHeaders?: HeadersInit): HeadersInit {
    const baseHeaders: Record<string, string> = {
      'Content-Type': 'application/json',
    };

    // Add auth token if available
    if (this.authToken) {
      baseHeaders['Authorization'] = `Bearer ${this.authToken}`; // Use the internally stored token
    }

    // Merge with custom headers. Custom headers will override baseHeaders if keys conflict.
    // This is important if customHeaders also tries to set 'Content-Type' or 'Authorization'
    return {
      ...baseHeaders,
      ...(customHeaders as Record<string, string>), // Cast customHeaders for proper merging
    };
  }

  // Generic request method with error handling
  private async request<T>(
    endpoint: string, 
    method: string = 'GET', 
    body?: any,
    customHeaders?: HeadersInit // Add customHeaders parameter
  ): Promise<ApiResponse<T>> {
    const url = `${this.baseUrl}${endpoint}`;
    
    try {
      const response = await fetch(url, {
        method,
        headers: this.getHeaders(customHeaders), // Pass customHeaders to getHeaders
        body: body instanceof URLSearchParams ? body : (body ? JSON.stringify(body) : undefined),
        mode: 'cors', // Explicitly set CORS mode
        credentials: 'include', // Include cookies and authentication headers
      });

      const data = await response.json();
      
      if (!response.ok) {
        // Handle FastAPI error format which typically returns a "detail" field
        const errorMessage = data.detail || data.message || `Error ${response.status}: ${response.statusText}`;
        return {
          success: false,
          data: null as any,
          error: errorMessage,
        };
      }

      // Handle FastAPI success response which may not have a message field
      return {
        success: true,
        data,
        message: data.message || 'Operation completed successfully',
      };
    } catch (error) {
      console.error('API request failed:', error);
      return {
        success: false,
        data: null as any,
        error: error instanceof Error ? error.message : 'Unknown error occurred',
      };
    }
  }
  
  // HTTP method wrappers
  async get<T>(endpoint: string): Promise<ApiResponse<T>> {
    return this.request<T>(endpoint, 'GET');
  }
  
  async post<T>(endpoint: string, data?: any, options?: { headers?: HeadersInit }): Promise<ApiResponse<T>> {
    return this.request<T>(endpoint, 'POST', data, options?.headers);
  }
  
  async put<T>(endpoint: string, data?: any, options?: { headers?: HeadersInit }): Promise<ApiResponse<T>> {
    return this.request<T>(endpoint, 'PUT', data, options?.headers);
  }
  
  async delete<T>(endpoint: string, options?: { headers?: HeadersInit }): Promise<ApiResponse<T>> {
    return this.request<T>(endpoint, 'DELETE', undefined, options?.headers);
  }

  // Set a new base URL (useful for switching environments or services)
  setBaseUrl(url: string): void {
    this.baseUrl = url;
  }
}

export default ApiClient;