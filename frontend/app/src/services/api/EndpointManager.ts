import AuthService from './AuthService';
import ApiClient from './ApiClient';

/**
 * EndpointManager provides a centralized access point to different microservices
 * Each service can have its own base URL and specific endpoints
 */
class EndpointManager {
  private static instance: EndpointManager;
  
  // Store service clients by their name
  private services: Record<string, ApiClient> = {};
  
  // Auth service is handled separately since it has a more complex interface
  private authService = AuthService;
  
  private constructor() {
    // Initialize any default services here
    this.initializeDefaultServices();
  }

  // Singleton pattern
  public static getInstance(): EndpointManager {
    if (!EndpointManager.instance) {
      EndpointManager.instance = new EndpointManager();
    }
    return EndpointManager.instance;
  }

  private initializeDefaultServices(): void {
    // Initialize default service clients
    // You can customize these based on your environment configuration
    
    // Main API (for general purposes)
    const apiBaseUrl = import.meta.env.VITE_API_BASE_URL || 'http://localhost:3000/api';
    this.registerService('main', apiBaseUrl);
    
    // You can add more services here as needed
    // Example:
    // this.registerService('users', import.meta.env.VITE_USERS_API_URL);
    // this.registerService('products', import.meta.env.VITE_PRODUCTS_API_URL);
    // this.registerService('payments', import.meta.env.VITE_PAYMENTS_API_URL);
  }

  /**
   * Register a new service with its base URL
   */
  public registerService(serviceName: string, baseUrl: string): void {
    this.services[serviceName] = new ApiClient(baseUrl);
  }

  /**
   * Get an API client for a specific service
   */
  public getService(serviceName: string): ApiClient {
    if (!this.services[serviceName]) {
      throw new Error(`Service "${serviceName}" is not registered`);
    }
    return this.services[serviceName];
  }

  /**
   * Get the auth service
   */
  public getAuthService() {
    return this.authService;
  }
}

// Export a singleton instance
export default EndpointManager.getInstance();