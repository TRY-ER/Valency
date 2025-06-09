# Authentication & Token Refresh Implementation

This document describes the comprehensive token refresh system implemented across all API services in the application.

## Overview

The authentication system now automatically handles expired access tokens by attempting to refresh them using the refresh token. If the refresh token is also expired, the system will automatically log out the user and clear all authentication data.

## Components

### 1. AuthService (Core Authentication)
- **Location**: `src/services/api/AuthService.ts`
- **Features**:
  - Singleton pattern for consistent state management
  - Automatic token refresh with subscriber pattern to prevent multiple simultaneous refresh attempts
  - Token storage and retrieval from localStorage
  - Event dispatching for successful token refresh

### 2. API Services with Token Refresh

#### MCP Tools Service
- **Location**: `src/services/api/mcpToolsService.js`
- **Features**:
  - Axios response interceptor for 401 error handling
  - Automatic token refresh and request retry
  - Custom logout event dispatch on refresh failure

#### Agent Service
- **Location**: `src/services/api/agentService.js`
- **Features**:
  - Similar axios-based token refresh implementation
  - Handles streaming and regular API calls

#### ApiClient (Fetch-based)
- **Location**: `src/services/api/ApiClient.ts`
- **Features**:
  - Fetch-based token refresh handling
  - Dynamic AuthService import to avoid circular dependencies
  - Excludes auth endpoints from refresh attempts

### 3. AuthContext
- **Location**: `src/contexts/AuthContext.tsx`
- **Features**:
  - Listens for storage changes (multi-tab synchronization)
  - Listens for custom logout events from services
  - Automatic state synchronization

### 4. Utility Functions
- **Location**: `src/utils/authUtils.ts`
- **Features**:
  - Global logout trigger function
  - Token validation utilities
  - API call wrapper with automatic token refresh
  - Auth event listener creator

### 5. React Hooks
- **Location**: `src/hooks/useAuth.ts`
- **Features**:
  - Enhanced auth hook with event handling
  - Specialized hook for auth events only

## How It Works

### Token Refresh Flow

1. **API Request**: Any API service makes a request with the current access token
2. **401 Response**: If the server returns 401 (Unauthorized), the interceptor catches it
3. **Refresh Attempt**: The service calls `AuthService.refreshToken()`
4. **Refresh Success**: If successful, the original request is retried with the new token
5. **Refresh Failure**: If refresh fails, auth data is cleared and logout event is dispatched

### Event System

The system uses custom events for cross-component communication:

- `auth:logout` - Dispatched when logout occurs (manual or automatic)
- `auth:token_refreshed` - Dispatched when token refresh succeeds

### Multi-Tab Synchronization

The AuthContext listens for `storage` events to synchronize auth state across browser tabs.

## Usage Examples

### Using the Enhanced Auth Hook

```tsx
import { useAuth } from '../hooks/useAuth.ts';

const MyComponent = () => {
  const { isAuthenticated, user, logout } = useAuth();
  
  // Component automatically handles logout events
  return (
    <div>
      {isAuthenticated ? `Welcome ${user?.username}` : 'Please log in'}
    </div>
  );
};
```

### Handling Auth Events

```tsx
import { useAuthEvents } from '../hooks/useAuth.ts';

const MyComponent = () => {
  useAuthEvents({
    onLogout: (detail) => {
      console.log('User logged out:', detail.reason);
      // Handle logout (e.g., show notification, redirect)
    },
    onTokenRefresh: () => {
      console.log('Token refreshed successfully');
      // Handle successful refresh (e.g., update UI state)
    }
  });
  
  return <div>My Component</div>;
};
```

### Manual Token Operations

```typescript
import { 
  triggerGlobalLogout, 
  validateCurrentToken, 
  withTokenRefresh 
} from '../utils/authUtils.ts';

// Manual logout
triggerGlobalLogout('user_action', 'logout_button');

// Validate current token
const isValid = await validateCurrentToken();

// Wrap API calls with automatic refresh
const wrappedApiCall = withTokenRefresh(someApiFunction);
const result = await wrappedApiCall(params);
```

## Configuration

### Service Base URLs

- Auth Service: `http://localhost:8001`
- MCP Tools Service: `http://localhost:8000`
- Agent Service: `http://localhost:8015`

### Token Storage Keys

- `access_token` - Current access token
- `refresh_token` - Refresh token for getting new access tokens
- `user` - Serialized user data
- `isAuthenticated` - Boolean auth status
- `isUserVerified` - Boolean email verification status
- `mfa_required` - Boolean MFA requirement status

## Error Handling

### Automatic Scenarios

1. **401 on API Call**: Automatic refresh attempt
2. **Refresh Token Expired**: Automatic logout and data clearing
3. **Network Errors**: Proper error propagation with context

### Manual Scenarios

1. **User Logout**: Manual trigger via AuthContext
2. **Security Logout**: Manual trigger via utility functions

## Security Considerations

1. **Single Refresh Attempt**: Each request only attempts refresh once to prevent loops
2. **Concurrent Refresh Protection**: Subscriber pattern prevents multiple simultaneous refresh attempts
3. **Secure Storage**: Tokens stored in localStorage with proper cleanup
4. **Cross-Tab Sync**: Logout in one tab affects all tabs immediately

## Troubleshooting

### Common Issues

1. **Infinite Refresh Loops**: Prevented by `_retry` flag and `isRetry` parameter
2. **Circular Dependencies**: Resolved with dynamic imports in ApiClient
3. **Race Conditions**: Managed with subscriber pattern in AuthService
4. **Storage Sync**: Handled with storage event listeners

### Debug Logging

All services include comprehensive console logging for debugging:
- Token refresh attempts
- Success/failure status
- Error details
- Event dispatching

## Future Enhancements

1. **Token Expiry Prediction**: Proactive refresh before expiration
2. **Retry Strategies**: Exponential backoff for network failures
3. **Security Headers**: Additional security headers for API requests
4. **Token Rotation**: Enhanced security with token rotation policies
