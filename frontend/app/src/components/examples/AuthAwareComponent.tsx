import React, { useState } from 'react';
import { useAuth, useAuthEvents } from '../../hooks/useAuth.ts';
import { uniprotSearchUniprotkb } from '../../services/api/mcpToolsService.js';

/**
 * Example component demonstrating how to use the enhanced auth system
 * This component will automatically handle token refresh when making API calls
 */
const AuthAwareComponent: React.FC = () => {
  const { isAuthenticated, user, logout } = useAuth();
  const [searchResult, setSearchResult] = useState<any>(null);
  const [loading, setLoading] = useState(false);
  const [error, setError] = useState<string | null>(null);

  // Handle auth events
  useAuthEvents({
    onLogout: (detail) => {
      console.log('Component received logout event:', detail);
      setSearchResult(null);
      setError(`Logged out due to: ${detail.reason}`);
    },
    onTokenRefresh: () => {
      console.log('Token refreshed, clearing any auth-related errors');
      if (error?.includes('authentication')) {
        setError(null);
      }
    }
  });

  const handleSearchUniProt = async () => {
    if (!isAuthenticated) {
      setError('Please log in to search UniProt');
      return;
    }

    setLoading(true);
    setError(null);

    try {
      // This call will automatically handle token refresh if needed
      const result = await uniprotSearchUniprotkb({
        query_string: 'insulin',
        size: 5
      });
      
      setSearchResult(result);
    } catch (err: any) {
      setError(err.message || 'Search failed');
    } finally {
      setLoading(false);
    }
  };

  if (!isAuthenticated) {
    return (
      <div className="p-4 border rounded-lg">
        <p className="text-gray-600">Please log in to use this component</p>
      </div>
    );
  }

  return (
    <div className="p-4 border rounded-lg space-y-4">
      <div className="flex justify-between items-center">
        <h3 className="text-lg font-semibold">Auth-Aware UniProt Search</h3>
        <div className="flex items-center space-x-2">
          <span className="text-sm text-gray-600">
            Welcome, {user?.username || user?.email}
          </span>
          <button
            onClick={logout}
            className="text-sm text-red-600 hover:text-red-800"
          >
            Logout
          </button>
        </div>
      </div>

      <button
        onClick={handleSearchUniProt}
        disabled={loading}
        className="px-4 py-2 bg-blue-500 text-white rounded hover:bg-blue-600 disabled:opacity-50"
      >
        {loading ? 'Searching...' : 'Search UniProt for Insulin'}
      </button>

      {error && (
        <div className="p-3 bg-red-100 border border-red-400 text-red-700 rounded">
          {error}
        </div>
      )}

      {searchResult && (
        <div className="p-3 bg-green-100 border border-green-400 text-green-700 rounded">
          <p className="font-semibold">Search successful!</p>
          <p className="text-sm">Found {searchResult.results?.length || 0} results</p>
        </div>
      )}
    </div>
  );
};

export default AuthAwareComponent;
