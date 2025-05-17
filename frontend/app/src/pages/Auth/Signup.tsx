import React, { useState } from 'react';
import AuthLayout from '../../components/Auth/AuthLayout.tsx';
import GradientButton from '../../components/Auth/GradientButton.tsx';
import PasswordInput from '../../components/Auth/PasswordInput.tsx';
import { Link, useNavigate } from 'react-router-dom';
import AuthService from '../../services/api/AuthService.ts'; // Corrected import path
import { useAuthContext } from '../../contexts/AuthContext.tsx'; // Import the auth context hook

const Signup: React.FC = () => {
  const [username, setUsername] = useState('');
  const [email, setEmail] = useState('');
  const [password, setPassword] = useState('');
  const [confirmPassword, setConfirmPassword] = useState('');
  const [isSubmitting, setIsSubmitting] = useState(false);
  const [error, setError] = useState<string | null>(null);
  const [successMessage, setSuccessMessage] = useState<string | null>(null);
  const navigate = useNavigate();
  const { login } = useAuthContext(); // Get the login function from context

  const handleSubmit = async (e: React.FormEvent) => {
    e.preventDefault();
    
    // Reset messages
    setError(null);
    setSuccessMessage(null);

    // Form validations
    if (!username.trim()) {
      setError("Username is required");
      return;
    }
    
    if (password !== confirmPassword) {
      setError("Passwords don't match");
      return;
    }
    
    if (password.length < 8) {
      setError("Password should be at least 8 characters");
      return;
    }
    
    setIsSubmitting(true);
    
    try {
      // Configure the client to use your localhost endpoint
      const apiBaseUrl = 'http://localhost:8001';
      const client = AuthService['client'];
      client.setBaseUrl(apiBaseUrl);
      
      // Call signup endpoint with username, email, and password
      const response = await AuthService.signup({
        username,
        email,
        password,
      });
      
      if (!response.success || !response.data) {
        throw new Error(response.error || 'Failed to create account');
      }
      
      // Handle successful signup
      setSuccessMessage(
        "Account created successfully! Please check your email for verification instructions."
      );
      
      // Log the user in (they are authenticated but not verified yet)
      // Assuming the signup response includes tokens and user data similar to login
      // And a flag like `is_verified` which would be false initially.
      // login(response.data.access_token, response.data.refresh_token, response.data.user, false);
        
      // Wait a moment for the success message to be visible, then navigate
      setTimeout(() => {
        navigate('/login');
      }, 2000);

    } catch (error) {
      console.error('Signup error:', error);
      setError(
        error instanceof Error
          ? error.message
          : 'Failed to create account. Please try again.'
      );
    } finally {
      setIsSubmitting(false);
    }
  };

  return (
    <AuthLayout title="Create Account">
      <form onSubmit={handleSubmit} className="space-y-6 text-left">
        {error && (
          <div className="p-3 mb-4 text-sm text-theme-alert bg-theme-alert/20 border border-theme-alert rounded-lg" role="alert">
            <span className="block sm:inline">{error}</span>
          </div>
        )}
        
        {successMessage && (
          <div className="p-3 mb-4 text-sm text-theme-success bg-theme-success/20 border border-theme-success rounded-lg" role="alert">
            <span className="block sm:inline">{successMessage}</span>
          </div>
        )}
        
        <div>
          <label htmlFor="username" className="block text-sm font-medium text-theme-text-secondary text-left">Username</label>
          <input
            id="username"
            name="username"
            type="text"
            autoComplete="username"
            required
            value={username}
            onChange={(e) => setUsername(e.target.value)}
            className="mt-1 block w-full px-3 py-2 border border-theme-text-secondary/50 rounded-md shadow-sm placeholder-theme-text-secondary focus:outline-none focus:ring-theme-success focus:border-theme-success sm:text-sm bg-theme-glassy text-theme-text-primary"
          />
        </div>
        
        <div>
          <label htmlFor="email" className="block text-sm font-medium text-theme-text-secondary text-left">Email address</label>
          <input
            id="email"
            name="email"
            type="email"
            autoComplete="email"
            required
            value={email}
            onChange={(e) => setEmail(e.target.value)}
            className="mt-1 block w-full px-3 py-2 border border-theme-text-secondary/50 rounded-md shadow-sm placeholder-theme-text-secondary focus:outline-none focus:ring-theme-success focus:border-theme-success sm:text-sm bg-theme-glassy text-theme-text-primary"
          />
        </div>
        
        <PasswordInput
          id="password"
          name="password"
          label="Password"
          value={password}
          onChange={(e) => setPassword(e.target.value)}
          autoComplete="new-password"
          disabled={isSubmitting}
        />
        
        <PasswordInput
          id="confirm-password"
          name="confirm-password"
          label="Confirm Password"
          value={confirmPassword}
          onChange={(e) => setConfirmPassword(e.target.value)}
          autoComplete="new-password"
          disabled={isSubmitting}
        />
        
        <GradientButton
          type="submit"
          isLoading={isSubmitting}
          loadingText="Signing up..."
          disabled={isSubmitting}
        >
          Sign Up
        </GradientButton>
        
        <p className="text-sm text-center text-theme-text-secondary">
          Already have an account? <Link to="/login" className="font-medium text-theme-success hover:text-accent-light">Log In</Link>
        </p>
      </form>
    </AuthLayout>
  );
};

export default Signup;