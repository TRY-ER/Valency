import React, { useState } from 'react';
import AuthLayout from '../../components/Auth/AuthLayout.tsx';
import GradientButton from '../../components/Auth/GradientButton.tsx';
import PasswordInput from '../../components/Auth/PasswordInput.tsx';
import { Link, useNavigate } from 'react-router-dom';
import AuthService from '../../services/api/AuthService.ts'; // Corrected import path
import { useAuthContext } from '../../contexts/AuthContext.tsx'; // Import the auth context hook

const Login: React.FC = () => {
  const [email, setEmail] = useState('');
  const [password, setPassword] = useState('');
  const [isSubmitting, setIsSubmitting] = useState(false);
  const [error, setError] = useState<string | null>(null); // For displaying login errors
  const navigate = useNavigate(); // For navigation
  const { login } = useAuthContext(); // Get the login function from context

  const handleSubmit = async (e: React.FormEvent) => {
    e.preventDefault();
    setIsSubmitting(true);
    setError(null); // Clear previous errors

    try {
      const response = await AuthService.login({ email, password });
      
      if (response.success && response.data) {
        // Use the login function from AuthContext
        login(response.data.access_token, response.data.refresh_token, response.data.user, response.data.user.email_verified);
        
        // console.log('Login successful:', response.data);
        if (!response.data.user.email_verified) {
          navigate('/verify-email'); // Redirect to verify email page if not verified
        } else if (response.data.mfa_required) {
          navigate('/mfa-verification'); // Redirect to MFA page if MFA is required
        } else {
          navigate('/dashboard'); // Redirect to dashboard or home page
        }
      } else {
        setError(response.error || 'Login failed. Please try again.');
      }
    } catch (err) {
      // Handle unexpected errors (e.g., network issues)
      setError(err instanceof Error ? err.message : 'An unexpected error occurred.');
    } finally {
      setIsSubmitting(false);
    }
  };

  return (
    <AuthLayout title="Login to Your Account">
      <form onSubmit={handleSubmit} className="space-y-6 text-left">
        {error && (
          <div className="p-3 mb-4 text-sm text-theme-alert bg-theme-alert/20 border border-theme-alert rounded-lg" role="alert">
            {error}
          </div>
        )}
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
            disabled={isSubmitting}
            className="mt-1 block w-full px-3 py-2 border border-theme-text-secondary/50 rounded-md shadow-sm placeholder-theme-text-secondary focus:outline-none focus:ring-theme-success focus:border-theme-success sm:text-sm bg-theme-glassy text-theme-text-primary"
          />
        </div>
        
        <PasswordInput
          id="password"
          name="password"
          label="Password"
          value={password}
          onChange={(e) => setPassword(e.target.value)}
          autoComplete="current-password"
          disabled={isSubmitting}
        />
        
        <div className="flex items-center justify-between">
          <div className="text-sm">
            <Link to="/forgot-password" className="font-medium text-theme-success hover:text-accent-light">Forgot your password?</Link>
          </div>
        </div>
        
        <GradientButton 
          type="submit" 
          isLoading={isSubmitting} 
          loadingText="Logging in..."
          disabled={isSubmitting}
        >
          Log In
        </GradientButton>
        
        <p className="text-sm text-center text-theme-text-secondary">
          Don't have an account? <Link to="/signup" className="font-medium text-theme-success hover:text-accent-light">Sign Up</Link>
        </p>
      </form>
    </AuthLayout>
  );
};

export default Login;