import React, { useState } from 'react';
import { Link } from 'react-router-dom';
import AuthLayout from '../../components/Auth/AuthLayout.tsx';
import GradientButton from '../../components/Auth/GradientButton.tsx';
import AuthService from '../../services/api/AuthService.ts'; // Corrected import path

const ForgotPassword: React.FC = () => {
  const [email, setEmail] = useState('');
  const [isSubmitting, setIsSubmitting] = useState(false);
  const [message, setMessage] = useState<string | null>(null); // For success/error messages
  const [error, setError] = useState<string | null>(null); // For error messages

  const handleSubmit = async (e: React.FormEvent) => {
    e.preventDefault();
    setIsSubmitting(true);
    setMessage(null);
    setError(null);
    
    try {
      const response = await AuthService.forgotPassword({ email });
      
      if (response.success) {
        setMessage(response.data?.message || "If an account with this email exists, a password reset link has been sent.");
      } else {
        setError(response.error || 'Failed to send reset instructions. Please try again.');
      }
    } catch (err: any) {
      console.error('Error requesting password reset:', err);
      // Check if the error has a response and data with a detail message
      if (err.response && err.response.data && err.response.data.detail) {
        setError(err.response.data.detail);
      } else if (err instanceof Error) {
        setError(err.message);
      } else {
        setError('An unexpected error occurred. Please try again.');
      }
    } finally {
      setIsSubmitting(false);
    }
  };

  return (
    <AuthLayout title="Reset Your Password">
      <form onSubmit={handleSubmit} className="space-y-6 text-left">
        <p className="text-sm text-theme-text-secondary">
          Enter your email address and we'll send you instructions to reset your password.
        </p>
        
        {message && (
          <div className="p-3 mb-4 text-sm text-theme-success bg-theme-success/20 border border-theme-success rounded-lg" role="alert">
            <span className="block sm:inline">{message}</span>
          </div>
        )}

        {error && (
          <div className="p-3 mb-4 text-sm text-theme-alert bg-theme-alert/20 border border-theme-alert rounded-lg" role="alert">
            <span className="block sm:inline">{error}</span>
          </div>
        )}
        
        <div>
          <label htmlFor="email" className="block text-sm font-medium text-theme-text-secondary text-left">
            Email address
          </label>
          <input
            id="email"
            name="email"
            type="email"
            autoComplete="email"
            required
            value={email}
            onChange={(e) => setEmail(e.target.value)}
            className="mt-1 block w-full px-3 py-2 border border-theme-text-secondary/50 rounded-md shadow-sm placeholder-theme-text-secondary focus:outline-none focus:ring-theme-success focus:border-theme-success sm:text-sm bg-theme-glassy text-theme-text-primary"
            placeholder="Enter your email"
            disabled={isSubmitting}
          />
        </div>
        
        <GradientButton
          type="submit"
          isLoading={isSubmitting}
          loadingText="Sending..."
        >
          Send Reset Instructions
        </GradientButton>
        
        <p className="text-sm text-center text-theme-text-secondary">
          Remember your password? <Link to="/login" className="font-medium text-theme-success hover:text-accent-light">Back to Login</Link>
        </p>
      </form>
    </AuthLayout>
  );
};

export default ForgotPassword;