import React, { useState, useEffect } from 'react';
import { useNavigate, useSearchParams, Link } from 'react-router-dom';
import AuthLayout from '../../components/Auth/AuthLayout.tsx';
import GradientButton from '../../components/Auth/GradientButton.tsx';
import PasswordInput from '../../components/Auth/PasswordInput.tsx';
import type { ResetPasswordRequest } from '../../services/api/AuthService'; // Type import remains named
import AuthService from '../../services/api/AuthService.ts'; // AuthService is a default import

const ResetPassword: React.FC = () => {
  const [password, setPassword] = useState('');
  const [confirmPassword, setConfirmPassword] = useState('');
  const [isSubmitting, setIsSubmitting] = useState(false);
  const [error, setError] = useState<string | null>(null);
  const [message, setMessage] = useState<string | null>(null);
  
  const navigate = useNavigate();
  const [searchParams] = useSearchParams();
  const [token, setToken] = useState<string | null>(null);

  useEffect(() => {
    const tokenFromQuery = searchParams.get('token');
    // console.log('Token from query:', tokenFromQuery);
    if (tokenFromQuery) {
      setToken(tokenFromQuery);
    } else {
      setError("No reset token found. Please use the link from your email.");
    }
  }, [searchParams]);

  const handleSubmit = async (e: React.FormEvent) => {
    e.preventDefault();
    
    if (!token) {
      setError("Password reset token is missing. Please use the link provided in your email.");
      return;
    }

    if (password !== confirmPassword) {
      setError("Passwords don't match");
      return;
    }
    
    if (password.length < 8) {
      setError("Password should be at least 8 characters long.");
      return;
    }
    
    setError(null);
    setMessage(null);
    setIsSubmitting(true);
    
    try {
      const payload: ResetPasswordRequest = { token, password };
      const response = await AuthService.resetPassword(payload);
      
      if (response.success) {
        setMessage(response.data?.message || 'Password has been reset successfully!');
        setTimeout(() => navigate('/login'), 3000);
      } else {
        setError(response.error || 'Failed to reset password. The link may be invalid or expired.');
      }
    } catch (err: any) {
      console.error('Error resetting password:', err);
      if (err.response && err.response.data && err.response.data.detail) {
        setError(err.response.data.detail);
      } else if (err instanceof Error) {
        setError(err.message);
      } else {
        setError('An unexpected error occurred. Please try again or request a new reset link.');
      }
    } finally {
      setIsSubmitting(false);
    }
  };

  return (
    <AuthLayout title="Set New Password">
      <form onSubmit={handleSubmit} className="space-y-6 text-left">
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
        
        <PasswordInput
          id="password"
          name="password"
          label="New Password"
          value={password}
          onChange={(e) => setPassword(e.target.value)}
          autoComplete="new-password"
          disabled={isSubmitting || !!message}
        />
        
        <PasswordInput
          id="confirmPassword"
          name="confirmPassword"
          label="Confirm New Password"
          value={confirmPassword}
          onChange={(e) => setConfirmPassword(e.target.value)}
          autoComplete="new-password"
          disabled={isSubmitting || !!message}
        />
        
        <GradientButton
          type="submit"
          isLoading={isSubmitting}
          loadingText="Resetting..."
          disabled={isSubmitting || !token || !!message}
        >
          Reset Password
        </GradientButton>
        
        <p className="text-sm text-center text-theme-text-secondary">
          {message ? (
            <Link to="/login" className="font-medium text-theme-success hover:text-accent-light">Proceed to Login</Link>
          ) : (
            <>Remember your password? <Link to="/login" className="font-medium text-theme-success hover:text-accent-light">Back to Login</Link></>
          )}
        </p>
      </form>
    </AuthLayout>
  );
};

export default ResetPassword;