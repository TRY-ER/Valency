import React, { useEffect, useState } from 'react';
import AuthLayout from '../../components/Auth/AuthLayout.tsx';
import { motion } from 'framer-motion';
import { Link, useSearchParams, useNavigate } from 'react-router-dom';
import AuthService from '../../services/api/AuthService.ts'; // Corrected import path
import { useAuthContext } from '../../contexts/AuthContext.tsx'; // Import the auth context hook

const VerifyEmail: React.FC = () => {
  const [searchParams] = useSearchParams();
  const navigate = useNavigate();
  const [status, setStatus] = useState<'loading' | 'success' | 'error'>('loading');
  const [message, setMessage] = useState<string>('Verifying your email...');
  const { setIsVerifiedState, isAuthenticated } = useAuthContext(); // Get setIsVerifiedState and isAuthenticated

  useEffect(() => {
    const token = searchParams.get('token');

    if (!isAuthenticated) {
      // If user is not even authenticated (e.g. landed here directly without signup/login flow)
      // or if they logged out and tried to use an old link.
      setStatus('error');
      setMessage('You need to be logged in to verify your email. Please login and try again or use the link from your email after logging in.');
      // Optionally, redirect to login after a delay or provide a login link
      // setTimeout(() => navigate('/login'), 5000);
      return;
    }

    if (!token) {
      setStatus('error');
      setMessage('No verification token found. Please check your email link or request a new one.');
      return;
    }

    const verifyToken = async () => {
      try {
        const response = await AuthService.verifyEmail(token);
        if (response.success) {
          setStatus('success');
          setMessage(response.data?.message || 'Email verified successfully!');
          setIsVerifiedState(true); // Update the verification status in the context
          // Navigate to dashboard or next step after a short delay
          setTimeout(() => navigate('/dashboard'), 2000); 
        } else {
          setStatus('error');
          // Use the error message from the response if available
          setMessage(response.error || 'Failed to verify email. The link may be invalid or expired.');
        }
      } catch (error: any) { // Added :any to access error.response
        setStatus('error');
        // Check if the error has a response and data with a detail message (for 400 errors)
        if (error.response && error.response.data && error.response.data.detail) {
          setMessage(error.response.data.detail);
        } else if (error instanceof Error) {
          setMessage(error.message);
        } else {
          setMessage('An unexpected error occurred during email verification.');
        }
      }
    };

    verifyToken();
  }, [searchParams, navigate, setIsVerifiedState, isAuthenticated]);

  // Handle resend verification email
  const handleResendEmail = async () => {
    setStatus('loading');
    setMessage('Sending a new verification email...');
    try {
      const response = await AuthService.requestEmailVerification();
      // Check if the response is successful
      // console.log('Resend email response:', response);
      if (response.success) {
        setStatus('success');
        setMessage(response.data?.message || 'A new verification email has been sent. Please check your inbox.');
      } else {
        setStatus('error');
        setMessage(response.error || 'Failed to send verification email. Please try again.');
      }
    } catch (error: any) {
      setStatus('error');
      if (error.response && error.response.data && error.response.data.detail) {
        setMessage(error.response.data.detail);
      } else if (error instanceof Error) {
        setMessage(error.message);
      } else {
        setMessage('An unexpected error occurred while resending the verification email.');
      }
    }
  };

  return (
    <AuthLayout title="Verify Your Email">
      <div className="text-center">
        {status === 'loading' && (
          <div className="flex justify-center mb-4">
            <div className="animate-spin rounded-full h-10 w-10 border-b-2 border-theme-success"></div> {/* Use theme-success for spinner */}
          </div>
        )}

        <p className={`text-theme-text-secondary mb-6 ${status === 'error' ? 'text-theme-alert' : status === 'success' ? 'text-theme-success' : ''}`}>
          {message}
        </p>

        {status === 'success' && (
          <motion.button
            whileHover={{ scale: 1.05 }}
            whileTap={{ scale: 0.95 }}
            onClick={() => navigate('/login')}
            className="w-full flex justify-center py-2 px-4 border border-transparent rounded-md shadow-sm text-sm font-medium text-white bg-major hover:bg-opacity-80 focus:outline-none focus:ring-2 focus:ring-offset-2 focus:ring-major mb-4"
          >
            Proceed to Login
          </motion.button>
        )}

        {status === 'error' && (
          <motion.button
            whileHover={{ scale: 1.05 }}
            whileTap={{ scale: 0.95 }}
            onClick={handleResendEmail}
            className="w-full flex justify-center py-2 px-4 border border-transparent rounded-md shadow-sm text-sm font-medium text-white bg-major hover:bg-opacity-80 focus:outline-none focus:ring-2 focus:ring-offset-2 focus:ring-major mb-4"
          >
            Resend Verification Email
          </motion.button>
        )}

        <p className="text-sm text-theme-text-secondary">
          Already verified? <Link to="/login" className="font-medium text-theme-success hover:text-accent-light">Log In</Link>
        </p>
      </div>
    </AuthLayout>
  );
};

export default VerifyEmail;