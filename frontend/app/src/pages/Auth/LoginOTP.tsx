import React, { useState } from 'react';
import AuthLayout from '../../components/Auth/AuthLayout.tsx';
import { motion } from 'framer-motion';
import GradientButton from '../../components/Auth/GradientButton.tsx';
import { Link, useLocation, useNavigate } from 'react-router-dom';
import AuthService from '../../services/api/AuthService.ts'; // Corrected import path

const LoginOTP: React.FC = () => {
  const [otp, setOtp] = useState('');
  const [isVerifying, setIsVerifying] = useState(false);
  const [error, setError] = useState<string | null>(null);
  const navigate = useNavigate();
  const location = useLocation();
  const email = location.state?.email;

  const handleSubmit = async (e: React.FormEvent) => {
    e.preventDefault();
    
    if (!otp || otp.length < 6) {
      setError('Please enter a valid OTP code');
      return;
    }

    if (!email) {
      setError('Email not found. Please restart the OTP login process.');
      return;
    }
    
    setIsVerifying(true);
    setError(null);
    
    try {
      // Call our AuthService to verify the OTP
      const response = await AuthService.verifyOtp({ email, otp });
      
      if (!response.success) {
        throw new Error(response.error || 'Failed to verify OTP');
      }
      
      // Store authentication data
      if (response.data) {
        AuthService.storeAuthToken(response.data);
        
        // Redirect to home page
        navigate('/');
      } else {
        throw new Error('Invalid response from server');
      }
    } catch (error) {
      console.error('OTP verification failed:', error);
      setError(error instanceof Error ? error.message : 'Invalid OTP. Please try again.');
    } finally {
      setIsVerifying(false);
    }
  };

  const handleResendOTP = async () => {
    if (!email) {
      setError('Email not found. Please go back to the previous step.');
      return;
    }
    
    try {
      // Call our AuthService to request a new OTP
      const response = await AuthService.requestOtp(email);
      
      if (!response.success) {
        throw new Error(response.error || 'Failed to resend OTP');
      }
      
      // Show success message
      alert(response.message || `OTP resent to ${email}`);
    } catch (error) {
      console.error('Failed to resend OTP:', error);
      setError(error instanceof Error ? error.message : 'Failed to resend OTP. Please try again.');
    }
  };

  return (
    <AuthLayout title="Enter OTP">
      <form onSubmit={handleSubmit} className="space-y-6 text-left">
        {email ? (
          <p className="text-sm text-center text-theme-text-secondary">
            Enter the 6-digit code sent to <span className="font-medium text-theme-text-primary">{email}</span>
          </p>
        ) : (
          <p className="text-sm text-center text-theme-alert">
            No email provided. Please <Link to="/request-otp-email" className="underline hover:text-accent-light">restart the process</Link>.
          </p>
        )}
        
        {error && (
          <div className="p-3 mb-4 text-sm text-theme-alert bg-theme-alert/20 border border-theme-alert rounded-lg" role="alert">
            <span className="block sm:inline">{error}</span>
          </div>
        )}
        
        <div>
          <label htmlFor="otp" className="block text-sm font-medium text-theme-text-secondary text-left">One-Time Password</label>
          <input
            id="otp"
            name="otp"
            type="text"
            inputMode="numeric"
            autoComplete="one-time-code"
            required
            value={otp}
            onChange={(e) => {
              const value = e.target.value.replace(/[^0-9]/g, '');
              setOtp(value);
            }}
            className="mt-1 block w-full px-3 py-2 border border-theme-text-secondary/50 rounded-md shadow-sm placeholder-theme-text-secondary focus:outline-none focus:ring-theme-success focus:border-theme-success sm:text-sm bg-theme-glassy text-theme-text-primary text-center tracking-widest"
            maxLength={6}
            placeholder="••••••"
            disabled={isVerifying}
          />
        </div>
        
        <GradientButton 
          type="submit"
          isLoading={isVerifying}
          loadingText="Verifying..."
          disabled={isVerifying}
        >
          Verify OTP
        </GradientButton>
        
        <div className="flex flex-col sm:flex-row justify-between items-center gap-2">
          <motion.button 
            type="button"
            onClick={handleResendOTP}
            whileTap={{ scale: 0.95 }}
            className="text-sm font-medium py-2 px-4 rounded-md border border-theme-success/50 bg-theme-glassy text-theme-success hover:bg-theme-success/10 focus:outline-none focus:ring-2 focus:ring-offset-2 focus:ring-theme-success transition-colors"
            disabled={isVerifying}
          >
            Resend OTP
          </motion.button>
          
          <Link to="/login" className="text-sm font-medium text-theme-text-secondary hover:text-theme-success">
            Login with Password
          </Link>
        </div>
      </form>
    </AuthLayout>
  );
};

export default LoginOTP;