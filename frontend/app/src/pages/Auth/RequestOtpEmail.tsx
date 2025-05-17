import React, { useState } from 'react';
import { useNavigate, Link } from 'react-router-dom';
import AuthLayout from '../../components/Auth/AuthLayout.tsx';
import GradientButton from '../../components/Auth/GradientButton.tsx';

const RequestOtpEmail: React.FC = () => {
  const [email, setEmail] = useState('');
  const [isSubmitting, setIsSubmitting] = useState(false);
  const navigate = useNavigate();

  const handleSubmit = async (e: React.FormEvent) => {
    e.preventDefault();
    setIsSubmitting(true);
    
    try {
      // In a real app, you would call your API endpoint to send OTP
      // const response = await fetch('/api/auth/send-otp', {
      //   method: 'POST',
      //   headers: { 'Content-Type': 'application/json' },
      //   body: JSON.stringify({ email }),
      // });
      // 
      // if (!response.ok) throw new Error('Failed to send OTP');
      
      // Simulating API call delay
      await new Promise(resolve => setTimeout(resolve, 1000));
      
      // Success - navigate to OTP verification page with email in state
      // console.log('OTP sent to:', email);
      navigate('/login-otp', { state: { email } });
    } catch (error) {
      console.error('Error sending OTP:', error);
      alert('Failed to send OTP. Please try again.');
    } finally {
      setIsSubmitting(false);
    }
  };

  return (
    <AuthLayout title="Login with OTP">
      <form onSubmit={handleSubmit} className="space-y-6 text-left">
        <div>
          <label htmlFor="email" className="block text-sm font-medium text-gray-700 dark:text-gray-300 text-left">
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
            className="mt-1 block w-full px-3 py-2 border border-gray-300 dark:border-gray-600 rounded-md shadow-sm placeholder-gray-400 dark:placeholder-gray-500 focus:outline-none focus:ring-major focus:border-major sm:text-sm bg-white/50 dark:bg-gray-800/50 text-gray-900 dark:text-white"
            placeholder="Enter your email"
            disabled={isSubmitting}
          />
        </div>
        
        <GradientButton 
          type="submit"
          isLoading={isSubmitting}
          loadingText="Sending..."
        >
          Send OTP
        </GradientButton>
        
        <p className="text-sm text-center text-gray-600 dark:text-gray-400">
          Remember your password? <Link to="/login" className="font-medium text-major hover:underline dark:text-major-light">Login with Password</Link>
        </p>
      </form>
    </AuthLayout>
  );
};

export default RequestOtpEmail;