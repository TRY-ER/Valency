import React, { useState } from 'react';
import AuthLayout from '../../components/Auth/AuthLayout.tsx';
import { motion } from 'framer-motion';

const MFAVerification: React.FC = () => {
  const [mfaCode, setMfaCode] = useState('');

  const handleSubmit = (e: React.FormEvent) => {
    e.preventDefault();
    // Handle MFA verification logic here
    // console.log('MFA code submitted:', { mfaCode });
    // On successful MFA, navigate to dashboard or appropriate next step
  };

  return (
    <AuthLayout title="Verify Your Identity">
      <form onSubmit={handleSubmit} className="space-y-6">
        <div>
          <label htmlFor="mfa-code" className="block text-sm font-medium text-theme-text-secondary">
            Enter the code from your authenticator app
          </label>
          <input
            id="mfa-code"
            name="mfa-code"
            type="text"
            autoComplete="one-time-code"
            required
            value={mfaCode}
            onChange={(e) => setMfaCode(e.target.value)}
            className="mt-1 block w-full px-3 py-2 border border-theme-text-secondary/50 rounded-md shadow-sm placeholder-theme-text-secondary focus:outline-none focus:ring-theme-success focus:border-theme-success sm:text-sm bg-theme-glassy text-theme-text-primary"
            maxLength={6}
          />
        </div>
        <motion.button
          whileHover={{ scale: 1.05 }}
          whileTap={{ scale: 0.95 }}
          type="submit"
          className="w-full flex justify-center py-2 px-4 border border-transparent rounded-md shadow-sm text-sm font-medium text-white bg-major hover:bg-opacity-80 focus:outline-none focus:ring-2 focus:ring-offset-2 focus:ring-major"
        >
          Verify Code
        </motion.button>
        <p className="text-sm text-center text-theme-text-secondary">
          Having trouble? <a href="#" className="font-medium text-theme-success hover:text-accent-light">Get help</a>
        </p>
      </form>
    </AuthLayout>
  );
};

export default MFAVerification;