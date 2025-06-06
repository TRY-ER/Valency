import React from 'react';
import { motion } from 'framer-motion';

interface AuthLayoutProps {
  children: React.ReactNode;
  title: string;
}

const AuthLayout: React.FC<AuthLayoutProps> = ({ children, title }) => {
  return (
    <div className="min-h-screen w-full flex flex-col items-center justify-center bg-theme-bg-secondary p-4">
      <motion.div
        initial={{ opacity: 0, y: -20 }}
        animate={{ opacity: 1, y: 0 }}
        transition={{ duration: 0.5 }}
        className="w-full max-w-lg"
      >
        <div className="bg-theme-glassy shadow-glassy backdrop-blur-xl rounded-xl p-8">
          <div className="flex justify-center mb-6">
            <img src="/images/logo_rendered_main.png" alt="Valency Logo" className="w-24 h-24" />
          </div>
          <h2 className="text-3xl font-bold text-center text-theme-text-primary mb-8">{title}</h2>
          {children}
        </div>
      </motion.div>
    </div>
  );
};

export default AuthLayout;