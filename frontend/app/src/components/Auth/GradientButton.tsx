import React, { useState, useEffect } from 'react';
import { motion } from 'framer-motion';

interface GradientButtonProps extends React.ButtonHTMLAttributes<HTMLButtonElement> {
  children: React.ReactNode;
  className?: string;
  isLoading?: boolean;
  loadingText?: string;
}

const GradientButton: React.FC<GradientButtonProps> = ({
  children,
  className = '',
  isLoading = false,
  loadingText = 'Processing...',
  disabled = false,
  onClick,
  type = 'button',
  ...props
}) => {
  const [isHovering, setIsHovering] = useState(false);
  const [mousePosition, setMousePosition] = useState({ x: 0, y: 0 });
  const [currentButtonStyle, setCurrentButtonStyle] = useState({});

  const handleMouseMove = (event: React.MouseEvent<HTMLButtonElement>) => {
    const rect = event.currentTarget.getBoundingClientRect();
    setMousePosition({
      x: event.clientX - rect.left,
      y: event.clientY - rect.top,
    });
  };

  useEffect(() => {
    const isDarkMode = document.documentElement.classList.contains('dark');
    let style = {};
    if (isHovering) {
      style = {
        backgroundImage: `radial-gradient(circle at ${mousePosition.x}px ${mousePosition.y}px, var(${isDarkMode ? '--color-gradient-hover-start-dark' : '--color-gradient-hover-start'}) 0%, var(${isDarkMode ? '--color-gradient-hover-end-dark' : '--color-gradient-hover-end'}) 100%)`,
      };
    } else {
      style = {
        backgroundImage: `linear-gradient(to right, var(${isDarkMode ? '--color-gradient-start-dark' : '--color-gradient-start'}), var(${isDarkMode ? '--color-gradient-end-dark' : '--color-gradient-end'}))`,
      };
    }
    setCurrentButtonStyle(style);
  }, [isHovering, mousePosition]);

  useEffect(() => {
    const observer = new MutationObserver(() => {
      const isDarkMode = document.documentElement.classList.contains('dark');
      let style = {};
      if (isHovering) {
        style = {
          backgroundImage: `radial-gradient(circle at ${mousePosition.x}px ${mousePosition.y}px, var(${isDarkMode ? '--color-gradient-hover-start-dark' : '--color-gradient-hover-start'}) 0%, var(${isDarkMode ? '--color-gradient-hover-end-dark' : '--color-gradient-hover-end'}) 100%)`,
        };
      } else {
        style = {
          backgroundImage: `linear-gradient(to right, var(${isDarkMode ? '--color-gradient-start-dark' : '--color-gradient-start'}), var(${isDarkMode ? '--color-gradient-end-dark' : '--color-gradient-end'}))`,
        };
      }
      setCurrentButtonStyle(style);
    });
    observer.observe(document.documentElement, { attributes: true, attributeFilter: ['class'] });
    return () => observer.disconnect();
  }, [isHovering, mousePosition]);

  const handleClick = (e: React.MouseEvent<HTMLButtonElement>) => {
    if (onClick && !isLoading && !disabled) {
      onClick(e);
    }
  };

  const baseClasses = "w-full flex justify-center py-2 px-4 border border-transparent rounded-md shadow-sm text-sm font-medium text-white focus:outline-none focus:ring-2 focus:ring-offset-2 focus:ring-major transition-all duration-100 ease-in-out";
  
  return (
    <motion.button
      whileTap={{ scale: 0.95 }}
      onHoverStart={() => setIsHovering(true)}
      onHoverEnd={() => setIsHovering(false)}
      onMouseMove={handleMouseMove}
      onClick={handleClick}
      type={type}
      style={currentButtonStyle}
      className={`${baseClasses} ${className}`}
      disabled={disabled || isLoading}
      {...props}
    >
      {isLoading ? loadingText : children}
    </motion.button>
  );
};

export default GradientButton;