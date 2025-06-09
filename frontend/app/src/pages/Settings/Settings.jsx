import React from 'react';
import { useTheme } from '../../contexts/ThemeContext.js';
import { motion } from 'framer-motion';
import { fadeInUpVariants } from '../../components/animations/framerAnim.jsx';
import './Settings.css';
import CustomThemeDropdown from '../../components/UI/CustomThemeDropdown/CustomThemeDropdown.jsx'; // Import the custom dropdown

const Settings = () => {
  const { theme, setTheme } = useTheme();

  const themeOptions = [
    { value: 'light', label: 'Light' },
    { value: 'dark', label: 'Dark' },
    { value: 'system', label: 'System' },
  ];

  return (
    <motion.div 
      className="base-page-container"
      initial="hidden"
      animate="visible"
      variants={fadeInUpVariants}
    >
      <motion.h1
        variants={fadeInUpVariants}
        custom={1}
      >
        Settings
      </motion.h1>
      <motion.div 
        className="appearance-section glassy-feel"
        variants={fadeInUpVariants}
        custom={2}
      >
        <motion.h2
          variants={fadeInUpVariants}
          custom={3}
        >
          Appearance
        </motion.h2>
        <motion.div
          variants={fadeInUpVariants}
          custom={4}
        >
          <label htmlFor="theme-select">Theme: </label>
          <CustomThemeDropdown
            options={themeOptions}
            selectedValue={theme}
            onChange={(value) => setTheme(value)}
          />
        </motion.div>
      </motion.div>
      {/* Add other settings content here */}
    </motion.div>
  );
};

export default Settings;