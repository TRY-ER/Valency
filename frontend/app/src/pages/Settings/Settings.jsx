import React from 'react';
import { useTheme } from '../../contexts/ThemeContext.js';
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
    <div className="base-page-container">
      <h1>Settings</h1> {/* Changed title here */}
      <div className="appearance-section glassy-feel">
        <h2>Appearance</h2>
        <div>
          <label htmlFor="theme-select">Theme: </label>
          <CustomThemeDropdown
            options={themeOptions}
            selectedValue={theme}
            onChange={(value) => setTheme(value)}
          />
        </div>
      </div>
      {/* Add other settings content here */}
    </div>
  );
};

export default Settings;