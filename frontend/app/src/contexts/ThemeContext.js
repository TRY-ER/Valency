import React, { createContext, useContext, useEffect, useState } from 'react';

const ThemeContext = createContext();

export const useTheme = () => useContext(ThemeContext);

export const ThemeProvider = ({ children }) => {
  const [theme, setTheme] = useState(() => {
    const storedTheme = localStorage.getItem('theme');
    return storedTheme || 'system';
  });

  useEffect(() => {
    localStorage.setItem('theme', theme);
    applyTheme();
  }, [theme]);

  const applyTheme = () => {
    const root = window.document.documentElement;
    const isDarkQuery = window.matchMedia('(prefers-color-scheme: dark)');

    const applyActualTheme = (actualTheme) => {
      if (actualTheme === 'dark') {
        root.classList.add('dark');
      } else {
        root.classList.remove('dark');
      }
    };

    if (theme === 'system') {
      applyActualTheme(isDarkQuery.matches ? 'dark' : 'light');
      const mediaQueryListener = (e) => applyActualTheme(e.matches ? 'dark' : 'light');
      isDarkQuery.addEventListener('change', mediaQueryListener);
      return () => isDarkQuery.removeEventListener('change', mediaQueryListener);
    } else {
      applyActualTheme(theme);
    }
  };
  
  // Initial theme application on load
  useEffect(() => {
    applyTheme();
  }, []);


  return (
    <ThemeContext.Provider value={{ theme, setTheme }}>
      {children}
    </ThemeContext.Provider>
  );
};
