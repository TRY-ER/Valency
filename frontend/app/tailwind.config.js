/** @type {import('tailwindcss').Config} */
export default {
  content: [
    "./index.html",
    "./src/**/*.{js,ts,jsx,tsx}",
  ],
  darkMode: 'class', // or 'media' if you prefer system setting directly
  theme: {
    extend: {
      colors: {
        // Core theme colors from App.css
        'theme-bg-primary': 'var(--color-bg-primary)', // black
        'theme-bg-secondary': 'var(--color-bg-secondary)', // rgb(36, 36, 36)
        'theme-text-primary': 'var(--color-text-primary)', // white
        'theme-text-secondary': 'var(--color-text-secondary)', // rgb(179, 179, 179)
        'theme-success': 'var(--color-success)', // rgb(5, 173, 5)
        'theme-alert': 'var(--color-alert)', // red
        'theme-glassy': 'var(--glassy-color)', // rgba(255, 255, 255, 0.2)

        // Accent colors based on theme-success (green)
        accent: 'var(--color-success)', // rgb(5, 173, 5)
        'accent-light': 'rgb(10, 200, 10)', // Lighter green for hover/gradients
        'accent-dark': 'rgb(0, 150, 0)',   // Darker green for active/gradients

        // Gradient specific colors - now green based
        'gradient-start': 'var(--color-accent-light)', 
        'gradient-end': 'var(--color-accent)',   
        'gradient-start-dark': 'var(--color-accent-light)', // Using same for dark mode, adjust if needed
        'gradient-end-dark': 'var(--color-accent)',   
        'gradient-hover-start': 'var(--color-accent)', 
        'gradient-hover-end': 'var(--color-accent-dark)',   
        'gradient-hover-start-dark': 'var(--color-accent)', 
        'gradient-hover-end-dark': 'var(--color-accent-dark)',
      },
      boxShadow: {
        'glassy': '0 4px 30px rgba(0, 0, 0, 0.1)', // from .glassy-feel in App.css
      },
      backgroundImage: {
        // Gradients using the new accent (green) colors
        'button-gradient-light': 'linear-gradient(to right, theme(colors.gradient-start), theme(colors.gradient-end))',
        'button-gradient-dark': 'linear-gradient(to right, theme(colors.gradient-start-dark), theme(colors.gradient-end-dark))',
        'button-gradient-light-hover': 'linear-gradient(to right, theme(colors.gradient-hover-start), theme(colors.gradient-hover-end))',
        'button-gradient-dark-hover': 'linear-gradient(to right, theme(colors.gradient-hover-start-dark), theme(colors.gradient-hover-end-dark))',
      }
    },
  },
  plugins: [],
}