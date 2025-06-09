import React, { useState, useEffect } from 'react';
import { useAuthContext } from '../../contexts/AuthContext.tsx';
import AuthService from '../../services/api/AuthService.ts';
import { useNavigate } from 'react-router-dom';
import { motion } from 'framer-motion';
import { fadeInUpVariants } from '../../components/animations/framerAnim.jsx';
import './Profile.css';

const Profile = () => {
  const { user, logout, isAuthenticated } = useAuthContext();
  const [currentUser, setCurrentUser] = useState(null);
  const [loading, setLoading] = useState(true);
  const [error, setError] = useState(null);
  const navigate = useNavigate();

  useEffect(() => {
    const fetchUserData = async () => {
      if (!isAuthenticated) {
        navigate('/login');
        return;
      }

      try {
        setLoading(true);
        const response = await AuthService.getCurrentUser();
        
        if (response.success && response.data) {
          setCurrentUser(response.data);
        } else {
          setError(response.error || 'Failed to fetch user data');
        }
      } catch (err) {
        setError('An error occurred while fetching user data');
        console.error('Error fetching user data:', err);
      } finally {
        setLoading(false);
      }
    };

    fetchUserData();
  }, [isAuthenticated, navigate]);

  const handleLogout = () => {
    logout();
    navigate('/login');
  };

  const formatDate = (dateString) => {
    return new Date(dateString).toLocaleDateString('en-US', {
      year: 'numeric',
      month: 'long',
      day: 'numeric',
      hour: '2-digit',
      minute: '2-digit'
    });
  };

  if (loading) {
    return (
      <motion.div 
        className="profile-loading-container"
        initial="hidden"
        animate="visible"
        variants={fadeInUpVariants}
      >
        <div className="profile-loading-content">
          <div className="profile-loading-spinner"></div>
          <p className="profile-loading-text">Loading profile...</p>
        </div>
      </motion.div>
    );
  }

  if (error) {
    return (
      <motion.div 
        className="profile-error-container"
        initial="hidden"
        animate="visible"
        variants={fadeInUpVariants}
      >
        <div className="profile-error-content">
          <div className="profile-error-card">
            <h2 className="profile-error-title">Error</h2>
            <p className="profile-error-message">{error}</p>
            <button 
              onClick={() => window.location.reload()}
              className="profile-error-retry-btn"
            >
              Retry
            </button>
          </div>
        </div>
      </motion.div>
    );
  }

  const userData = currentUser || user;

  return (
    <motion.div 
      className="profile-container"
      initial="hidden"
      animate="visible"
      variants={fadeInUpVariants}
    >
      <div className="profile-content">
        <motion.div 
          className="profile-card"
          variants={fadeInUpVariants}
          custom={1}
        >
          {/* Header */}
          <motion.div 
            className="profile-header"
            variants={fadeInUpVariants}
            custom={2}
          >
            <h1 className="profile-title">Profile</h1>
            <p className="profile-subtitle">Manage your account information</p>
          </motion.div>

          {/* User Information */}
          <motion.div 
            className="profile-info-section"
            variants={fadeInUpVariants}
            custom={3}
          >
            <h2 className="profile-info-title">User Information</h2>
            
            {userData ? (
              <div className="profile-info-list">
                <motion.div 
                  className="profile-info-item"
                  variants={fadeInUpVariants}
                  custom={4}
                >
                  <span className="profile-info-label">Username:</span>
                  <span className="profile-info-value">{userData.username}</span>
                </motion.div>
                
                <motion.div 
                  className="profile-info-item"
                  variants={fadeInUpVariants}
                  custom={5}
                >
                  <span className="profile-info-label">Email Address:</span>
                  <span className="profile-info-value">{userData.email}</span>
                </motion.div>
                
                <motion.div 
                  className="profile-info-item"
                  variants={fadeInUpVariants}
                  custom={6}
                >
                  <span className="profile-info-label">Account Created:</span>
                  <span className="profile-info-value">{formatDate(userData.created_at)}</span>
                </motion.div>
              </div>
            ) : (
              <div className="profile-no-data">
                <p className="profile-no-data-text">No user data available</p>
              </div>
            )}
          </motion.div>

          {/* Logout Button */}
          <motion.div 
            className="profile-logout-section"
            variants={fadeInUpVariants}
            custom={7}
          >
            <div className="profile-logout-container">
              <button
                onClick={handleLogout}
                className="profile-logout-btn"
              >
                <svg className="profile-logout-icon" fill="none" stroke="currentColor" viewBox="0 0 24 24">
                  <path strokeLinecap="round" strokeLinejoin="round" strokeWidth="2" d="M17 16l4-4m0 0l-4-4m4 4H7m6 4v1a3 3 0 01-3 3H6a3 3 0 01-3-3V7a3 3 0 013-3h4a3 3 0 013 3v1" />
                </svg>
                Logout
              </button>
            </div>
          </motion.div>
        </motion.div>
      </div>
    </motion.div>
  );
};

export default Profile;