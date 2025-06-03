import React, { useState, useEffect } from 'react';
import './ActivitiesPage.css';
import { motion } from 'framer-motion';
import { useNavigate } from 'react-router-dom';
import { 
    FaClock, FaCheckCircle, FaExclamationTriangle, FaHourglass, FaSpinner, 
    FaSortAmountDown, FaSortAmountUp, FaQuestionCircle, FaRobot, 
    FaCodeBranch, FaTrashAlt
} from 'react-icons/fa';
import { listUserSessions, deleteSession } from '../../services/api/agentService';

const ActivitiesPage = () => {
    const [activities, setActivities] = useState([]);
    const [isLoading, setIsLoading] = useState(true);
    const [sortDirection, setSortDirection] = useState('desc'); // 'desc' for newest first, 'asc' for oldest first
    const navigate = useNavigate(); // Hook for navigation
    const [feedback, setFeedback] = useState({ message: '', type: '' }); // For success/error messages

    // Fetch actual user sessions
    useEffect(() => {
        setIsLoading(true);
        listUserSessions()
            .then((response) => {
                console.log('Fetched activities:', response);
                // Convert session data to activity format
                const formattedActivities = Array.isArray(response) ? 
                    response.map((session, index) => {
                        // Determine status based on activity counts
                        return {
                            id: session.session_id || index,
                            name: session.session_identifier || 'Unnamed Session',
                            timestamp: session.last_update_time,
                            type: 'session',
                            rawData: session
                        };
                    }) : [];
                
                setActivities(formattedActivities);
                setIsLoading(false);
            })
            .catch(error => {
                console.error('Error fetching sessions:', error);
                setIsLoading(false);
            });
    }, []);

    // Sort activities by timestamp
    const sortedActivities = [...activities].sort((a, b) => {
        const dateA = new Date(a.timestamp);
        const dateB = new Date(b.timestamp);
        return sortDirection === 'asc' ? dateA - dateB : dateB - dateA;
    });

    const toggleSortDirection = () => {
        setSortDirection(sortDirection === 'desc' ? 'asc' : 'desc');
    };

    const formatDate = (dateString) => {
        if (!dateString) return 'Unknown date';
        
        // Handle the format "2025-06-03 06:34:04"
        const options = {
            year: 'numeric',
            month: 'short',
            day: 'numeric',
            hour: '2-digit',
            minute: '2-digit'
        };
        
        try {
            return new Date(dateString).toLocaleString(undefined, options);
        } catch (error) {
            console.error('Error formatting date:', error);
            return dateString; // Return original string if parsing fails
        }
    };

    // Handle deletion of a session
    const handleDeleteSession = (sessionId, event) => {
        event.stopPropagation(); // Prevent event from bubbling up
        console.log('Delete session requested for ID:', sessionId);
        
        // Confirm before deletion
        if (window.confirm('Are you sure you want to delete this session? This action cannot be undone.')) {
            setIsLoading(true); // Show loading state
            
            deleteSession(sessionId)
                .then((response) => {
                    console.log('Session deleted successfully:', response);
                    // Update the UI by removing the deleted session
                    setActivities((prevActivities) => 
                        prevActivities.filter((activity) => activity.id !== sessionId)
                    );
                    // Show success message
                    setFeedback({
                        message: 'Session deleted successfully!',
                        type: 'success'
                    });
                    
                    // Clear feedback message after 3 seconds
                    setTimeout(() => {
                        setFeedback({ message: '', type: '' });
                    }, 3000);
                })
                .catch((error) => {
                    console.error('Error deleting session:', error);
                    // Show error message
                    setFeedback({
                        message: 'Failed to delete session. Please try again later.',
                        type: 'error'
                    });
                    
                    // Clear feedback message after 3 seconds
                    setTimeout(() => {
                        setFeedback({ message: '', type: '' });
                    }, 3000);
                })
                .finally(() => {
                    setIsLoading(false); // Hide loading state
                });
        }
    };

    // Handle click on a session card
    const handleCardClick = (sessionId) => {
        navigate(`/chatbot/?session_id=${sessionId}`);
    };

    return (
        <div className="activities-page">
            <div className="activities-header">
                <h1>Session Activities Dashboard</h1>
                <p>Track and manage your scientific activities here</p>

                {feedback.message && (
                    <div className={`feedback-message ${feedback.type}`}>
                        {feedback.message}
                    </div>
                )}

                <div className="sort-control">
                    <button
                        className="sort-btn"
                        onClick={toggleSortDirection}
                    >
                        {sortDirection === 'desc' ? (
                            <>
                                <FaSortAmountDown className="sort-icon" />
                                <span>Newest First</span>
                            </>
                        ) : (
                            <>
                                <FaSortAmountUp className="sort-icon" />
                                <span>Oldest First</span>
                            </>
                        )}
                    </button>
                </div>
            </div>

            <div className="activities-content">
                {isLoading ? (
                    <div className="loading-container">
                        <div className="spinner"></div>
                        <p>Loading activities...</p>
                    </div>
                ) : sortedActivities.length === 0 ? (
                    <div className="no-activities">
                        <p>No activities found.</p>
                    </div>
                ) : (
                    <div className="activities-list">
                        {sortedActivities.map((activity) => (
                            <motion.div
                                key={activity.id}
                                className="activity-card"
                                initial={{ opacity: 0, y: 20 }}
                                animate={{ opacity: 1, y: 0 }}
                                transition={{ duration: 0.3 }}
                                onClick={() => handleCardClick(activity.id)}
                                style={{ cursor: 'pointer' }}
                            >
                                <div className="activity-header">
                                    <h3>{activity.name}</h3>
                                    <div className="delete-action">
                                        <button 
                                            className="delete-btn"
                                            onClick={(e) => handleDeleteSession(activity.id, e)}
                                            title="Delete session"
                                        >
                                            <FaTrashAlt className="delete-icon" />
                                        </button>
                                    </div>
                                </div>
                                <p className="metrics-label">Session Statistics</p>
                                <div className="activity-metrics">
                                    <div className="metric-item">
                                        <FaQuestionCircle className="metric-icon" />
                                        <span className="metric-count">{activity.rawData.query_count || 0}</span>
                                        <div className="metric-tooltip">Queries submitted</div>
                                    </div>
                                    <div className="metric-item">
                                        <FaRobot className="metric-icon" />
                                        <span className="metric-count">{activity.rawData.agent_call_count || 0}</span>
                                        <div className="metric-tooltip">AI agent interactions</div>
                                    </div>
                                    <div className="metric-item">
                                        <FaCodeBranch className="metric-icon" />
                                        <span className="metric-count">{activity.rawData.function_response_count || 0}</span>
                                        <div className="metric-tooltip">Function responses</div>
                                    </div>
                                </div>
                                <div className="activity-footer">
                                    <span className="activity-session-id">
                                        ID: {typeof activity.id === 'string' && activity.id.length > 8 
                                            ? `${activity.id.slice(0, 8)}...` 
                                            : activity.id}
                                    </span>
                                    <span className="activity-time">{formatDate(activity.timestamp)}</span>
                                </div>
                            </motion.div>
                        ))}
                    </div>
                )}
            </div>
        </div>
    );
};

export default ActivitiesPage;
