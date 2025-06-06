import React, { useState } from 'react';
import { NavLink, useLocation } from 'react-router-dom';
import menuContent from '../../contents/menuContent';
import { motion } from 'framer-motion';
import { fadeInLeftVariants, fadeInStatic } from '../animations/framerAnim';
import './sidebar.css';
import { FaAngleLeft, FaAngleRight, FaUserCircle, FaCog } from 'react-icons/fa';

const Sidebar = ({ }) => {
    const [isCollapsed, setIsCollapsed] = useState(false);
    const location = useLocation();

    const toggleSidebar = () => {
        setIsCollapsed(!isCollapsed);
    };

    return (
        <div className={`sidebar ${isCollapsed ? 'collapsed' : ''}`}>
            <div className={`toggle-button-container ${isCollapsed ? 'collapsed' : ''}`} onClick={toggleSidebar}>
                {isCollapsed ? <FaAngleRight className="toggle-btn-icon" /> : <FaAngleLeft className="toggle-btn-icon" />}
            </div>
            <motion.div
                variants={fadeInStatic}
                initial="hidden"
                animate="visible"
                className={`logo ${isCollapsed ? 'collapsed' : ''}`}>
                {/* Placeholder for Logo */}
                <img src="/images/logo_rendered_main.png" alt="Logo" />
                {!isCollapsed && (
                    <motion.h1
                        variants={fadeInLeftVariants}
                        initial="hidden"
                        animate="visible"
                        custom={0}
                        className="logo-text"
                    >
                        VALENCY
                    </motion.h1>
                )}
            </motion.div>

            <div className={`sidebar-menu-container ${isCollapsed ? 'collapsed' : ''}`}>
                {
                    menuContent.map((item, index) => {
                        const itemPath = item.link === "" ? "/" : `/${item.link}`;

                        return <NavLink key={item.id}
                            to={itemPath}
                            end={itemPath === "/"} // Exact match for root, prefix match for others.
                            className={({ isActive: navLinkIsActive }) => {
                                return navLinkIsActive ? "sidebar-menu-item active" : "sidebar-menu-item";
                            }}
                        >
                            <motion.div
                                variants={fadeInLeftVariants}
                                initial="hidden"
                                animate="visible"
                                custom={index}
                                className={`sidebar-menu-content ${isCollapsed ? 'collapsed' : ''}`}>
                                <span className='sidebar-icon'>{item.icon}</span>
                                {
                                    isCollapsed ? (
                                        <span className="sidebar-tooltip">{item.title}</span>
                                    ) : (
                                        <motion.h2
                                            variants={fadeInLeftVariants}
                                            initial="hidden"
                                            animate="visible"
                                            custom={index}
                                        >{item.title}</motion.h2>
                                    )
                                }
                            </motion.div>
                        </NavLink>
                    })
                }
            </div>

            {/* Added bottom navigation icons */}
            <div className={`sidebar-bottom-nav ${isCollapsed ? 'collapsed' : ''}`}>
                <NavLink
                    to="/profile"
                    className={({ isActive }) => isActive ? "sidebar-bottom-nav-item active" : "sidebar-bottom-nav-item"}
                    title="Profile"
                >
                    <span className='sidebar-icon'><FaUserCircle /></span>
                    {!isCollapsed && <span className="sidebar-bottom-nav-text">Profile</span>}
                    {isCollapsed && <span className="sidebar-tooltip">Profile</span>}
                </NavLink>
                <NavLink
                    to="/settings"
                    className={({ isActive }) => isActive ? "sidebar-bottom-nav-item active" : "sidebar-bottom-nav-item"}
                    title="Settings"
                >
                    <span className='sidebar-icon'><FaCog /></span>
                    {!isCollapsed && <span className="sidebar-bottom-nav-text">Settings</span>}
                    {isCollapsed && <span className="sidebar-tooltip">Settings</span>}
                </NavLink>
            </div>
        </div>
    );
};

export default Sidebar;
