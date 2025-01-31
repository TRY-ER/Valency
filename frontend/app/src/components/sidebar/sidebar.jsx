import React, { useState } from 'react';
import { Link } from 'react-router-dom';
import { FiDroplet } from 'react-icons/fi';  // Example icon
import { FaBars, FaArrowRight, FaArrowLeft } from 'react-icons/fa';  // Collapsible button icon
import menuContent from '../../contents/menuContent';
import './sidebar.css';

const Sidebar = ({activeMenu, setActiveMenu}) => {
    const [isCollapsed, setIsCollapsed] = useState(false);

    const toggleSidebar = () => {
        setIsCollapsed(!isCollapsed);
    };

    return (
        <div className={`sidebar ${isCollapsed ? 'collapsed' : ''}`}>
            <div className={`logo ${isCollapsed ? 'collapsed' : ''}`}>
                {/* Placeholder for Logo */}
                <img src="images/valency_logo_light_600x600.png" alt="Logo" />
            </div>
            {
                isCollapsed ? <>
                    <FaArrowRight color={"white"} className={`toggle-btn ${isCollapsed ? 'collapsed' : ''}`} onClick={toggleSidebar} />
                </>
                    :
                    <>
                        <FaArrowLeft color={"white"} className={`toggle-btn ${isCollapsed ? 'collapsed' : ''}`} onClick={toggleSidebar} />
                    </>
            }

            <div className={`sidebar-menu-container ${isCollapsed ? 'collapsed' : ''}`}>
                {
                    menuContent.map((item) => (
                        <Link key={item.id} 
                              to={item.link} 
                              className="sidebar-menu-item"
                                onClick={() => setActiveMenu(item.id)} 
                              >
                            <div className={`sidebar-menu-content ${isCollapsed ? 'collapsed' : ''} ${activeMenu === item.id ? 'active' : ''}`}>
                                <img src={item.iconPath} alt="Icon" />
                                {
                                    isCollapsed ? null : <>
                                        <h2>{item.title}</h2>
                                    </>
                                }
                            </div>
                        </Link>
                    ))
                }

            </div>
        </div>
    );
};

export default Sidebar;
