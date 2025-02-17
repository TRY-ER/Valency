import React, { useState } from 'react';
import { NavLink, Outlet, useLocation } from 'react-router-dom';
import { FiDroplet } from 'react-icons/fi';  // Example icon
import { FaBars, FaArrowRight, FaArrowLeft } from 'react-icons/fa';  // Collapsible button icon
import menuContent from '../../contents/menuContent';
import { motion } from 'framer-motion';
import { fadeInLeftVariants, fadeInStatic } from '../animations/framerAnim';
import './sidebar.css';

const constructAllPaths = (item) => {
    let allPaths = [];
    if (item.link == "") {
        allPaths.push("/");
        if (item.subElements) {
            for (let i = 0; i < item.subElements.length; i++) {
                allPaths.push(`/${item.subElements[i].link}`);
                if (item.subElements[i].subElements) {
                    item.subElements[i].subElements.map((subSubItem) => {
                        if (subSubItem.link == "") {
                            return allPaths.push(`/${item.subElements[i].link}/${item.subElements[i].link}`);
                        }
                        else {
                            return allPaths.push(`/${item.subElements[i].link}/${item.subElements[i].link}/${subSubItem.link}`);
                        }
                    })
                }
            }
        }
    }
    else {
        allPaths.push(`/${item.link}`);
        // if (item.subElements) {
        //     item.subElements.map((subItem) => {
        //         allPaths.push(`/${item.link}/${subItem.link}`);
        //     })
        // }
        if (item.subElements) {
            var pre = "";
            for (let i = 0; i < item.subElements.length; i++) {
                if (item.subElements[i].link !== "") {
                    allPaths.push(`/${item.link}/${item.subElements[i].link}`);
                    pre = `/${item.link}/${item.subElements[i].link}`
                }
                else {
                    pre = `/${item.link}`
                }
                if (item.subElements[i].subElements) {
                    // console.log("sub elemems >>", item.subElements[i].subElements);
                    item.subElements[i].subElements.map((subSubItem) => {
                        if (subSubItem.link !== "") {
                            allPaths.push(`${pre}/${subSubItem.link}`);
                        }
                    })
                }
            }
        }
    }

    // console.log("all paths >>", allPaths);

    return allPaths;
}

const Sidebar = ({ }) => {
    const [isCollapsed, setIsCollapsed] = useState(false);
    const location = useLocation();

    const toggleSidebar = () => {
        setIsCollapsed(!isCollapsed);
    };

    return (
        <div className={`sidebar ${isCollapsed ? 'collapsed' : ''}`}>
            <motion.div
                variants={fadeInStatic}
                initial="hidden"
                animate="visible"
                className={`logo ${isCollapsed ? 'collapsed' : ''}`}>
                {/* Placeholder for Logo */}
                <img src="/images/valency_logo_light_600x600.png" alt="Logo" />
            </motion.div>
            {
                isCollapsed ? <>
                    <motion.div
                        variants={fadeInStatic}
                        initial="hidden"
                        animate="visible"
                        className='temp-arrow-sec'
                    >
                        <FaArrowRight color={"white"} className={`toggle-btn ${isCollapsed ? 'collapsed' : ''}`} onClick={toggleSidebar} />
                    </motion.div>
                </>
                    :
                    <>
                        <motion.div
                            variants={fadeInStatic}
                            initial="hidden"
                            animate="visible"
                            className='temp-arrow-sec'
                        >
                            <FaArrowLeft color={"white"} className={`toggle-btn ${isCollapsed ? 'collapsed' : ''}`} onClick={toggleSidebar} />
                        </motion.div>
                    </>
            }

            <div className={`sidebar-menu-container ${isCollapsed ? 'collapsed' : ''}`}>
                {
                    menuContent.map((item, index) => {
                        const allNestPaths = constructAllPaths(item);
                        return <NavLink key={item.id}
                            //   to={`/${item.link}`} 
                            to={item.link === "" ? "/" : `/${item.link}`}
                            // onClick={() => setActiveMenu(item.id)} 
                            // className={({ isActive }) => `sidebar-menu-item ${isActive ? 'active' : ''}`}
                            className={() =>
                                allNestPaths.includes(location.pathname)
                                    ? "sidebar-menu-item active"
                                    : "sidebar-menu-item"
                            }
                        >
                            <motion.div
                                variants={fadeInLeftVariants}
                                initial="hidden"
                                animate="visible"
                                custom={index}
                                className={`sidebar-menu-content ${isCollapsed ? 'collapsed' : ''}`}>
                                <img src={item.iconPath} alt="Icon" />
                                {
                                    isCollapsed ? null : <>
                                        <motion.h2
                                            variants={fadeInLeftVariants}
                                            initial="hidden"
                                            animate="visible"
                                            custom={index}
                                        >{item.title}</motion.h2>
                                    </>
                                }
                            </motion.div>
                        </NavLink>
                    })
                }
            </div>
        </div>
    );
};

export default Sidebar;
