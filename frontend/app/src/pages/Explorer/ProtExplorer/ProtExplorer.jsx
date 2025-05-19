import React from "react"; // Removed useContext
import { motion } from "framer-motion";
import { fadeInUpVariantStatic } from "../../../components/animations/framerAnim";
import Divider from "../../../components/divider";
import { useLocation, Outlet, NavLink } from 'react-router-dom';
import { fadeInDownVariants } from "../../../components/animations/framerAnim";
// Removed ThemeContext import

const ProtEComponent = ({
    tabContent, basePath = ""
}) => {
    const location = useLocation();
    // Removed theme consumption from ThemeContext

    return (
        <>
            <div className="tab-container">
                {
                    tabContent.map((tab, index) => {
                        return (
                           <NavLink
                                key={tab.id}
                                to={ tab.link === "" ? `/${basePath}` : `/${basePath}/${tab.link}`}
                                className={() =>{
                                    if (location.pathname === `/${basePath}` && tab.link === ""){
                                        return `tab-link active`
                                    }
                                    else if(location.pathname === `/${basePath}/${tab.link}`){
                                        return `tab-link active`
                                    }
                                    else{
                                        return `tab-link`
                                    }
                                }} 
                            >
                                <motion.div 
                                    initial="hidden"
                                    animate="visible"
                                    variants={fadeInDownVariants} 
                                    custom={index}
                                    className="tab-tag glassy-feel"> {/* Reverted className */}
                                    <p className="tab-tag-text">{tab.title}</p>
                                </motion.div>
                            </NavLink>
                        )
                    })
                }
            </div>
            <Divider />
            <Outlet /> 
        </>
    ) 
}


export default ProtEComponent;