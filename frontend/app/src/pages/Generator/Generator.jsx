import React  from 'react';
import './Generator.css';
import Divider from '../../components/divider';
// import { useLocation, useNavigate } from 'react-router-dom';
import { useLocation, Outlet, NavLink} from 'react-router-dom';
import { motion } from 'framer-motion';
import { fadeInDownVariants } from '../../components/animations/framerAnim';

export default function Generator({
    tabContent,
    basePath = "generators"
}) {

    const location = useLocation();

    const constructNestPath = (tab) => {
        let allPaths = [];
        if (tab.link === "") {
            allPaths.push(`/${basePath}`);
            if (tab.subElements) {
                for (let subTab of tab.subElements) {
                    if (subTab.link === "") {
                        allPaths.push(`/${basePath}`);
                    }
                    else {
                        allPaths.push(`/${basePath}/${subTab.link}`);
                    }
                }
            }
        } else {
            allPaths.push(`/${basePath}/${tab.link}`);
            if (tab.subElements) {
                for (let subTab of tab.subElements) {
                    if (subTab.link === "") {
                        allPaths.push(`/${basePath}/${tab.link}`);
                    }
                    else {
                        allPaths.push(`/${basePath}/${tab.link}/${subTab.link}`);
                    }
                }
            }
        }
        // console.log("disc nest paths >>", allPaths);
        return allPaths;
    }

    return (
        <div className="base-page-container">
             <div className="tab-container">
                {
                    tabContent.map((tab, index) => {
                        const allNestPaths = constructNestPath(tab);
                        // console.log("current path >>", location.pathname);
                        // console.log("allNestPaths >>", allNestPaths);
                        // console.log("bool res >>", allNestPaths.includes(location.pathname));
                        return (
                            <NavLink
                                key={tab.id}
                                to={tab.link === "" ? `/${basePath}` : `/${basePath}/${tab.link}`}
                                className={() =>
                                    allNestPaths.includes(location.pathname)
                                        ? "tab-link active"
                                        : "tab-link"
                                }
                            >
                                <motion.div 
                                initial="hidden"
                                animate="visible"
                                variants={fadeInDownVariants}
                                custom={index}  
                                className={`tab-tag glassy-feel`}>
                                    <p className="tab-tag-text">{tab.title}</p>
                                </motion.div>
                            </NavLink>
                        )
                    })
                }
            </div>
            <Divider />
            <Outlet />
            {/* <Routes>
                {
                    GeneratorTabContents.map((tab, index) => {
                        return (
                            <Route path={`generators${tab.link}`} element={<FunctionalSection
                                docElem={tab.docs}
                                funcElem={tab.component}
                            />} />
                        )
                    })
                }
            </Routes> */}
            {/* <FunctionalSection 
                docElem={GeneratorTabContents[activeTabId - 1].docs} 
                funcElem={GeneratorTabContents[activeTabId - 1].component}
            /> */}
        </div>
    )
}