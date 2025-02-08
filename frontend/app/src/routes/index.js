import React from 'react';
import { BrowserRouter as Router, Route, Routes } from 'react-router-dom';
import Home from '../pages/Home';
import About from '../pages/About';
import menuContent from '../contents/menuContent';
import ExploreTabContent from '../contents/tag_content/ExploreTags';
import FunctionalSection from '../components/functional_section/Functional';
import Explorer from '../pages/Explorer/Explorer';
import { sub } from 'framer-motion/client';

function ReRoutes() {
    return (
        <Routes>
            {
                menuContent.map((item) => (
                    <Route key={item.id} path={`${item.link}`} element={item.component}>
                        {
                            item.subElements && item.subElements.map((subItem) => {
                                return item.includeDocs ? <Route key={subItem.id} path={`${subItem.link}`} element={<FunctionalSection
                                    docElem={subItem.docs}
                                    funcElem={subItem.component}
                                />}>{
                                        subItem.subElements && subItem.subElements.map((subSubItem) => {
                                            return subItem.includeDocs ? <Route key={subSubItem.id} path={`${subSubItem.link}`} element={<FunctionalSection
                                                docElem={subSubItem.docs}
                                                funcElem={subSubItem.component}
                                            />} /> :
                                                <Route key={subSubItem.id} path={`${subSubItem.link}`} element={subSubItem.component} />
                                        })
                                    }
                                </Route> :
                                    <Route key={subItem.id} path={`${subItem.link}`} element={subItem.component}>
                                        {
                                            subItem.subElements && subItem.subElements.map((subSubItem) => {
                                                return subItem.includeDocs ? <Route key={subSubItem.id} path={`${subSubItem.link}`} element={<FunctionalSection
                                                    docElem={subSubItem.docs}
                                                    funcElem={subSubItem.component}
                                                />} /> :
                                                    <Route key={subSubItem.id} path={`${subSubItem.link}`} element={subSubItem.component} />
                                            })
                                        }
                                    </Route>

                            })
                        }
                    </Route>
                ))
            }
            {/* <Route key="0" path={``} element={<Explorer />}>
            <Route index element={
                    <FunctionalSection
                        docElem={ExploreTabContent[0].docs}
                        funcElem={ExploreTabContent[0].component}
                    />
                } />
                {
                    ExploreTabContent.map((tab, index) => {
                        return tab.link && ( 
                            <Route path={tab.link} element={<FunctionalSection
                                docElem={tab.docs}
                                funcElem={tab.component}
                            />} />
                        )
                    })
                }
          </Route> */}
            {/* <Route path="" element={<Explorer />} >
                <Route path="" element={<FunctionalSection
                    docElem={ExploreTabContent[0].docs}
                    funcElem={ExploreTabContent[0].component}/>} />
                <Route path="proe" element={<FunctionalSection
                    docElem={ExploreTabContent[1].docs}
                    funcElem={ExploreTabContent[1].component}/>} />
        </Route> */}
        </Routes>
    );
}

export default ReRoutes;