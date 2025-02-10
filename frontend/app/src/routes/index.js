import React from 'react';
import { BrowserRouter as Router, Route, Routes } from 'react-router-dom';
import menuContent from '../contents/menuContent';
import FunctionalSection from '../components/functional_section/Functional';

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
                                                <Route key={subSubItem.id} path={`${subSubItem.link}`} element={<FunctionalSection
                                                    funcElem={subSubItem.component}
                                                />} />
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
                                                    <Route key={subSubItem.id} path={`${subSubItem.link}`} element={<FunctionalSection
                                                        funcElem={subSubItem.component}
                                                    />} />
                                            })
                                        }
                                    </Route>

                            })
                        }
                    </Route>
                ))
            }
        </Routes>
    );
}

export default ReRoutes;