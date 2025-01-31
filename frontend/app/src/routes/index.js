import React from 'react';
import { BrowserRouter as Router, Route, Routes} from 'react-router-dom';
import Home from '../pages/Home';
import About from '../pages/About';
import menuContent from '../contents/menuContent';

function ReRoutes() {
  return (
      <Routes>
        {
          menuContent.map((item) => (
            <Route key={item.id} path={item.link} element={item.component} />
          ))
        }
      </Routes>
  );
}

export default ReRoutes;