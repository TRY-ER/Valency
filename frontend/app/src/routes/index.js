import React from 'react';
import { Route, Routes } from 'react-router-dom';
import menuContent from '../contents/menuContent';
import FunctionalSection from '../components/functional_section/Functional';

// Auth components
import ProtectedRoute from '../components/Auth/ProtectedRoute.tsx';
import AuthRedirectRoute from '../components/Auth/AuthRedirectRoute.tsx';

// Page components
import Home from '../pages/Home'; // Main application home page
import About from '../pages/About'; // Public about page

// Auth Page Components
import Login from '../pages/Auth/Login.tsx';
import Signup from '../pages/Auth/Signup.tsx';
import VerifyEmail from '../pages/Auth/VerifyEmail.tsx';
import ForgotPassword from '../pages/Auth/ForgotPassword.tsx';
import ResetPassword from '../pages/Auth/ResetPassword.tsx';
import LoginOTP from '../pages/Auth/LoginOTP.tsx';
import RequestOtpEmail from '../pages/Auth/RequestOtpEmail.tsx';
import MFAVerification from '../pages/Auth/MFAVerification.tsx';
import NotFoundRedirect from '../components/Auth/NotFoundRedirect.tsx'; // Import the new component
import Profile from '../pages/Profile/Profile'; // Added
import Settings from '../pages/Settings/Settings'; // Added

// Protected Layout
import ProtectedLayout from '../components/layout/ProtectedLayout';


//default to explorer
import MolEComponent from '../pages/Explorer/MolExplorer/MolExplorer.jsx';
import DocRenderer from '../contents/doc_content/DocRenderer.jsx';


function ReRoutes() {
    return (
        <Routes>
            {/* Authentication routes (e.g., login, signup, verify, password reset) */}
            {/* Accessible only if the user is NOT authenticated. */}
            {/* If authenticated, AuthRedirectRoute redirects to '/' */}
            <Route element={<AuthRedirectRoute />}>
                <Route path="/login" element={<Login />} />
                <Route path="/signup" element={<Signup />} />
                <Route path="/verify-email" element={<VerifyEmail />} /> {/* Token is read from query params in the component */}
                <Route path="/forgot-password" element={<ForgotPassword />} />
                <Route path="/reset-password" element={<ResetPassword />} /> {/* Token is read from query params in the component */}
                <Route path="/request-otp-email" element={<RequestOtpEmail />} />
                <Route path="/login-otp" element={<LoginOTP />} />
                {/* MFAVerification is typically part of the protected flow or post-login, not under AuthRedirectRoute */}
            </Route>

            {/* Protected application routes */}
            {/* Accessible only if the user IS authenticated. */}
            {/* If not authenticated, ProtectedRoute redirects to '/login'. */}
            {/* ProtectedLayout includes the Sidebar and an Outlet for the page content. */}
            <Route element={<ProtectedRoute />}>
                {/* MFA Verification route - should be accessed after login if MFA is required */}
                <Route path="/mfa-verification" element={<MFAVerification />} /> 
                <Route element={<ProtectedLayout />}>
                    {
                        menuContent.map((item) => (
                            <Route key={item.id} path={`${item.link}`} element={item.component}>
                                {
                                    item.subElements && item.subElements.map((subItem) => {
                                        return item.includeDocs ? (
                                            <Route key={subItem.id} path={`${subItem.link}`} element={<FunctionalSection
                                                docElem={subItem.docs}
                                                funcElem={subItem.component}
                                            />}>
                                                {
                                                    subItem.subElements && subItem.subElements.map((subSubItem) => {
                                                        return subItem.includeDocs ? (
                                                            <Route key={subSubItem.id} path={`${subSubItem.link}`} element={<FunctionalSection
                                                                docElem={subSubItem.docs}
                                                                funcElem={subSubItem.component}
                                                            />} />
                                                        ) : (
                                                            <Route key={subSubItem.id} path={`${subSubItem.link}`} element={<FunctionalSection
                                                                funcElem={subSubItem.component}
                                                            />} />
                                                        );
                                                    })
                                                }
                                            </Route>
                                        ) : (
                                            <Route key={subItem.id} path={`${subItem.link}`} element={subItem.component}>
                                                {
                                                    subItem.subElements && subItem.subElements.map((subSubItem) => {
                                                        return subItem.includeDocs ? (
                                                            <Route key={subSubItem.id} path={`${subSubItem.link}`} element={<FunctionalSection
                                                                docElem={subSubItem.docs}
                                                                funcElem={subSubItem.component}
                                                            />} />
                                                        ) : (
                                                            <Route key={subSubItem.id} path={`${subSubItem.link}`} element={<FunctionalSection
                                                                funcElem={subSubItem.component}
                                                            />} />
                                                        );
                                                    })
                                                }
                                            </Route>
                                        );
                                    })
                                }
                            </Route>
                        ))
                    }
                    {/* Add other top-level protected routes here if not in menuContent */}
                    <Route path="/profile" element={<Profile />} /> {/* Added */}
                    <Route path="/settings" element={<Settings />} /> {/* Added */}
                </Route>
            </Route>

            {/* Publicly accessible routes */}
            <Route path="/about" element={<About />} />
            {/* Add other public routes like a general landing page here if needed */}

            {/* Catch-all route for 404 Not Found - uses NotFoundRedirect */}
            <Route path="*" element={<NotFoundRedirect />} />
        </Routes>
    );
}

export default ReRoutes;