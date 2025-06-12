import React from 'react';
import { Route, Routes, Navigate } from 'react-router-dom';
import menuContent from '../contents/menuContent';
import FunctionalSection from '../components/functional_section/Functional';

// Auth components
import ProtectedRoute from '../components/Auth/ProtectedRoute.tsx';
import VerifiedRoute from '../components/Auth/VerifiedRoute.tsx';
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
                <Route path="/forgot-password" element={<ForgotPassword />} />
                <Route path="/reset-password" element={<ResetPassword />} /> {/* Token is read from query params in the component */}
                <Route path="/request-otp-email" element={<RequestOtpEmail />} />
                <Route path="/login-otp" element={<LoginOTP />} />
                {/* MFAVerification is typically part of the protected flow or post-login, not under AuthRedirectRoute */}
            </Route>

            {/* Protected application routes */}
            {/* Accessible only if the user IS authenticated. */}
            {/* If not authenticated, ProtectedRoute redirects to '/login'. */}
            <Route element={<ProtectedRoute />}>
                {/* MFA Verification route - should be accessed after login if MFA is required */}
                <Route path="/mfa-verification" element={<MFAVerification />} />
                {/* Email verification route - should be accessible to authenticated but unverified users */}
                <Route path="/verify-email" element={<VerifyEmail />} />
                
                {/* Routes that require both authentication AND email verification */}
                <Route element={<VerifiedRoute />}>
                    <Route element={<ProtectedLayout />}>
                        {
                            menuContent.map((item) => {
                                // Handle the root route (Home page with empty link)
                                if (item.link === '') {
                                    return (
                                        <Route key={item.id} index element={item.component} />
                                    );
                                }
                                
                                // Determine if the current item is the 'Structure Analysis' (explorer) section
                                const isExplorer = item.link === 'explorer';
                                // For the explorer route, we want its direct component to render, and its subElements to define nested routes.
                                // The Explorer component itself contains an <Outlet/> for these nested routes.
                                const elementToRender = item.component;

                                return (
                                    <Route key={item.id} path={`${item.link}`} element={elementToRender}>
                                        {/* Nested routes for sub-elements (e.g., tabs within Explorer) */}
                                        {item.subElements && item.subElements.map((subItem) => {
                                            const subElementPath = subItem.link; // subItem.link is already relative
                                            const subElementComponent = subItem.includeDocs ? (
                                                <FunctionalSection
                                                    docElem={subItem.docs}
                                                    funcElem={subItem.component}
                                                />
                                            ) : subItem.component;

                                            return (
                                                <Route key={subItem.id} path={subElementPath} element={subElementComponent}>
                                                    {/* Handling third-level nesting if necessary (e.g., Protein Explorer's own tabs) */}
                                                    {subItem.subElements && subItem.subElements.map((subSubItem) => {
                                                        const subSubElementPath = subSubItem.link;
                                                        const subSubElementComponent = subSubItem.includeDocs ? (
                                                            <FunctionalSection
                                                                docElem={subSubItem.docs}
                                                                funcElem={subSubItem.component}
                                                            />
                                                        ) : subSubItem.component;
                                                        return (
                                                            <Route key={subSubItem.id} path={subSubElementPath} element={subSubElementComponent} />
                                                        );
                                                    })}
                                                </Route>
                                            );
                                        })}
                                        {/* If it's the explorer route and it has subElements, 
                                            create an index route for the first subElement if its link is empty. 
                                            This makes /explorer default to its first tab.*/}
                                        {isExplorer && item.subElements && item.subElements.find(se => se.link === '') && (
                                            <Route 
                                                index 
                                                element={(
                                                    item.subElements.find(se => se.link === '').includeDocs ? (
                                                        <FunctionalSection
                                                            docElem={item.subElements.find(se => se.link === '').docs}
                                                            funcElem={item.subElements.find(se => se.link === '').component}
                                                        />
                                                    ) : item.subElements.find(se => se.link === '').component
                                                )}
                                            />
                                        )}
                                    </Route>
                                );
                            })
                        }
                        {/* Add other top-level protected routes here if not in menuContent */}
                        <Route path="/profile" element={<Profile />} />
                        <Route path="/settings" element={<Settings />} />
                    </Route>
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