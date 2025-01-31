import React, { useState, useRef, useEffect } from 'react';
import GlassyContainer from '../glassy_container/gc';
import "./Docs.css";
import { FaArrowDown, FaArrowUp } from "react-icons/fa";

const DocsContainer = ({ docElem }) => {
    const [isExpanded, setIsExpanded] = useState(false);
    const [showButton, setShowButton] = useState(false);
    const contentRef = useRef(null);

    const toggleExpand = () => {
        setIsExpanded(!isExpanded);
    };

    useEffect(() => {
        const contentHeight = contentRef.current?.scrollHeight; // Get the full height of the content
        const viewportHeight = window.innerHeight * 0.5; // 50vh of the viewport

        // Show the button only if the content height is greater than 50vh
        if (contentHeight > viewportHeight) {
            setShowButton(true);
        } else {
            setShowButton(false);
        }
    }, []);

    return (
        <>
        <GlassyContainer expandable={true}>
            <div
                ref={contentRef}
                className={`doc-container ${isExpanded ? 'expanded' : 'collapsed'}`}
            >
               {docElem} 
            </div>

            {/* Only show the button if content is more than 50vh */}
            <div className="doc-drop-img">
                {showButton && (
                <>
                    {isExpanded ? 
                        <FaArrowUp onClick={toggleExpand}/> : <FaArrowDown onClick={toggleExpand}/>
                    }
                </>
                )}
            </div>
            
        </GlassyContainer>

        <style jsx>{`
            .doc-container::after {
                content: '';
                position: absolute;
                bottom: 0;
                left: 0;
                width: 100%;
                height: 50px;
                background: linear-gradient(to bottom, transparent, rgba(255, 255, 255, 0.2)); /* Fade effect */
                pointer-events: none;
                opacity: ${isExpanded || !showButton ? 0 : 1}; /* Hide the fade when expanded */
                transition: opacity 0.3s ease;
            }
        `}</style>
        </>
    );
};

export default DocsContainer;