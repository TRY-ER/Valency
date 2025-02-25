import React from 'react';
import "./Functional.css";
import DocsContainer from '../docs_section/Docs';
import FunctionContainer from '../function_container/FunctionContainer';

const FunctionalSection = ({ docElem, funcElem, customClassName = null }) => {
    return (
        <div className={`${customClassName ? customClassName : "functional-container"}`}>
            <div className="functional-item">
                <FunctionContainer functionalComponents={funcElem} />
            </div>
            {
                docElem && <div
                    className="functional-doc">
                        <DocsContainer docElem={docElem} />
                </div>
            }
        </div>
    )
}

export default FunctionalSection;