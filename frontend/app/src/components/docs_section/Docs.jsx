import React from 'react';
import GlassyContainer from '../glassy_container/gc';
import "./Docs.css";

const DocsContainer = ({ docElem }) => {
    return (
        <>
            {/* <motion.div
                initial="hidden"
                animate="visible"
                variants={fadeInRightVariantStatic}> */}
                <GlassyContainer expandable={true}>
                    <div
                        className={`doc-container`}
                    >
                        {docElem}
                    </div>
                </GlassyContainer>
        </>
    );
};

export default DocsContainer;