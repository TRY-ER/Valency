import React from 'react';
import './ProgressLoader.css';
import { motion } from 'framer-motion';
import { fadeInUpVariantStatic } from '../animations/framerAnim';
import GlassyContainer from '../glassy_container/gc';

const ProgressLoader = ({ progressNum, maxNum, status }) => {
    return (
        <motion.div
            variants={fadeInUpVariantStatic}
            initial="hidden"
            animate="visible"
        >
            <GlassyContainer>
                <div className="progress-container">
                    <div className="progress-track">
                        <div className="progress-unfil"></div>
                        {
                            maxNum === 0 ?
                                <div
                                    className="progress-fill"
                                    style={{
                                        width: `0%`
                                    }} /> :
                                <div
                                    className="progress-fill"
                                    style={{
                                        width: `calc(${(progressNum / maxNum) * 100}% - 20px)`
                                    }} />
                        }
                    </div>
                    <div className="progress-label">
                        <p className='progress-label-left'>Status: {status}</p>
                        {
                            progressNum > 0 && <p className='progress-label-right'>
                                Progress: {(progressNum / maxNum) * 100}% ( {progressNum} / {maxNum})
                            </p>
                        }

                    </div>
                </div>
            </GlassyContainer>
        </motion.div>
    );
};

export default ProgressLoader;