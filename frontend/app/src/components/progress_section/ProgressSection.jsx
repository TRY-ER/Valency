import React from 'react';
import './ProgressSection.css';
import { motion } from 'framer-motion';
import { fadeInUpVariantStatic } from '../animations/framerAnim';
import GlassyContainer from '../glassy_container/gc';

const ProgressSection = ({ processSteps }) => {
    // const currentStep = processSteps.findIndex((step) => step.status === "completed");
    const currentStep = processSteps.reduce((acc, step, idx) => {
        return step.status === 'completed' ? idx : acc;
      }, -1);
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
                    <div
                        className="progress-fill"
                        style={{
                            width: `calc(${(currentStep / (processSteps.length - 1)) * 100}% - 20px)`
                        }}
                    />
                    <div className="bullet-container">
                        {processSteps.map((step, index) => (
                            <div
                                key={index}
                                className={`progress-bullet ${step.status === "init" ?  'disabled' : ''} ${step.status === "pending" ?  'pending' : ''} ${step.status === "completed" ?  'completed' : ''} ${step.status === "failed" ?  'failed' : ''}`}
                            >
                                <div className="bullet-point" />
                            </div>
                        ))}
                    </div>
                    <div className="label-container">
                        {processSteps.map((step, index) => (
                            <span
                                className={`step-label ${step.status === "init" ?  'disabled' : ''} ${step.status === "pending" ?  'pending' : ''} ${step.status === "completed" ?  'completed' : ''}`}
                            >{step.title}</span>
                        ))}
                    </div>
                </div>
            </div>
        </GlassyContainer>
        </motion.div>
    );
};

export default ProgressSection;