import React, { useState, useEffect } from "react";
import "./ADMET.css";
import MolInputBox from "../../../components/UI/InputBox/InputBox";
import { motion } from "framer-motion";
import { fadeInUpVariantStatic } from "../../../components/animations/framerAnim";
import GlassyContainer from "../../../components/glassy_container/gc";
import InterButtons from "../../../components/buttons/InterButtons";
import { getAdmetPrediction } from "../../../services/api/mcpToolsService";
import DataViewer from "../../../components/UI/DataViewer/DataViewer";

const ADMETComponent = ({ toolData = null, initialSmiles = "" }) => {
    const [activeMol, setActiveMol] = useState(initialSmiles);
    const [isValidMol, setIsValidMol] = useState(false);
    const [predictionResults, setPredictionResults] = useState(toolData);
    const [isLoading, setIsLoading] = useState(false);

    // Handle initial prediction results from props
    useEffect(() => {
        if (toolData) {
            setPredictionResults(toolData);
        }
    }, [toolData]);

    // Handle initial SMILES from props
    useEffect(() => {
        if (initialSmiles) {
            setActiveMol(initialSmiles);
            // You might want to validate the SMILES here
            // For now, assuming it's valid if provided
            setIsValidMol(true);
        }
    }, [initialSmiles]);

    const handlePredict = async () => {
        if (!activeMol || !isValidMol) return;

        setIsLoading(true);
        try {
            // Use the MCP ADMET prediction service
            console.log('Calling ADMET prediction with SMILES:', activeMol);
            const response = await getAdmetPrediction({ smiles: activeMol });
            const result = response?.result
            if (!result) {
                console.error('No prediction result returned');
                setIsLoading(false);
                return;
            }

            let actualResults = null;
            if (typeof result === "object") {
                // Assuming result is already an object
                const result_data_json = JSON.parse(result["0"]);
                if (result_data_json.error === null) {
                    actualResults = result_data_json.data;
                } 
            }

            // Set the actual prediction results or fallback to demo data
            setPredictionResults(actualResults || { Status: "Prediction results not available" });
            setIsLoading(false);
        } catch (error) {
            console.error('Error predicting ADMET properties:', error);
            setIsLoading(false);
        }
    };

    const handleReset = () => {
        // If there were initial prediction results, reset to those instead of null
        setPredictionResults(null);
        setActiveMol("");
        setIsValidMol(false);
    };

    const labels = {
        "enabled": "Predict ADMET Properties",
        "disabled": "Enter Valid Molecule to Predict",
        "cancel": toolData ? "Reset to Initial" : "Reset"
    };

    const actors = {
        "enabled": handlePredict,
        "disabled": null,
        "cancel": handleReset
    };

    const getButtonStatus = () => {
        if (predictionResults) return "cancel";
        if (activeMol && isValidMol && !isLoading) return "enabled";
        return "disabled";
    };

    const renderPredictionResults = () => {
        if (!predictionResults) return null;

        return (
            <>
                <DataViewer 
                    data={predictionResults} 
                    title="Complete ADMET Prediction Data"
                    initiallyExpanded={false}
                />
            </>
        );
    };

    return (
        <motion.div
            initial="hidden"
            animate="visible"
            variants={fadeInUpVariantStatic}
            className="admet-container"
        >
            <div className="admet-row-1">
                <MolInputBox
                    activeMol={activeMol}
                    setActiveMol={setActiveMol}
                    isValidMol={isValidMol}
                    setIsValidMol={setIsValidMol}
                    inputType={"MOL"}
                    header="Enter Molecular SMILES for ADMET Prediction"
                    placeholder="Try CC(C)CC1=CC=C(C=C1)C(C)C(=O)O for demo"
                />
            </div>

            <div className="admet-row-2">
                <GlassyContainer>
                    <div className="predict-section">
                        <p className="predict-description">
                            Predict Absorption, Distribution, Metabolism, Excretion, and Toxicity (ADMET) properties
                            of your molecule using machine learning models.
                        </p>
                        {isLoading && (
                            <div className="loading-indicator">
                                <div className="loading-spinner"></div>
                                <p>Predicting ADMET properties...</p>
                            </div>
                        )}
                        <InterButtons
                            labels={labels}
                            actors={actors}
                            status={getButtonStatus()}
                        />
                    </div>
                </GlassyContainer>
            </div>

            {predictionResults && (
                <div className="admet-row-3">
                    {renderPredictionResults()}
                </div>
            )}
        </motion.div>
    );
};

export default ADMETComponent;