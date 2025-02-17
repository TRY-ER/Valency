import React from "react";
import "./FileUploader.css"; // or your styling
import GlassyContainer from "../glassy_container/gc";
import { motion } from "framer-motion";
import { fadeInUpVariantStatic } from "../animations/framerAnim";

export default function FileDownloader({ status,
    generatedData,
    consoleContent,
    handleReset,
    expandable = false,
    showGenerated = true}) {
    // if (!generatedData || generatedData.length === 0) {
    //     return (
    //         <div className="file-download-section">
    //             <GlassyContainer >
    //             </GlassyContainer>
    //         </div>
    //     )
    // }

    const handleDownload = () => {
        const blob = new Blob([generatedData.join("\n")], { type: "text/plain" });
        const url = URL.createObjectURL(blob);
        const link = document.createElement("a");
        link.href = url;
        link.download = "generated_output.txt";
        document.body.appendChild(link);
        link.click();
        document.body.removeChild(link);
    };

    return (
        <motion.div 
            variants={fadeInUpVariantStatic}
            initial="hidden"
            animate="visible" 
        className="file-download-section">
            <GlassyContainer expandable={expandable}>
                <div className="file-download-inner-container">
                    <div className="console-container">
                        {consoleContent?.map((line, idx) => (
                            line
                        ))}
                        {
                            status === "completed" && generatedData.length === 0 &&
                            <p className="console-alert">No Candidates generated</p>
                        }
                        {
                            status === "completed" && generatedData.length > 0 && showGenerated &&
                            <>
                                <br />
                                <p>Generated Sample Candidates</p>
                                <p>-------------------------------------</p>
                                {generatedData?.slice(0, 5).map((line, idx) => (
                                    <p key={idx}>{idx + 1}. {line}</p>
                                ))}
                                {generatedData?.length > 5 && <p>...</p>}
                            </>
                        }


                    </div>
                    {
                        status === "completed" && generatedData.length > 0 &&
                        <button onClick={handleDownload} className="download-btn">Download</button>
                    }
                    {
                        status === "completed" || status === "failed" ?
                            <button className="download-btn" onClick={handleReset}>Reset</button>
                            : ""
                    }
                </div>
            </GlassyContainer>
        </motion.div>
    );
}