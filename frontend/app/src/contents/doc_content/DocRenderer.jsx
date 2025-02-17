//// filepath: /home/kalki/imported/career/Competitions/2025/Solutions Challenge GDG/source_code/frontend/app/src/contents/doc_content/explorer_content/MoEDocs.js
import React, { useEffect, useState } from "react";
import { motion } from "framer-motion";
import ReactMarkdown from "react-markdown";
import { fadeInRightVariantStatic } from "../../components/animations/framerAnim";
import "./DocRenderer.css";

const DocRenderer = ({ filePath }) => {
    const [docContent, setDocContent] = useState("");

    useEffect(() => {
        console.log('docContent:', docContent)
        console.log('filePath:', filePath)
    }, [docContent])

    useEffect(() => {
        fetch(filePath)
            .then((res) => res.text())
            .then((text) => {
                console.log("text >>", text)
                setDocContent(text);
            })
            .catch((err) => console.error("Error loading markdown:", err));
    }, [filePath]);

    return (
        <motion.div
            initial="hidden"
            animate="visible"
            variants={fadeInRightVariantStatic}
        >
            <div className="doc-renderer-wrapper">
                <ReactMarkdown>{docContent}</ReactMarkdown>
            </div>
        </motion.div>
    );
};

export default DocRenderer;