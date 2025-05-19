import React, { useEffect, useRef } from "react";
import * as $3Dmol from "3dmol";

const ProtViewer = ({
    activeMol,
    setIsValidMol,
}) => {
    const viewerRef = useRef(null);
    const currentViewer = useRef(null); // To hold the viewer instance across renders for cleanup

    useEffect(() => {
        if (!viewerRef.current) {
            return;
        }

        // Clear previous viewer content and instance
        if (currentViewer.current) {
            currentViewer.current.clear(); // Or specific 3Dmol.js cleanup if available
        }
        viewerRef.current.innerHTML = ""; 

        if (!activeMol || activeMol.trim() === "") {
            setIsValidMol(false);
            return;
        }

        const viewerElement = viewerRef.current;
        const viewer = $3Dmol.createViewer(viewerElement, { backgroundColor: 'white' });
        currentViewer.current = viewer; // Store the new viewer instance

        let isValidInput = false;
        let isUrl = false;
        let urlToFetch = "";
        let pdbIdToLoad = "";

        if (activeMol.startsWith("http://") || activeMol.startsWith("https://") || activeMol.startsWith("ftp://")) {
            isValidInput = true;
            isUrl = true;
            urlToFetch = activeMol;
        } else {
            const query = activeMol.toUpperCase();
            if (query.match(/^[1-9][A-Za-z0-9]{3}$/)) {
                isValidInput = true;
                pdbIdToLoad = query;
            } else {
                isValidInput = false;
            }
        }

        setIsValidMol(isValidInput);

        if (isValidInput) {
            if (isUrl) {
                console.log("ProtViewer: Attempting to load model from URL:", urlToFetch);
                fetch(urlToFetch)
                    .then(response => {
                        if (!response.ok) {
                            throw new Error(`HTTP error! status: ${response.status}`);
                        }
                        return response.text();
                    })
                    .then(data => {
                        let fileType = urlToFetch.split('.').pop().toLowerCase();
                        if (fileType !== 'pdb' && fileType !== 'cif' && fileType !== 'mmcif') {
                            console.warn(`ProtViewer: Assuming PDB format for URL with extension: ${fileType}`);
                            fileType = 'pdb'; // Default to PDB if extension is not recognized
                        }
                        viewer.addModel(data, fileType);
                        viewer.setStyle({}, { cartoon: { color: 'spectrum' } });
                        viewer.zoomTo();
                        viewer.render();
                        console.log("ProtViewer: Model from URL rendered:", urlToFetch);
                    })
                    .catch(error => {
                        console.error("ProtViewer: Failed to load PDB from URL", urlToFetch, ":", error);
                        setIsValidMol(false);
                    });
            } else { // It's a PDB ID
                console.log("ProtViewer: Attempting to load model for PDB ID:", pdbIdToLoad);
                $3Dmol.download(`pdb:${pdbIdToLoad}`, viewer, { multimodel: true, frames: true }, function (loadedModel) {
                    if (loadedModel) {
                        console.log("ProtViewer: Model data for PDB ID successfully downloaded and processed.");
                        viewer.setStyle({}, { cartoon: { color: "spectrum" } });
                        viewer.zoomTo();
                        viewer.render();
                        console.log("ProtViewer: Model for PDB ID rendered:", pdbIdToLoad);
                    } else {
                        console.error("ProtViewer: Failed to load 3D model for PDB ID:", pdbIdToLoad);
                        setIsValidMol(false);
                    }
                });
            }
        } else {
             if (currentViewer.current) { // Ensure viewer is cleared if input becomes invalid
                currentViewer.current.clear();
            }
            viewerRef.current.innerHTML = ""; 
        }

        return () => {
            // Cleanup when the component unmounts or before re-running the effect
            if (currentViewer.current) {
                // 3Dmol.js viewers don't have a dedicated unmount/destroy. Clearing should suffice.
                // Or, if specific cleanup methods are found in 3Dmol.js docs, use them.
                currentViewer.current.clear(); 
                currentViewer.current = null;
            }
            if (viewerRef.current) {
                 viewerRef.current.innerHTML = ""; // Ensure DOM is clean
            }
        };

    }, [activeMol, setIsValidMol]);

    return (
        <div>
            <div
                id="protviewer-unique-id" // Changed ID to be more specific
                ref={viewerRef}
                style={{ width: '100%', height: '400px', position: 'relative' }}
            ></div>
        </div>
    );
};

export default ProtViewer;