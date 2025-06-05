import React from 'react';
import AlphafoldAnnotationsViewer from '../pages/Explorer/ProtExp// With specific annotation type
<AlphafoldAnnotationsViewer 
    initialAccessionKey="Q5VSL9"
    initialAnnotationType="MUTAGEN"
/>

// With pre-loaded data
<AlphafoldAnnotationsViewer 
    toolData={annotationsData}
    hideInputBox={true}
/>

// Handling data from API
const handleAnnotationsReceived = async () => {
    const response = await getAlphafoldAnnotations({ 
        qualifier: 'Q5VSL9',
        annotation_type: 'MUTAGEN'
    });ationsViewer';

/**
 * Demo component showing usage examples of the AlphaFold Annotations Viewer
 */
const AlphafoldAnnotationsDemo = () => {
    // Example data that might be returned from the API
    const sampleAnnotationsData = {
        annotations: [
            {
                type: "Domain",
                start: 10,
                end: 50,
                description: "DNA-binding domain",
                label: "DBD"
            },
            {
                type: "Active Site",
                start: 25,
                end: 25,
                description: "Catalytic residue",
                label: "Active"
            },
            {
                type: "Secondary Structure",
                start: 60,
                end: 80,
                description: "Alpha helix",
                label: "Helix"
            },
            {
                type: "Binding Site",
                start: 100,
                end: 110,
                description: "ATP binding site",
                label: "ATP"
            }
        ]
    };

    return (
        <div style={{ padding: '20px' }}>
            <h1>AlphaFold Annotations Viewer - Demo</h1>
            
            <h2>Basic Usage Examples</h2>
            
            <div style={{ marginBottom: '40px' }}>
                <h3>1. Default Viewer (with input box)</h3>
                <p>Users can enter a UniProt accession and fetch annotations:</p>
                <AlphafoldAnnotationsViewer />
            </div>

            <div style={{ marginBottom: '40px' }}>
                <h3>2. Pre-loaded with Initial Accession</h3>
                <p>Viewer with an initial UniProt accession (Q5VSL9):</p>
                <AlphafoldAnnotationsViewer 
                    initialAccessionKey="Q5VSL9"
                />
            </div>

            <div style={{ marginBottom: '40px' }}>
                <h3>3. Pre-loaded with Data (Hidden Input)</h3>
                <p>Viewer with sample data and hidden input box:</p>
                <AlphafoldAnnotationsViewer 
                    toolData={sampleAnnotationsData}
                    hideInputBox={true}
                />
            </div>

            <div style={{ marginBottom: '40px' }}>
                <h3>4. With Predefined Annotation Type</h3>
                <p>Viewer with predefined accession and annotation type:</p>
                <AlphafoldAnnotationsViewer 
                    initialAccessionKey="Q5VSL9"
                    initialAnnotationType="MUTAGEN"
                />
            </div>

            <div style={{ marginBottom: '40px' }}>
                <h3>Usage in Your Components</h3>
                <pre style={{ 
                    backgroundColor: '#f4f4f4', 
                    padding: '15px', 
                    borderRadius: '5px',
                    overflow: 'auto'
                }}>
{`import AlphafoldAnnotationsViewer from './AlphafoldAnnotationsViewer';

// Basic usage
<AlphafoldAnnotationsViewer />

// With initial data
<AlphafoldAnnotationsViewer 
    initialAccessionKey="Q5VSL9"
    startPosition={1}
    endPosition={100}
/>

// With pre-loaded data
<AlphafoldAnnotationsViewer 
    toolData={annotationsData}
    hideInputBox={true}
/>

// Handling data from API
const handleAnnotationsReceived = async () => {
    const response = await getAlphafoldAnnotations({ 
        uniprot_accession: 'Q5VSL9',
        start: 1,
        end: 100
    });
    
    return (
        <AlphafoldAnnotationsViewer 
            toolData={response.result}
            hideInputBox={true}
        />
    );
};`}
                </pre>
            </div>
        </div>
    );
};

export default AlphafoldAnnotationsDemo;
