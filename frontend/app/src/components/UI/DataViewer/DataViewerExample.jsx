import React from 'react';
import DataViewer from './DataViewer';

// Example usage of the DataViewer component
const DataViewerExample = () => {
    // Sample complex data structures
    const sampleObjectData = {
        user: {
            id: 12345,
            name: "John Doe",
            email: "john.doe@example.com",
            preferences: {
                theme: "dark",
                notifications: {
                    email: true,
                    push: false,
                    sms: null
                },
                languages: ["en", "es", "fr"]
            },
            metadata: {
                lastLogin: "2025-06-04T10:30:00Z",
                accountType: "premium",
                isActive: true
            }
        },
        settings: {
            privacy: {
                profileVisible: true,
                dataSharing: false
            },
            features: ["analytics", "reporting", "api_access"]
        }
    };

    const sampleArrayData = [
        {
            id: 1,
            type: "absorption",
            value: 0.85,
            confidence: 0.92,
            details: {
                method: "ML",
                model_version: "v1.2.3",
                parameters: {
                    temperature: 37,
                    pH: 7.4
                }
            }
        },
        {
            id: 2,
            type: "distribution",
            value: 0.72,
            confidence: 0.88,
            details: {
                method: "DL",
                model_version: "v2.1.0",
                parameters: {
                    tissue_type: "plasma",
                    binding_affinity: 0.65
                }
            }
        },
        {
            id: 3,
            type: "metabolism",
            value: 0.68,
            confidence: 0.79,
            details: null
        }
    ];

    const admetPredictionData = {
        molecule: {
            smiles: "CC(C)CC1=CC=C(C=C1)C(C)C(=O)O",
            molecular_weight: 206.28,
            formula: "C13H18O2"
        },
        predictions: {
            absorption: {
                value: 0.85,
                confidence: 0.92,
                model_info: {
                    name: "ADMET_ABS_v1.0",
                    version: "1.0.0",
                    training_data_size: 50000
                }
            },
            distribution: {
                value: 0.72,
                confidence: 0.88,
                model_info: {
                    name: "ADMET_DIST_v1.0",
                    version: "1.0.0",
                    training_data_size: 45000
                }
            },
            metabolism: {
                value: 0.68,
                confidence: 0.79,
                model_info: {
                    name: "ADMET_MET_v1.0",
                    version: "1.0.0",
                    training_data_size: 42000
                }
            },
            excretion: {
                value: 0.91,
                confidence: 0.94,
                model_info: {
                    name: "ADMET_EXC_v1.0",
                    version: "1.0.0",
                    training_data_size: 38000
                }
            },
            toxicity: {
                value: 0.23,
                confidence: 0.87,
                model_info: {
                    name: "ADMET_TOX_v1.0",
                    version: "1.0.0",
                    training_data_size: 55000
                }
            }
        },
        metadata: {
            prediction_time: "2025-06-04T10:30:00Z",
            processing_time_ms: 1250,
            api_version: "2.1.0",
            request_id: "req_abc123def456"
        }
    };

    return (
        <div style={{ padding: '20px', maxWidth: '1200px', margin: '0 auto' }}>
            <h1>DataViewer Component Examples</h1>
            
            <DataViewer 
                data={sampleObjectData}
                title="User Profile Data"
                initiallyExpanded={true}
            />
            
            <DataViewer 
                data={sampleArrayData}
                title="Prediction Results Array"
                initiallyExpanded={false}
            />
            
            <DataViewer 
                data={admetPredictionData}
                title="ADMET Prediction Results"
                initiallyExpanded={true}
                maxDepth={8}
            />
            
            <DataViewer 
                data={null}
                title="Empty Data Example"
            />
            
            <DataViewer 
                data="Simple string value"
                title="Simple Value Example"
            />
            
            <DataViewer 
                data={[1, 2, 3, "hello", true, null, { nested: "object" }]}
                title="Mixed Array Example"
            />
        </div>
    );
};

export default DataViewerExample;
