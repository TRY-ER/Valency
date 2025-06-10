import React from 'react';
import './SimilarityCardContainer.css';
import GlassyContainer from '../../glassy_container/gc';

export default function SimilarityCardContainer({ data, status }) {
    if (!data || !Array.isArray(data)) return null;
    
    return (
        <div className="similarity-card-container">
            {status === "loading" && <GlassyContainer mode="loading" />}
            {status === "no-data" && (
                <GlassyContainer>
                    <div style={{ 
                        textAlign: 'center', 
                        padding: '2rem',
                        color: 'var(--color-text-secondary)' 
                    }}>
                        <h4 style={{ 
                            marginBottom: '0.5rem',
                            color: 'var(--color-text-primary)' 
                        }}>
                            No Matching Data Found
                        </h4>
                        <p style={{ margin: 0 }}>
                            No similar molecules were found for the given input. Try adjusting your search parameters.
                        </p>
                    </div>
                </GlassyContainer>
            )}
            {status === "completed" && data.map((item, index) => (
                <GlassyContainer key={index}>
                    <div className="similarity-card">
                        <img
                            src={`data:image/png;base64,${item.image}`}
                            alt={item.identifier}
                        />
                        <div className="card-details">
                            <p><strong>{item.identifier}</strong></p>
                            <p>Score: {item.score.toFixed(3)}</p>
                        </div>
                    </div>
                </GlassyContainer>
            ))}
        </div>
    );
}