import React from 'react';
import './SimilarityCardContainer.css';
import GlassyContainer from '../../glassy_container/gc';

export default function SimilarityCardContainer({ data, status }) {
    if (!data || !Array.isArray(data)) return null;
    return (
        <div className="similarity-card-container">
            { status === "loading" && <GlassyContainer mode="loading" />}
            {data.map((item, index) => (
                <GlassyContainer>
                    <div className="similarity-card" key={index}>
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