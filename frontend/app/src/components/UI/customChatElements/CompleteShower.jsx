import React from 'react';
import { BsCheck2All } from 'react-icons/bs'; // Double tick icon

const CompleteShower = ({ timestamp }) => {
    const options = {
        day: '2-digit',
        month: '2-digit',
        year: 'numeric',
        hour: '2-digit',
        minute: '2-digit',
        timeZoneName: 'short'
    };
    const formattedDate = timestamp ? new Date(timestamp).toLocaleString(undefined, options) : 'N/A';

    return (
        <div style={{ display: 'flex', alignItems: 'center', color: 'green', marginTop: '10px', marginBottom: '10px' }}>
            <BsCheck2All style={{ marginRight: '8px', fontSize: '1.2em' }} />
            <span>Stream completed at: {formattedDate}</span>
        </div>
    );
};

export default CompleteShower;
