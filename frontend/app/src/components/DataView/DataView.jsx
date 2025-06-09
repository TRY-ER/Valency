import React from 'react';
import PropTypes from 'prop-types';

/**
 * Placeholder DataView component.
 * Replace this with your actual DataView component.
 * It's expected to receive data and display it.
 */
const DataView = ({ data }) => {
  if (!data) {
    return <p>No data to display.</p>;
  }

  // Simple JSON string representation for now
  return (
    <div style={{ marginTop: '20px', padding: '10px', border: '1px solid #ccc', borderRadius: '4px', backgroundColor: '#f9f9f9' }}>
      <h4>Data View</h4>
      <pre>{JSON.stringify(data, null, 2)}</pre>
    </div>
  );
};

DataView.propTypes = {
  data: PropTypes.any,
};

export default DataView;
