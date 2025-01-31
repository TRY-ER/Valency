import React, { useEffect } from 'react';

export default function InterButtons({ labels, actors, status = 'enabled'}) {

  const enabledStyle = {
    backgroundColor: 'var(--color-success, #28a745)',
    color: '#fff',
    border: 'none',
    borderRadius: '4px',
    padding: '10px 16px',
    cursor: 'pointer',
    width: 'auto'
  };

  const disabledStyle = {
    backgroundColor: '#6c757d',
    color: '#fff',
    border: 'none',
    borderRadius: '4px',
    padding: '10px 16px',
    cursor: 'not-allowed',
  };

  const cancelStyle = {
    backgroundColor: 'var(--color-alert, #dc3545)',
    color: '#fff',
    border: 'none',
    borderRadius: '4px',
    padding: '10px 16px',
    cursor: 'pointer',
  };

  let buttonStyle = enabledStyle;
  let label = labels["enabled"];

  if (status === 'disabled') {
    buttonStyle = disabledStyle;
    label = labels["disabled"]; 
  } else if (status === 'cancel') {
    buttonStyle = cancelStyle;
    label = labels["cancel"];
  }

  return (
    <div style={{display: 'flex', justifyContent: 'flex-start'}}>
    <button style={buttonStyle} onClick={actors[status]} disabled={status === 'disabled'}>
      {label}
    </button>
    </div>
  );
}