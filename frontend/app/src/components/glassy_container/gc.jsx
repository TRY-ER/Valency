import React, { useState } from 'react';
import './gc.css';
import { FaExpand } from "react-icons/fa6";
import { IoClose } from "react-icons/io5";

export default function GlassyContainer({
  children,
  mode = 'enabled',
  expandable = false, // Add this prop to enable or disable the expand feature
  className = ""
}) {
  const [isModalOpen, setIsModalOpen] = useState(false);

  const isLoading = mode === 'loading';
  const isDisabled = mode === 'disabled';
  const containerClassName = `glassy-container ${isDisabled ? 'disabled' : ''} ${isLoading ? 'loading' : ''
    }`.trim();

  const handleExpandClick = () => {
    setIsModalOpen(true);
  };

  const handleCloseModal = () => {
    setIsModalOpen(false);
  };

  return (
    <>
      <div className={`${containerClassName} ${className}`}>
        {expandable && (
          <FaExpand className="expand-button" onClick={handleExpandClick} />
        )}
        {isLoading ? null : children}
      </div>

      {isModalOpen && (
        <div className="modal-overlay">
          <div className="modal-content">
              <IoClose className="close-button" onClick={handleCloseModal} />
            {children}
          </div>
        </div>
      )}
    </>
  );
}