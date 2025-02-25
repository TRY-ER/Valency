import React, { useContext } from "react";
import "./ToolInterface.css";
import GlassyContainer from "../glassy_container/gc";
import { MdOutlineKeyboardArrowDown, MdOutlineKeyboardArrowUp } from "react-icons/md";
import { MasterToolContext } from "../../contexts/MasterToolContexts";
import { motion } from "framer-motion";
import { fadeInRightVariants } from "../animations/framerAnim";

const ToolInterface = () => {
  const { masterTool, setMasterTool } = useContext(MasterToolContext);

  // Toggle expansion for top-level tool items
  const toggleToolExpansion = (toolId) => {
    setMasterTool((prevTools) =>
      prevTools.map((tool) =>
        tool.id === toolId ? { ...tool, isExpanded: !tool.isExpanded } : tool
      )
    );
  };

  // Toggle expansion for sub tool items
  const toggleSubToolExpansion = (toolId, subToolId) => {
    setMasterTool((prevTools) =>
      prevTools.map((tool) => {
        if (tool.id === toolId) {
          const updatedSubElements = tool.subElements.map((subTool) =>
            subTool.id === subToolId ? { ...subTool, isExpanded: !subTool.isExpanded } : subTool
          );
          return { ...tool, subElements: updatedSubElements };
        }
        return tool;
      })
    );
  };

  // Handle description change for top-level tool
  const handleToolDescriptionChange = (toolId, value) => {
    setMasterTool((prevTools) =>
      prevTools.map((tool) =>
        tool.id === toolId ? { ...tool, description: value } : tool
      )
    );
  };

  // Handle description change for sub tool
  const handleSubToolDescriptionChange = (toolId, subToolId, value) => {
    setMasterTool((prevTools) =>
      prevTools.map((tool) => {
        if (tool.id === toolId) {
          const updatedSubElements = tool.subElements.map((subTool) =>
            subTool.id === subToolId ? { ...subTool, description: value } : subTool
          );
          return { ...tool, subElements: updatedSubElements };
        }
        return tool;
      })
    );
  };

  // Handle description change for sub-sub tool
  const handleSubSubToolDescriptionChange = (toolId, subToolId, subSubToolId, value) => {
    setMasterTool((prevTools) =>
      prevTools.map((tool) => {
        if (tool.id === toolId) {
          const updatedSubElements = tool.subElements.map((subTool) => {
            if (subTool.id === subToolId) {
              const updatedSubSubElements = subTool.subElements.map((subSubTool) =>
                subSubTool.id === subSubToolId
                  ? { ...subSubTool, description: value }
                  : subSubTool
              );
              return { ...subTool, subElements: updatedSubSubElements };
            }
            return subTool;
          });
          return { ...tool, subElements: updatedSubElements };
        }
        return tool;
      })
    );
  };

  const formTool = (tool, index) => {
    return (
      <motion.div
        variants={fadeInRightVariants}
        initial="hidden"
        animate="visible"
        custom={index}
      >
        <React.Fragment key={tool.id}>
          <GlassyContainer>
            <div className="tool-container">
              <div className="tool-header">
                <p className="tool-title">
                  Name: <span className="tool-tag">{tool.title}</span>
                </p>
                <p className="tool-about">
                  Link: <span className="tool-tag">`{tool.link}`</span>
                </p>
              </div>
              <div className="tool-desc">
                <p className="tool-desc-text">Description:</p>
                <textarea
                  className="tool-desc-textarea"
                  value={tool.description}
                  onChange={(e) => handleToolDescriptionChange(tool.id, e.target.value)}
                ></textarea>
              </div>
              {tool.subElements && tool.subElements.length > 0 && (
                <>
                  {tool.isExpanded ? (
                    <MdOutlineKeyboardArrowUp
                      className="custom-down-arrow"
                      onClick={() => toggleToolExpansion(tool.id)}
                    />
                  ) : (
                    <MdOutlineKeyboardArrowDown
                      className="custom-down-arrow"
                      onClick={() => toggleToolExpansion(tool.id)}
                    />
                  )}
                </>
              )}
            </div>
          </GlassyContainer>
          {tool.isExpanded && (
            <div className="tool-sub-elements">
              {tool.subElements.map((subTool, index) => (
                <React.Fragment key={subTool.id}>
                  <motion.div
                    variants={fadeInRightVariants}
                    initial="hidden"
                    animate="visible"
                    custom={index}
                  >
                    <GlassyContainer>
                      <div className="tool-container">
                        <div className="tool-header">
                          <p className="tool-title">
                            Name: <span className="tool-tag">{subTool.title}</span>
                          </p>
                          <p className="tool-about">
                            Link: <span className="tool-tag">`{subTool.link}`</span>
                          </p>
                        </div>
                        <div className="tool-desc">
                          <p className="tool-desc-text">Description:</p>
                          <textarea
                            className="tool-desc-textarea"
                            value={subTool.description}
                            onChange={(e) =>
                              handleSubToolDescriptionChange(tool.id, subTool.id, e.target.value)
                            }
                          ></textarea>
                        </div>
                        {subTool.subElements && subTool.subElements.length > 0 && (
                          <>
                            {subTool.isExpanded ? (
                              <MdOutlineKeyboardArrowUp
                                className="custom-down-arrow"
                                onClick={() => toggleSubToolExpansion(tool.id, subTool.id)}
                              />
                            ) : (
                              <MdOutlineKeyboardArrowDown
                                className="custom-down-arrow"
                                onClick={() => toggleSubToolExpansion(tool.id, subTool.id)}
                              />
                            )}
                          </>
                        )}
                      </div>
                    </GlassyContainer>
                  </motion.div>
                  {subTool.isExpanded && (
                    <div className="tool-sub-sub-elements">
                      {subTool.subElements.map((subSubTool, index) => (
                        <motion.div
                          variants={fadeInRightVariants}
                          initial="hidden"
                          animate="visible"
                          custom={index} 
                        >
                          <GlassyContainer key={subSubTool.id}>
                            <div className="tool-container">
                              <div className="tool-header">
                                <p className="tool-title">
                                  Name: <span className="tool-tag">{subSubTool.title}</span>
                                </p>
                                <p className="tool-about">
                                  Link: <span className="tool-tag">`{subSubTool.link}`</span>
                                </p>
                              </div>
                              <div className="tool-desc">
                                <p className="tool-desc-text">Description:</p>
                                <textarea
                                  className="tool-desc-textarea"
                                  value={subSubTool.description}
                                  onChange={(e) =>
                                    handleSubSubToolDescriptionChange(
                                      tool.id,
                                      subTool.id,
                                      subSubTool.id,
                                      e.target.value
                                    )
                                  }
                                ></textarea>
                              </div>
                            </div>
                          </GlassyContainer>
                        </motion.div>
                      ))}
                    </div>
                  )}
                </React.Fragment>
              ))}
            </div>
          )}
        </React.Fragment>
      </motion.div>
    );
  };

  return (
    <div className="tools-master-container">
      <div className="tools-section">
        {masterTool.map((tool, index) => formTool(tool, index))}
      </div>
    </div>
  );
};

export default ToolInterface;
