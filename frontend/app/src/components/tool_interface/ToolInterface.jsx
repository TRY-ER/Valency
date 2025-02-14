import React, { useEffect, useState } from "react";
import "./ToolInterface.css";
import menuContent from "../../contents/menuContent";
import GlassyContainer from "../glassy_container/gc";
import { MdOutlineKeyboardArrowDown, MdOutlineKeyboardArrowUp } from "react-icons/md";
import { chat_endpoint } from "../../endpoints/endpoints";
import { call_endpoint_async } from "../../endpoints/caller";

const ExcludeToolTitles = ["Chatbot", "Tool Config"];

const constructMasterTool = (masterMenu) => {
  let masterTool = [];
  for (let i = 0; i < masterMenu.length; i++) {
    if (ExcludeToolTitles.includes(masterMenu[i].title)) {
      continue;
    }
    const items = {};
    items["id"] = i;
    items["title"] = masterMenu[i].title;
    items["link"] = masterMenu[i].link;
    items["description"] = masterMenu[i]?.description;
    items["isExpanded"] = false;
    items["subElements"] = [];
    if (masterMenu[i].subElements) {
      for (let j = 0; j < masterMenu[i].subElements.length; j++) {
        const subItems = {};
        subItems["id"] = j;
        subItems["title"] = masterMenu[i].subElements[j].title;
        let pre = "";
        if (masterMenu[i].subElements[j].link !== "") {
          subItems["link"] = `${masterMenu[i].link}/${masterMenu[i].subElements[j].link}`;
          pre = `${masterMenu[i].link}/${masterMenu[i].subElements[j].link}`;
        } else {
          subItems["link"] = `${masterMenu[i].link}`;
          pre = `${masterMenu[i].link}`;
        }
        subItems["description"] = masterMenu[i].subElements[j].description;
        subItems["isExpanded"] = false;
        subItems["subElements"] = [];
        if (masterMenu[i].subElements[j].subElements) {
          for (let k = 0; k < masterMenu[i].subElements[j].subElements.length; k++) {
            const subSubItems = {};
            subSubItems["id"] = k;
            subSubItems["title"] = masterMenu[i].subElements[j].subElements[k].title;
            if (masterMenu[i].subElements[j].subElements[k].link !== "") {
              subSubItems["link"] = `${pre}/${masterMenu[i].subElements[j].subElements[k].link}`;
            } else {
              subSubItems["link"] = `${pre}`;
            }
            subSubItems["description"] = masterMenu[i].subElements[j].subElements[k]?.description;
            // You may add an expansion flag for further nested levels if needed.
            subSubItems["isExpanded"] = false;
            subItems["subElements"].push(subSubItems);
          }
        }
        items["subElements"].push(subItems);
      }
    }
    masterTool.push(items);
  }
  return masterTool;
};

const ToolInterface = () => {
  const backupMasterTool = JSON.parse(localStorage.getItem("masterTool"));
  var processible = ""
  if (backupMasterTool) {
    processible = backupMasterTool
  }
  else {
    processible = constructMasterTool(menuContent);
  }

  const [masterTool, setMasterTool] = useState(processible);

  useEffect(() => {
    console.log("master tool >>", masterTool);
    localStorage.setItem("masterTool", JSON.stringify(masterTool));
  }, [masterTool]);

  useEffect(() => {
    const setToolConfig = async () => {
      if (masterTool !== "") {
        console.log(" >>", masterTool);
        const response = await call_endpoint_async(chat_endpoint.config, {"data" : masterTool}) 
        console.log("tool config response >>", response);
      }
    }
    setToolConfig();
  }, [])


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

  const formTool = (tool) => {
    return (
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
            {tool.subElements.map((subTool) => (
              <React.Fragment key={subTool.id}>
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
                {subTool.isExpanded && (
                  <div className="tool-sub-sub-elements">
                    {subTool.subElements.map((subSubTool) => (
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
                    ))}
                  </div>
                )}
              </React.Fragment>
            ))}
          </div>
        )}
      </React.Fragment>
    );
  };

  return (
    <div className="tools-master-container">
      <div className="tools-section">
        {masterTool.map((tool) => formTool(tool))}
      </div>
    </div>
  );
};

export default ToolInterface;
