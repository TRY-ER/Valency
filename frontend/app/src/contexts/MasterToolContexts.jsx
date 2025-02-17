//// filepath: /home/kalki/imported/career/Competitions/2025/Solutions Challenge GDG/source_code/frontend/app/src/context/MasterToolContext.jsx
import React, { createContext, useState, useEffect } from "react";
import menuContent from "../contents/menuContent";
import { call_endpoint_async } from "../endpoints/caller";
import { chat_endpoint } from "../endpoints/endpoints";

// Copy your constructMasterTool function here or import it if defined elsewhere.
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
        subItems["title"] = masterMenu[i].title + " > " + masterMenu[i].subElements[j].title;
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
            subSubItems["title"] = masterMenu[i].title + " > " + masterMenu[i].subElements[j].title + " > " + masterMenu[i].subElements[j].subElements[k].title;
            if (masterMenu[i].subElements[j].subElements[k].link !== "") {
              subSubItems["link"] = `${pre}/${masterMenu[i].subElements[j].subElements[k].link}`;
            } else {
              subSubItems["link"] = `${pre}`;
            }
            subSubItems["description"] = masterMenu[i].subElements[j].subElements[k]?.description;
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

export const MasterToolContext = createContext();

export const MasterToolProvider = ({ children }) => {
  const backupMasterTool = JSON.parse(localStorage.getItem("masterTool"));
  const initialState = backupMasterTool ? backupMasterTool : constructMasterTool(menuContent);

  const [masterTool, setMasterTool] = useState(initialState);

  // Keep localStorage updated
  useEffect(() => {
    localStorage.setItem("masterTool", JSON.stringify(masterTool));
  }, [masterTool]);

  // Update tool config once after mount
  useEffect(() => {
    const setToolConfig = async () => {
      if (masterTool !== "" && masterTool.length) {
        await call_endpoint_async(chat_endpoint.config, { data: masterTool });
      }
    };
    setToolConfig();
  }, []);

  return (
    <MasterToolContext.Provider value={{ masterTool, setMasterTool }}>
      {children}
    </MasterToolContext.Provider>
  );
};