import React, { useEffect, useState } from "react";
import "./ToolInterface.css";
import menuContent from "../../contents/menuContent";

const  constructMasterTool = (masterMenu) =>{
    let masterTool = [];
    for (let i = 0; i < masterMenu.length; i++) {
      const items = {}
      items["id"] = i;
      items["title"] = masterMenu[i].title;
      items["link"] = masterMenu[i].link;
      items["description"] = masterMenu[i]?.description;
      items["subElements"] = [];
      if (masterMenu[i].subElements) {
        for (let j = 0; j < masterMenu[i].subElements.length; j++) {
          const subItems = {}
          subItems["id"] = j;
          subItems["title"] = masterMenu[i].subElements[j].title;
          subItems["link"] = masterMenu[i].subElements[j].link;
          subItems["description"] = masterMenu[i].subElements[j].description;
          subItems["subElements"] = [];
          if (masterMenu[i].subElements[j].subElements) {
            for (let k = 0; k < masterMenu[i].subElements[j].subElements.length; k++) {
              const subSubItems = {}
              subSubItems["id"] = k;
              subSubItems["title"] = masterMenu[i].subElements[j].subElements[k].title;
              subSubItems["link"] = masterMenu[i].subElements[j].subElements[k].link;
              subSubItems["description"] = masterMenu[i].subElements[j].subElements[k]?.description;
              subItems["subElements"].push(subSubItems);
            }
          }
          items["subElements"].push(subItems);
        }
      }
      masterTool.push(items);
    }
    return masterTool;
}

const ToolInterface = () => {
  const [masterTool, setMasterTool] = useState(constructMasterTool(menuContent));

  useEffect(() =>{
    console.log("master tool >>", masterTool);
  }, [masterTool])

  return (
    <div>
      <h1>Tool Interface</h1>
    </div>
  );
}

export default ToolInterface;