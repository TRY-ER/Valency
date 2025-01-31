import React, {useState, useEffect} from "react";
import TabContainer from "../../../components/tab_section/TabSection";
import Divider from "../../../components/divider";
import FunctionalSection from "../../../components/functional_section/Functional";
import { SimilaritySearchTabContents } from "../../../contents/tag_content/NestDiscriminatorTags";

const SSComponent = () => {
    const [activeTabId, setActiveTabId] = useState(1);

    return (
        <>
            <TabContainer
                tabDetails={SimilaritySearchTabContents}
                activeTabId={activeTabId}
                setActiveTabId={setActiveTabId}
            />
            <Divider />
            <FunctionalSection
                docElem={SimilaritySearchTabContents[activeTabId - 1].docs}
                funcElem={SimilaritySearchTabContents[activeTabId - 1].component}
                customClassName={"functional-container-nested"}
            />
        </>
    ) 
}

export default SSComponent;