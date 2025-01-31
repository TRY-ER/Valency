import React, { useEffect, useState } from 'react';
import TabContainer from '../../../components/tab_section/TabSection';
import FunctionalSection from '../../../components/functional_section/Functional';
import Divider from '../../../components/divider';
import { LSTMGeneratorTabContents } from '../../../contents/tag_content/NestGeneratorTags';
import "./LSTMGenerator.css";

export default function LSTMGenerator() {
    const [activeTabId, setActiveTabId] = useState(1);

    return (
        <>
            <TabContainer
                tabDetails={LSTMGeneratorTabContents}
                activeTabId={activeTabId}
                setActiveTabId={setActiveTabId}
            />
            <Divider />
            <FunctionalSection
                docElem={LSTMGeneratorTabContents[activeTabId - 1].docs}
                funcElem={LSTMGeneratorTabContents[activeTabId - 1].component}
                customClassName={"functional-container-nested"}
            />
            {/* {functionalContent}  */}
        </>
    )
}