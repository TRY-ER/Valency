import React, { useEffect, useState } from 'react';
import TabContainer from '../../../components/tab_section/TabSection';
import FunctionalSection from '../../../components/functional_section/Functional';
import Divider from '../../../components/divider';
import BRICSGeneratorTabContents from '../../../contents/tag_content/NestGeneratorTags';
import "./Generator.css";

export default function BRICSGenerator() {
    const [activeTabId, setActiveTabId] = useState(1);

    return (
        <>
            <TabContainer
                tabDetails={BRICSGeneratorTabContents}
                activeTabId={activeTabId}
                setActiveTabId={setActiveTabId}
            />
            <Divider />
            <FunctionalSection
                docElem={BRICSGeneratorTabContents[activeTabId - 1].docs}
                funcElem={BRICSGeneratorTabContents[activeTabId - 1].component}
                customClassName={"functional-container-nested"}
            />
            {/* {functionalContent}  */}
        </>
    )
}