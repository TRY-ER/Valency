import React, {useEffect, useState} from 'react';
import './Generator.css';
import GlassyContainer from '../../components/glassy_container/gc';
import TabContainer from '../../components/tab_section/TabSection';
import FunctionalSection from '../../components/functional_section/Functional';
import GeneratorTabContents from '../../contents/tag_content/GeneratorTags';
import Divider from '../../components/divider';
import MoEDocs from '../../contents/doc_content/explorer_content/MoEDocs';

export default function Generator() {
    const [activeTabId, setActiveTabId] = useState(1);

    return (
        <div className="base-page-container">
            <TabContainer 
                tabDetails={GeneratorTabContents}
                activeTabId={activeTabId}
                setActiveTabId={setActiveTabId}
            />
            <Divider />
            <FunctionalSection 
                docElem={GeneratorTabContents[activeTabId - 1].docs} 
                funcElem={GeneratorTabContents[activeTabId - 1].component}
            />
        </div>
    )
}