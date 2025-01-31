// setup react component

import React, {useEffect, useState} from 'react';
import './Explorer.css';
import GlassyContainer from '../../components/glassy_container/gc';
import TabContainer from '../../components/tab_section/TabSection';
import ExploreTabContent from '../../contents/tag_content/ExploreTags';
import FunctionalSection from '../../components/functional_section/Functional';
import Divider from '../../components/divider';
import MoEDocs from '../../contents/doc_content/explorer_content/MoEDocs';

export default function Explorer() {
    const [activeTabId, setActiveTabId] = useState(1);

    return (
        <div className="base-page-container">
            <TabContainer 
                tabDetails={ExploreTabContent}
                activeTabId={activeTabId}
                setActiveTabId={setActiveTabId}
            />
            <Divider />
            <FunctionalSection 
                docElem={ExploreTabContent[activeTabId - 1].docs} 
                funcElem={ExploreTabContent[activeTabId - 1].component}
            />
            {/* {functionalContent}  */}
        </div>
    )
}