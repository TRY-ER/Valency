import React, {useEffect, useState} from 'react';
import './Discriminator.css';
import GlassyContainer from '../../components/glassy_container/gc';
import TabContainer from '../../components/tab_section/TabSection';
import DiscriminatorTabContents from '../../contents/tag_content/DiscriminatorTags';
import FunctionalSection from '../../components/functional_section/Functional';
import Divider from '../../components/divider';
import MoEDocs from '../../contents/doc_content/explorer_content/MoEDocs';

export default function Discriminator() {
    const [activeTabId, setActiveTabId] = useState(1);

    return (
        <div className="base-page-container">
            <TabContainer 
                tabDetails={DiscriminatorTabContents}
                activeTabId={activeTabId}
                setActiveTabId={setActiveTabId}
            />
            <Divider />
            <FunctionalSection 
                docElem={DiscriminatorTabContents[activeTabId - 1].docs} 
                funcElem={DiscriminatorTabContents[activeTabId - 1].component}
            />
            {/* {functionalContent}  */}
        </div>
    )
}