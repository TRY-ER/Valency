import MoEDocs from "../doc_content/explorer_content/MoEDocs";
import ProDocs from "../doc_content/explorer_content/ProDocs";
import PolyDocs from "../doc_content/explorer_content/PolyDocs";

import MolEComponent from "../../pages/Explorer/MolExplorer/MolExplorer";
import ProtEComponent from "../../pages/Explorer/ProtExplorer/ProtExplorer";
import PolyEComponent from "../../pages/Explorer/PolyExplorer/PolyExplorer";


const ExplorerTabContent = [
    {
        id: 1,
        title: 'Molecule Explorer',
        link: '',
        description: 'Explore the Molecule Explorer',
        component: <MolEComponent />,
        docs: <MoEDocs />,
    },
    {
        id: 2,
        title: 'Protein Explorer',
        link: 'proe',
        description: 'Explore the Protein Explorer',
        component: <ProtEComponent />,
        docs: <ProDocs />,
    },
    {
        id: 3,
        title: 'Polymer Explorer',
        link: 'polye',
        description: 'Explore the Polymer Explorer',
        component: <PolyEComponent />,
        docs: <PolyDocs />,
    }
]

export default ExplorerTabContent;