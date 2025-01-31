import MoEDocs from "../doc_content/explorer_content/MoEDocs";
import ProDocs from "../doc_content/explorer_content/ProDocs";
import PolyDocs from "../doc_content/explorer_content/PolyDocs";

import SSComponent from "../../pages/Discriminator/SimilaritySearch/SSearch";

const DiscriminatorTabContents = [
    {
        id: 1,
        title: 'Similarity Search',
        link: '/ssearch',
        component: <SSComponent />, 
        docs: null,
    },
];

export default DiscriminatorTabContents;