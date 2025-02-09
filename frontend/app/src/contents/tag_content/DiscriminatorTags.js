import MoEDocs from "../doc_content/explorer_content/MoEDocs";
import ProDocs from "../doc_content/explorer_content/ProDocs";
import PolyDocs from "../doc_content/explorer_content/PolyDocs";
import SSComponent from "../../pages/Discriminator/SimilaritySearch/SSearch";
import { SimilaritySearchTabContents } from "./NestDiscriminatorTags";
import { desc } from "framer-motion/client";

const DiscriminatorTabContents = [
    {
        id: 1,
        title: 'Similarity Search',
        link: '',
        description: 'Search for similar molecules and proteins.',
        component: <SSComponent tabContent={SimilaritySearchTabContents} basePath="discriminators" />,
        docs: null,
        includeDocs: true,
        subElements: SimilaritySearchTabContents
    }
]

export default DiscriminatorTabContents;