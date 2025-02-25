import SSComponent from "../../pages/Discriminator/SimilaritySearch/SSearch";
import { SimilaritySearchTabContents } from "./NestDiscriminatorTags";

const DiscriminatorTabContents = [
    {
        id: 1,
        title: 'Similarity Search',
        link: '',
        description: `This tool uses vectors generated from SMILES, PSMILES, Protein sequences to conduct a similarity search
                      and get relevant candidates from those. This tool contains three sub-tools for similarity search with 
                      SMILES for molecules, PSMILES for polymers, and PDB-ID for proteins`,
        component: <SSComponent tabContent={SimilaritySearchTabContents} basePath="discriminators" />,
        docs: null,
        includeDocs: true,
        subElements: SimilaritySearchTabContents
    }
]

export default DiscriminatorTabContents;