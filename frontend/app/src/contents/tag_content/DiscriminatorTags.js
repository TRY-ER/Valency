import SSComponent from "../../pages/Discriminator/SimilaritySearch/SSearch";
import { SimilaritySearchTabContents, ChemBLTabContents } from "./NestDiscriminatorTags";
import ChemBLComponent from "../../pages/Discriminator/ChemBL/ChemBL";

const DiscriminatorTabContents = [
    {
        id: 1,
        title: 'Similarity Search',
        link: '',
        description: `This tool uses vectors generated from SMILES, PSMILES, Protein sequences to conduct a similarity search
                      and get relevant candidates from those. This tool contains three sub-tools for similarity search with 
                      SMILES for molecules, PSMILES for polymers, and PDB-ID for proteins`,
        component: <SSComponent tabContent={SimilaritySearchTabContents} basePath="identification" />,
        docs: null,
        includeDocs: true,
        subElements: SimilaritySearchTabContents
    },
    {
        id: 2,
        title: 'ChemBL',
        link: 'chembl',
        description: `This tool uses ChemBL database to conduct a similarity search`,
        component: <ChemBLComponent tabContent={ChemBLTabContents} basePath="identification/chembl" />,
        docs: null,
        includeDocs: true,
        subElements: ChemBLTabContents 
    }
]

export default DiscriminatorTabContents;