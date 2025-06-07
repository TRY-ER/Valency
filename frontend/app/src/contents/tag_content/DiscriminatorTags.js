import SSComponent from "../../pages/Discriminator/SimilaritySearch/SSearch";
import { SimilaritySearchTabContents, ChemBLTabContents, PubChemTabContents } from "./NestDiscriminatorTags";
import ChemBLComponent from "../../pages/Discriminator/ChemBL/ChemBL";
import PubChemComponent from "../../pages/Discriminator/PubChem/PubChem";

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
        title: 'ChEMBL',
        link: 'chembl',
        description: `This tool uses ChemBL database to conduct a similarity search`,
        component: <ChemBLComponent tabContent={ChemBLTabContents} basePath="identification/chembl" />,
        docs: null,
        includeDocs: true,
        subElements: ChemBLTabContents 
    },
    {
        id: 3,
        title: 'PubChem',
        link: 'pubchem',
        description: `This tool uses PubChem database to retrieve compound information, conduct similarity searches, and access various chemical utilities including substructure search, mass search, cross-reference lookup, and property calculation.`,
        component: <PubChemComponent tabContent={PubChemTabContents} basePath="identification/pubchem" />,
        docs: null,
        includeDocs: true,
        subElements: PubChemTabContents 
    }
]

export default DiscriminatorTabContents;