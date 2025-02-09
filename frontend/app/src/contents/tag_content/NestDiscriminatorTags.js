import PSMILESComponent from "../../pages/Generator/BRICS/PSMILESContent";
import ProDocs from "../doc_content/explorer_content/ProDocs";
import BRICSComponent from "../../pages/Generator/BRICS/PSMILESContent";
import SSContent from "../../pages/Discriminator/SimilaritySearch/SSContent";
import { desc } from "framer-motion/client";

const SimilaritySearchTabContents = [
    {
        id: 1,
        title: 'SMILES',
        link: '',
        description: 'Search for similar molecules using SMILES strings.',
        component: <SSContent key="smiles" inputType={"MOL"} />, 
        docs: <ProDocs />,
    },
    {
        id: 2,
        title: 'PSMILES',
        link: 'psmiles',
        description: 'Search for similar polymers using PSMILES strings.',
        component:  <SSContent key="psmiles" inputType={"POLY"} />,
        docs: <ProDocs />,
    },
    {
        id: 3,
        title: 'Protein',
        link: 'protein',
        description: 'Search for similar proteins using Protein Data Bank identifiers.',
        component:  <SSContent key="pdb" inputType={"PROT"} />,
        docs: <ProDocs />,
    }
]; 

export {
    SimilaritySearchTabContents,
}