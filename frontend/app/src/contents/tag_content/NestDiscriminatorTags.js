import PSMILESComponent from "../../pages/Generator/BRICS/PSMILESContent";
import ProDocs from "../doc_content/explorer_content/ProDocs";
import BRICSComponent from "../../pages/Generator/BRICS/PSMILESContent";
import SSContent from "../../pages/Discriminator/SimilaritySearch/SSContent";

const SimilaritySearchTabContents = [
    {
        id: 1,
        title: 'SMILES',
        link: '/proe',
        component: <SSContent />, 
        docs: <ProDocs />,
    },
    {
        id: 2,
        title: 'PSMILES',
        link: '/moe',
        component: <BRICSComponent key="psmiles" inputType="psmiles"/>, 
        docs: <ProDocs />,
    }
];

export {
    SimilaritySearchTabContents,
}