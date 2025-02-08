import PSMILESComponent from "../../pages/Generator/BRICS/PSMILESContent";
import ProDocs from "../doc_content/explorer_content/ProDocs";
import BRICSComponent from "../../pages/Generator/BRICS/PSMILESContent";
import SSContent from "../../pages/Discriminator/SimilaritySearch/SSContent";

const SimilaritySearchTabContents = [
    {
        id: 1,
        title: 'SMILES',
        link: '',
        component: <SSContent key="smiles" inputType={"MOL"} />, 
        docs: <ProDocs  />,
    },
    {
        id: 2,
        title: 'PSMILES',
        link: 'psmiles',
        component:  <SSContent key="psmiles" inputType={"POLY"} />,
        docs: <ProDocs />,
    },
    {
        id: 3,
        title: 'Protein',
        link: 'protein',
        component:  <SSContent key="pdb" inputType={"PROT"} />,
        docs: <ProDocs />,
    }
];

export {
    SimilaritySearchTabContents,
}