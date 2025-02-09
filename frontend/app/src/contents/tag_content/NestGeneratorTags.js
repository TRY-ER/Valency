import PSMILESComponent from "../../pages/Generator/BRICS/PSMILESContent";
import ProDocs from "../doc_content/explorer_content/ProDocs";
import BRICSComponent from "../../pages/Generator/BRICS/PSMILESContent";
import LSTMComponent from "../../pages/Generator/LSTM/LSTMContent";
import { desc } from "framer-motion/client";

const BRICSGeneratorTabContents = [
    {
        id: 1,
        title: 'SMILES',
        link: '',
        description: 'Generate molecules using SMILES strings.',
        component: <BRICSComponent key="smiles" inputType="smiles"/>, 
        docs: <ProDocs />,
    },
    {
        id: 2,
        title: 'PSMILES',
        link: 'psmiles',
        description: 'Generate polymers using PSMILES strings.',
        component: <BRICSComponent key="psmiles" inputType="psmiles"/>, 
        docs: <ProDocs />,
    }
];

const LSTMGeneratorTabContents = [
    {
        id: 1,
        title: 'PSMILES',
        link: '',
        description: 'Generate polymers using PSMILES strings.',
        component: <LSTMComponent key="psmiles" />, 
        includeDocs: true, 
        docs: <ProDocs />,
    },
    {
        id: 2,
        title: 'WDG',
        link: 'moe',
        description: 'Generate molecules using WDG strings.',
        component: <LSTMComponent inputType="wdg" key="wdg" />,
        docs: <ProDocs />,
    }
];



export {
    BRICSGeneratorTabContents,
    LSTMGeneratorTabContents
}