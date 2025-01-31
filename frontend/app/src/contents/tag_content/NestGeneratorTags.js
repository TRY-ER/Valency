import PSMILESComponent from "../../pages/Generator/BRICS/PSMILESContent";
import ProDocs from "../doc_content/explorer_content/ProDocs";
import BRICSComponent from "../../pages/Generator/BRICS/PSMILESContent";
import LSTMComponent from "../../pages/Generator/LSTM/LSTMContent";

const BRICSGeneratorTabContents = [
    {
        id: 1,
        title: 'SMILES',
        link: '/proe',
        component: <BRICSComponent key="smiles" inputType="smiles"/>, 
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

const LSTMGeneratorTabContents = [
    {
        id: 1,
        title: 'PSMILES',
        link: '/proe',
        component: <LSTMComponent key="psmiles" />, 
        docs: <ProDocs />,
    },
    {
        id: 2,
        title: 'WDG',
        link: '/moe',
        component: <LSTMComponent inputType="wdg" key="wdg" />,
        docs: <ProDocs />,
    }
];



export {
    BRICSGeneratorTabContents,
    LSTMGeneratorTabContents
}