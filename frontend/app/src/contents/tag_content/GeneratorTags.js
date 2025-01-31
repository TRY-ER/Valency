import MoEDocs from "../doc_content/explorer_content/MoEDocs";
import ProDocs from "../doc_content/explorer_content/ProDocs";
import PolyDocs from "../doc_content/explorer_content/PolyDocs";

import ProtEComponent from "../../pages/Explorer/ProtExplorer/ProtExplorer";
import BRICSGenerator from "../../pages/Generator/BRICS/BRICSGenerator";
import LSTMGenerator from "../../pages/Generator/LSTM/LSTMGenerator";

const GeneratorTabContents = [
    {
        id: 1,
        title: 'BRICS Generator',
        link: '/moe',
        component: <BRICSGenerator />, 
        docs: null,
    },
    {
        id: 2,
        title: 'LSTM Generator',
        link: '/proe',
        component: <LSTMGenerator />, 
        docs: null,
    },
    {
        id: 3,
        title: 'LMLF Generator',
        link: '/proe',
        component: null, 
        docs: null,
    },
];

export default GeneratorTabContents;