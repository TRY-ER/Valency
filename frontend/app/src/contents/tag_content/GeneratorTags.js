import MoEDocs from "../doc_content/explorer_content/MoEDocs";
import ProDocs from "../doc_content/explorer_content/ProDocs";
import PolyDocs from "../doc_content/explorer_content/PolyDocs";

import ProtEComponent from "../../pages/Explorer/ProtExplorer/ProtExplorer";
import BRICSGenerator from "../../pages/Generator/BRICS/BRICSGenerator";
import LSTMGenerator from "../../pages/Generator/LSTM/LSTMGenerator";
import { BRICSGeneratorTabContents, LSTMGeneratorTabContents } from "./NestGeneratorTags";
import { desc } from "framer-motion/client";

const GeneratorTabContents = [
    {
        id: 1,
        title: 'BRICS Generator',
        link: '',
        description: 'Generate molecules using the BRICS algorithm.',
        component: <BRICSGenerator tabContent={BRICSGeneratorTabContents} basePath="generators" />, 
        subElements: BRICSGeneratorTabContents,
        includeDocs: true,
        docs: null,
    },
    {
        id: 2,
        title: 'LSTM Generator',
        link: 'lstm',
        description: 'Generate molecules using the LSTM algorithm.',
        component: <LSTMGenerator tabContent={LSTMGeneratorTabContents} basePath="generators/lstm" />, 
        subElements: LSTMGeneratorTabContents,
        includeDocs: true,
        docs: null,
    },
    // {
    //     id: 3,
    //     title: 'LMLF Generator',
    //     link: 'lmlf',
    //     component: null, 
    //     docs: null,
    // },
];

export default GeneratorTabContents;