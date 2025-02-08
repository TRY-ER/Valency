import MoEDocs from "../doc_content/explorer_content/MoEDocs";
import ProDocs from "../doc_content/explorer_content/ProDocs";
import PolyDocs from "../doc_content/explorer_content/PolyDocs";

import ProtEComponent from "../../pages/Explorer/ProtExplorer/ProtExplorer";
import BRICSGenerator from "../../pages/Generator/BRICS/BRICSGenerator";
import LSTMGenerator from "../../pages/Generator/LSTM/LSTMGenerator";
import { BRICSGeneratorTabContents, LSTMGeneratorTabContents } from "./NestGeneratorTags";

const GeneratorTabContents = [
    {
        id: 1,
        title: 'BRICS Generator',
        link: '',
        component: <BRICSGenerator tabContent={BRICSGeneratorTabContents} basePath="generators" />, 
        subElements: BRICSGeneratorTabContents,
        includeDocs: false,
        docs: <ProDocs />,
    },
    {
        id: 2,
        title: 'LSTM Generator',
        link: 'lstm',
        component: <LSTMGenerator tabContent={LSTMGeneratorTabContents} basePath="generators/lstm" />, 
        subElements: LSTMGeneratorTabContents,
        includeDocs: false,
        docs: <ProDocs />,
    },
    {
        id: 3,
        title: 'LMLF Generator',
        link: 'lmlf',
        component: null, 
        docs: null,
    },
];

export default GeneratorTabContents;