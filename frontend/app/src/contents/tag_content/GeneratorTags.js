import BRICSGenerator from "../../pages/Generator/BRICS/BRICSGenerator";
import LSTMGenerator from "../../pages/Generator/LSTM/LSTMGenerator";
import { BRICSGeneratorTabContents, LSTMGeneratorTabContents } from "./NestGeneratorTags";

const GeneratorTabContents = [
    {
        id: 1,
        title: 'BRICS Generator',
        link: '',
        description: `Generate molecules using the BRICS algorithm. The BRICS algorithm is a rule-based algorithm that generates molecules based on the given input.
                      The input can be a SMILES string or a PSMILES string. The output will be a set of molecules that are generated based on the input.
                      This tool has two sub-tools for generation for SMILES and PSMILES format for molecules and polymers respectively.`,
        component: <BRICSGenerator tabContent={BRICSGeneratorTabContents} basePath="generators" />, 
        subElements: BRICSGeneratorTabContents,
        includeDocs: true,
        docs: null,
    },
    {
        id: 2,
        title: 'LSTM Generator',
        link: 'lstm',
        description: `Generate molecules using the LSTM algorithm. The LSTM algorithm is a deep learning algorithm that generates molecules based on the given input.
                      The input can be a SMILES string or a PSMILES string. The output will be a set of molecules that are generated based on the input.
                      This tool has two sub-tools for generation for SMILES and PSMILES format for molecules and polymers respectively.`,
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