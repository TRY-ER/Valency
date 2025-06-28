import BRICSGenerator from "../../pages/Generator/BRICS/BRICSGenerator";
import LSTMGenerator from "../../pages/Generator/LSTM/LSTMGenerator";
import ADMETComponent from "../../pages/Generator/ADMET/ADMET";
import { BRICSGeneratorTabContents, LSTMGeneratorTabContents } from "./NestGeneratorTags";
import DocRenderer from "../doc_content/DocRenderer";

const GeneratorTabContents = [
    {
        id: 1,
        title: 'BRICS Generator',
        link: '',
        description: `Generate molecules using the BRICS algorithm. The BRICS algorithm is a rule-based algorithm that generates molecules based on the given input.
                      The input can be a SMILES string or a PSMILES string. The output will be a set of molecules that are generated based on the input.
                      This tool has two sub-tools for generation for SMILES and PSMILES format for molecules and polymers respectively.`,
        component: <BRICSGenerator tabContent={BRICSGeneratorTabContents} basePath="optimization" />, 
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
        component: <LSTMGenerator tabContent={LSTMGeneratorTabContents} basePath="optimization/lstm" />, 
        subElements: LSTMGeneratorTabContents,
        includeDocs: true,
        docs: null,
    },
    {
        id: 3,
        title: 'ADMET',
        link: 'admet',
        description: `Predict Absorption, Distribution, Metabolism, Excretion, and Toxicity (ADMET) properties of molecules using machine learning models.
                      Enter a molecular SMILES string to get predictions for how well the molecule will be absorbed, distributed throughout the body,
                      metabolized, excreted, and its potential toxicity. This tool helps in early drug discovery by identifying promising compounds
                      with favorable pharmacokinetic and safety profiles.`,
        component: <ADMETComponent />,   
        docs: <DocRenderer filePath="/markdown_repo/ADMETPredictor.md"/>,
        includeDocs: true,
    },
];

export default GeneratorTabContents;