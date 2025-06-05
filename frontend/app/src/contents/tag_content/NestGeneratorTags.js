import BRICSComponent from "../../pages/Generator/BRICS/PSMILESContent";
import LSTMComponent from "../../pages/Generator/LSTM/LSTMContent";
import DocRenderer from "../doc_content/DocRenderer";
import PSMILESListComponent from "../../pages/Generator/BRICS/PSMILESListComponent";

const BRICSGeneratorTabContents = [
    {
        id: 1,
        title: 'SMILES',
        link: '',
        description: `This tool uses SMILES string to generate molecules using BRICS algorithm. The algorithm generates hypothetical molecules based on the given input.
                      The input can be a SMILES string in a text file each represented in a new line. The output will be a set of SMILES that are generated 
                      based on the input. The tool has three sections, the first section shows a progress bar showing the completed steps for the job.
                      The second section contains the file upload setup with a panel to show first few files of file to ensure the right text file has been selected.
                      The third section is the output section that consists of a output panel, download button and the reset button.
                      the output panel shows the steps details as well as completed steps and the molecules generated. This tool is specifically useful when you have
                      a set of molecules and you want to generate more molecules based on the given set.`,
        component: <BRICSComponent key="smiles" inputType="smiles"/>, 
        docs: <DocRenderer filePath="/markdown_repo/BRICSSMILES.md"/>,
    },
    {
        id: 2,
        title: 'PSMILES',
        link: 'psmiles',
        description: `This tool uses PSMILES string to generate molecules using BRICS algorithm. The algorithm generates hypothetical polymers based on the given input.
                      The input can be a PSMILES string in a text file each represented in a new line. The output will be a set of PSMILES that are generated 
                      based on the input. The tool has three sections, the first section shows a progress bar showing the completed steps for the job.
                      The second section contains the file upload setup with a panel to show first few files of file to ensure the right text file has been selected.
                      The third section is the output section that consists of a output panel, download button and the reset button.
                      the output panel shows the steps details as well as completed steps and the psmiles generated upon completed. This tool is specifically useful when you have
                      a set of polymers and you want to generate more polymers based on the given set.`,
        component: <BRICSComponent key="psmiles" inputType="psmiles"/>, 
        docs: <DocRenderer filePath="/markdown_repo/BRICSPSMILES.md"/>,
    },
    {
        id: 3,
        title: 'Advanced',
        link: 'adv',
        description: `This tool uses SMILES string to generate molecules using LSTM algorithm. The algorithm generates hypothetical molecules without any given input molecule string.`,
        component: <PSMILESListComponent />,
        docs: null,
    }
];

const LSTMGeneratorTabContents = [
    {
        id: 1,
        title: 'PSMILES',
        link: '',
        description: `This tool uses PSMILES string to generate molecules using LSTM algorithm. The algorithm generates hypothetical polymers without any given input polymer string.
                      It's very useful when you want to consider entirely new polymer rather refering from a specific kind. The tool takes input integer for number of generations.
                      Then it generates the polymers based on the given number of generations. The tool has three sections, the first section shows a progress bar showing the number 
                      of genertions recieved. The second section contains the input setup with a input box to take the number of generations. The third section is the output section that
                      consists of a output panel, download button and the reset button. the output panel shows the steps details as well as completed steps and the psmiles generated upon completed.`,
        component: <LSTMComponent key="psmiles" />, 
        includeDocs: true, 
        docs: <DocRenderer filePath="/markdown_repo/LSTMPSMILES.md"/>,
    },
    {
        id: 2,
        title: 'WDG',
        link: 'wdg',
        description: `This tool uses weighted directed graph string to generate molecules using LSTM algorithm. The algorithm generates hypothetical polymers without any given input polymer string.
                      It's very useful when you want to consider entirely new polymer rather refering from a specific kind. The tool takes input integer for number of generations.
                      Then it generates the polymers based on the given number of generations. The tool has three sections, the first section shows a progress bar showing the number 
                      of genertions recieved. The second section contains the input setup with a input box to take the number of generations. The third section is the output section that
                      consists of a output panel, download button and the reset button. the output panel shows the steps details as well as completed steps and the psmiles generated upon completed.`,
        component: <LSTMComponent inputType="wdg" key="wdg" />,
        docs: <DocRenderer filePath="/markdown_repo/LSTMWDG.md"/>,
    }
];


export {
    BRICSGeneratorTabContents,
    LSTMGeneratorTabContents
}