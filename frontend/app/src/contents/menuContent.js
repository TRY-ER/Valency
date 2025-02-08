import Discriminator from "../pages/Discriminator/Discriminator";
import Explorer from "../pages/Explorer/Explorer";
import Generator from "../pages/Generator/Generator";
import ChatInterface from "../components/chat_interface/ChatInterface";
import ToolInterface from "../components/tool_interface/ToolInterface";
import MolEComponent from "../pages/Explorer/MolExplorer/MolExplorer";
import ProtEComponent from "../pages/Explorer/ProtExplorer/ProtExplorer";
import PolyEComponent from "../pages/Explorer/PolyExplorer/PolyExplorer";
import MoEDocs from "./doc_content/explorer_content/MoEDocs";
import ProDocs from "./doc_content/explorer_content/ProDocs";
import PolyDocs from "./doc_content/explorer_content/PolyDocs";
import SSComponent from "../pages/Discriminator/SimilaritySearch/SSearch";
import SSContent from "../pages/Discriminator/SimilaritySearch/SSContent";
import GeneratorTabContents from "./tag_content/GeneratorTags";

const ExplorerTabContent = [
    {
        id: 1,
        title: 'Molecule Explorer',
        link: '',
        component: <MolEComponent />,
        docs: <MoEDocs />,
    },
    {
        id: 2,
        title: 'Protein Explorer',
        link: 'proe',
        component: <ProtEComponent />,
        docs: <ProDocs />,
    },
    {
        id: 3,
        title: 'Polymer Explorer',
        link: 'polye',
        component: <PolyEComponent />,
        docs: <PolyDocs />,
    }
]

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

const DiscriminatorTabContents = [
    {
        id: 1,
        title: 'Similarity Search',
        link: '',
        component: <SSComponent tabContent={SimilaritySearchTabContents} basePath="discriminators" />,
        docs: null,
        includeDocs: true,
        subElements: SimilaritySearchTabContents
    },
    {
        id: 2,
        title: 'Similarity Search',
        link: 'temp',
        component: <SSComponent tabContent={SimilaritySearchTabContents} basePath="discriminators/temp" />,
        docs: null,
        includeDocs: true,
        subElements: SimilaritySearchTabContents
    }
]

const menuContent = [
    {
        id: 1,
        title: 'Explorers',
        iconPath: 'images/explorer.png',
        link: '',
        component: <Explorer tabContent={ExplorerTabContent} />,
        includeDocs: true,
        subElements: ExplorerTabContent
    },
    {
        id: 2,
        title: 'Discriminators',
        iconPath: 'images/discriminator.png',
        link: 'discriminators',
        component: <Discriminator tabContent={DiscriminatorTabContents} basePath="discriminators" />,
        subElements: DiscriminatorTabContents 
    },
    {
        id: 3,
        title: 'Generators',
        iconPath: 'images/generator.png',
        link: 'generators',
        component: <Generator tabContent={GeneratorTabContents} basePath="generators" />,
        subElements: GeneratorTabContents
    },
    {
        id: 4,
        title: 'Chatbot',
        iconPath: 'images/chatbot.png',
        link: 'chatbot',
        component: <ChatInterface />
    },
    {
        id: 5,
        title: 'Tool Config',
        iconPath: 'images/tool_config.png',
        link: 'tool-config',
        component: <ToolInterface />
    }
]

export default menuContent;