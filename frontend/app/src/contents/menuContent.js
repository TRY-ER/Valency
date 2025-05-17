import Discriminator from "../pages/Discriminator/Discriminator";
import Explorer from "../pages/Explorer/Explorer";
import Generator from "../pages/Generator/Generator";
import ChatInterface from "../components/chat_interface/ChatInterface";
import ToolInterface from "../components/tool_interface/ToolInterface";
import GeneratorTabContents from "./tag_content/GeneratorTags";
import ExplorerTabContent from "./tag_content/ExploreTags";
import DiscriminatorTabContents from "./tag_content/DiscriminatorTags";
import { FaCompass, FaBalanceScale, FaFlask, FaComments, FaTools } from 'react-icons/fa';

const menuContent = [
    {
        id: 1,
        title: 'Explorers',
        icon: <FaCompass />,
        link: '',
        description: `Explorers are used to explore molecules, proteins, and polymers.
                      This is the main tool below which there are three more tools consisting
                      tools specfic to molecules, proteins and polymers`,
        component: <Explorer tabContent={ExplorerTabContent} />,
        includeDocs: true,
        subElements: ExplorerTabContent
    },
    {
        id: 2,
        title: 'Discriminators',
        icon: <FaBalanceScale />,
        link: 'discriminators',
        description: `Discriminators tool cosists of tools to discriminate molecules, proteins, and polymers.
                      But as of now we have similarity search implemented as a sub tool. This sub tool contains
                      structural similarity search for molecules, proteins, and polymers handled as sub-sub tools.`,
        includeDocs: true,
        component: <Discriminator tabContent={DiscriminatorTabContents} basePath="discriminators" />,
        subElements: DiscriminatorTabContents
    },
    {
        id: 3,
        title: 'Generators',
        icon: <FaFlask />,
        link: 'generators',
        description: `Generators are used to create hypothetical molecules, polymers (protein generation is not yet implemented).
                      It has two sub-tools name BRICS and LSTM generators. Both of these tools have sub tools for generation for SMILES
                      and PSMILES format for molecules and polymers respectively.`,
        includeDocs: true,
        component: <Generator tabContent={GeneratorTabContents} basePath="generators" />,
        subElements: GeneratorTabContents
    },
    {
        id: 4,
        title: 'Chatbot',
        icon: <FaComments />,
        link: 'chatbot',
        description: 'Chatbot is used to craftfully navigate in the system.',
        component: <ChatInterface />
    },
    {
        id: 5,
        title: 'Tool Config',
        icon: <FaTools />,
        description: 'Tool Config is used to configure the system.',
        link: 'tool-config',
        component: <ToolInterface />
    }
]

export default menuContent;