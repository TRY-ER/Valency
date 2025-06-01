import Discriminator from "../pages/Discriminator/Discriminator";
import Explorer from "../pages/Explorer/Explorer";
import Generator from "../pages/Generator/Generator";
import ChatInterface from "../components/chat_interface/ChatInterface";
import ToolInterface from "../components/tool_interface/ToolInterface";
import GeneratorTabContents from "./tag_content/GeneratorTags";
import ExplorerTabContent from "./tag_content/ExploreTags";
import DiscriminatorTabContents from "./tag_content/DiscriminatorTags";
import {
    FaCompass,
    FaBalanceScale,
    FaFlask,
    FaComments,
    FaTools,
    FaHome,
    FaTasks,
    FaRobot,
    FaUserSecret // Added FaUserSecret
} from 'react-icons/fa';

// Placeholder components - consider moving these to their own files
const HomeComponent = () => <div style={{ padding: '20px' }}><h2>Home Page</h2><p>Welcome to the Drug Discovery Toolkit. This is the central hub for all tools and features.</p></div>;
const ActivitiesPage = () => <div style={{ padding: '20px' }}><h2>Activities Page</h2><p>This page will be used to manage and track drug discovery activities. Content coming soon.</p></div>;
const AutomationPage = () => <div style={{ padding: '20px' }}><h2>Automation Page</h2><p>This page will be used for automating drug discovery workflows. Content coming soon.</p></div>;

const menuContent = [
    {
        id: 0,
        title: 'Home',
        icon: <FaHome />,
        link: '', // Root path
        description: 'Overview and starting point of the drug discovery toolkit.',
        component: <HomeComponent />,
    },
    {
        id: 1,
        title: 'Structure Analysis', // Renamed from 'Explorers'
        icon: <FaCompass />,
        link: 'explorer', // Changed link from '' to 'explorer'
        description: `Tools for analyzing the structure of molecules, proteins, and polymers.
                      This section includes sub-tools specific to molecules, proteins, and polymers.`,
        component: <Explorer tabContent={ExplorerTabContent} basePath="explorer" />, // Added basePath
        includeDocs: true,
        subElements: ExplorerTabContent
    },
    {
        id: 2,
        title: 'Identification', // Renamed from 'Discriminators'
        icon: <FaBalanceScale />,
        link: 'identification', // Changed link from 'discriminators'
        description: `Tools for identifying and discriminating molecules, proteins, and polymers.
                      Currently, this includes similarity search functionalities.`,
        includeDocs: true,
        component: <Discriminator tabContent={DiscriminatorTabContents} basePath="identification" />, // Changed basePath
        subElements: DiscriminatorTabContents
    },
    {
        id: 3,
        title: 'Optimization', // Renamed from 'Generators'
        icon: <FaFlask />,
        link: 'optimization', // Changed link from 'generators'
        description: `Tools for generating and optimizing hypothetical molecules and polymers.
                      It features sub-tools like BRICS and LSTM generators for SMILES and PSMILES formats.`,
        includeDocs: true,
        component: <Generator tabContent={GeneratorTabContents} basePath="optimization" />, // Changed basePath
        subElements: GeneratorTabContents
    },
    {
        id: 4,
        title: 'Activities', // New dummy page
        icon: <FaTasks />,
        link: 'activities',
        description: 'Placeholder for managing and tracking drug discovery activities.',
        component: <ActivitiesPage />
    },
    {
        id: 5,
        title: 'Automation', // New dummy page
        icon: <FaRobot />,
        link: 'automation',
        description: 'Placeholder for automating drug discovery workflows.',
        component: <AutomationPage />
    },
    {
        id: 6, // ID shifted from 4
        title: 'Master Agent',
        icon: <FaUserSecret />, // Changed from FaComments to FaUserSecret
        link: 'chatbot',
        description: 'Chatbot is used to craftfully navigate in the system.',
        component: <ChatInterface />
    },
    {
        id: 7, // ID shifted from 5
        title: 'Tool Config',
        icon: <FaTools />,
        description: 'Tool Config is used to configure the system.',
        link: 'tool-config',
        component: <ToolInterface />
    }
]

export default menuContent;