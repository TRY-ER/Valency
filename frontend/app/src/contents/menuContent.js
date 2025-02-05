import Discriminator from "../pages/Discriminator/Discriminator";
import Explorer from "../pages/Explorer/Explorer";
import Generator from "../pages/Generator/Generator";
import ChatInterface from "../components/chat_interface/ChatInterface";

const menuContent = [
    {
        id: 1,
        title: 'Explorers',
        iconPath: 'images/explorer.png',
        link: '/',
        component: <Explorer /> 
    },
    {
        id: 2,
        title: 'Discriminators',
        iconPath: 'images/discriminator.png',
        link: '/discriminators',
        component: <Discriminator /> 
    },
    {
        id: 3,
        title: 'Generators',
        iconPath: 'images/generator.png',
        link: '/generators',
        component: <Generator />
    },
    {
        id: 4,
        title: 'Chatbot',
        iconPath: 'images/chatbot.png',
        link: '/chatbot',
        component: <ChatInterface /> 
    }
]

export default menuContent;