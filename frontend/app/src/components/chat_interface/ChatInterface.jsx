import React, { useEffect, useState, useRef, use } from "react";
import { useNavigate, useLocation } from "react-router-dom";
import axios from "axios";
import { v4 as uuidv4 } from "uuid";
import { EventSourcePolyfill } from "event-source-polyfill";
import clipboard from "clipboard-copy";
import { removeHtmlTags, responseParser } from "../../components/parsers/Basic";
import { motion } from "framer-motion";
import "./ChatInterface.css";
import { IoSend } from "react-icons/io5";
import { FaTrash, FaUserCircle } from "react-icons/fa"; // Added import
import { chat_endpoint } from "../../endpoints/endpoints";
import { baseURL } from "../../endpoints/base";
import { call_endpoint_async, call_eventsource } from "../../endpoints/caller";

import { MODEL_OPTIONS_TYPES } from "./extendedVariables";
import SpinningLoader from "../UI/SpinLoader/SpinLoader";
import CustomToggler from "../UI/CustomToggler/CustomToggler";
import CustomDropdown from "../UI/CustomDropdown/CustomDropdown";
import GlassyContainer from "../glassy_container/gc";
import Markdown from "react-markdown";
import CustomMarkdownRenderer from "../UI/CustomMarkdown";
import FlexRenderer from "../UI/FlexRenderer"; // Added FlexRenderer import
import AuthService from "../../services/api/AuthService.ts"; // Added AuthService import
import {
    pingAgent,
    testAgentEndpoint,
    createUserSession,
    listUserSessions,
    getSessionDetails, // Added
    sendQueryToSession, // Added
    getSessionEventState // Added for loading session event state
} from '../../services/api/agentService';

const ChatInterface = () => {
    const [chatStream, setChatStream] = useState(''); // This state might become redundant
    const [messages, setMessages] = useState([]);
    const location = useLocation(); // Add this to access URL query parameters
    const [currentSessionId, setCurrentSessionId] = useState(null);
    const [sessions, setSessions] = useState([]);
    const [isLoadingSession, setIsLoadingSession] = useState(false);
    const [generationState, setGenerationState] = useState("default"); // init, complete, default, generating (four value can be given)
    const [isContentLoading, setIsContentLoading] = useState(true);
    const [query, setQuery] = useState("");
    const navigate = useNavigate();
    const bottomRef = useRef(null);
    const [serverStreamId, setServerStreamId] = useState(null); // This state might become redundant
    const [isLoaded, setIsLoaed] = useState(false);
    const [isStreaming, setIsStreaming] = useState(false); // This state might become redundant
    const [isCompleted, setIsCompleted] = useState(false); // This state might become redundant
    const [isQuaryParamLoaded, setIsQuaryParamLoaded] = useState(false); // This state might become redundant


    //ui changes
    // const [isDark, setIsDark] = useState(false); // REMOVED
    const [fontSize, setFontSize] = useState("14");


    //info-tag setup
    const [showTag, setShowTag] = useState(false);
    const [tagText, setTagText] = useState("test");
    const [tagPosition, setTagPosition] = useState({ x: 0, y: 0 });

    const [username, setUsername] = useState(""); // Added state for username

    // input config for model
    const [modelConfig, setModelConfig] = useState({
        "model_name": "Gemini 2 Flash",
        "model_code": "gemini-2.5-flash-exp",
    })

    useEffect(() => {
        console.log('current session id >>', currentSessionId);
    }, [currentSessionId])

    // Added useEffect to check for session_id in URL query parameters
    useEffect(() => {
        const queryParams = new URLSearchParams(location.search);
        const sessionIdFromUrl = queryParams.get('session_id');

        if (sessionIdFromUrl) {
            console.log('Found session_id in URL:', sessionIdFromUrl);
            setCurrentSessionId(sessionIdFromUrl);
            setIsQuaryParamLoaded(true); // Set to true when session_id is found
        } else {
            console.log('No session_id found in URL, using default null value');
            setCurrentSessionId(null);
        }
    }, [location.search]); // Re-run if URL query parameters change

    useEffect(() => { // Added useEffect to get username
        const currentUser = AuthService.getUser();
        if (currentUser && currentUser.username) {
            setUsername(currentUser.username);
        }
    }, [isLoaded]); // Re-run if isLoaded changes, ensuring user data might be available after initial load logic

    const handleShowTag = (e, infoText, horizontalOffset = 5) => {
        setShowTag(true);
        const currentRect = e?.target.getBoundingClientRect();
        setTagText(infoText);
        setTagPosition({ x: currentRect.left - (currentRect.width / 2) + horizontalOffset, y: currentRect.top - (currentRect.height / 2) - 5 })
    }

    const handleHideTag = () => {
        setTagText("");
        setShowTag(false);
    }

    //states to place the scroll to bottom button
    const chatContainerRef = useRef(null);
    const scrollBottomRef = useRef(null);

    // component states 
    const [components, setComponents] = useState([]);
    const [bottomScroller, setBottomScroller] = useState(false);
    const [activeUid, setActiveUid] = useState("");

    // Define suggestion queries
    const suggestionQueries = [
        "Can you provide a detailed explanation of how the BRICSPSMILES tool works and its primary applications in cheminformatics?",
        "What are the latest advancements in molecule exploration techniques, and how do they compare to traditional methods?",
        "Tell me about the different types of улыбки (smiles) used in chemical informatics and their specific use cases.",
        "Explain the concept of Polymer Similarity Search and discuss some of the common algorithms or tools used for this purpose.",
        "What are the key challenges and future directions in the field of Protein Explorer tools and technologies?",
        "Describe the process of using the UniProt Viewer for accessing and analyzing protein sequence and functional information.",
        "How can machine learning models be applied to predict molecular properties, and what are some example use cases?",
        "Discuss the importance of data visualization in bioinformatics and provide examples of effective visualization techniques.",
        "What are the ethical considerations when working with large-scale biological datasets and AI in drug discovery?",
        "Can you explain the role of SMILES strings in representing chemical structures and their advantages over other notations?",
        "Explore the different databases available for protein structures and how to query them effectively.",
        "What is the significance of the FASTA format in bioinformatics, and how is it typically used?"
    ];

    // useEffect(() => {
    //     if (chatContainerRef.current && scrollBottomRef.current) {
    //         const chatContainerRect = chatContainerRef.current.getBoundingClientRect();
    //         const scrollBottomRect = scrollBottomRef.current.getBoundingClientRect();
    //         const newPosition = {
    //             bottom: chatContainerRect.bottom - scrollBottomRect.height,
    //             right: chatContainerRect.right - scrollBottomRect.width,
    //         }
    //         scrollBottomRef.current.style.transform = `translate(${newPosition.right}px, ${newPosition.bottom}px)`;
    //     }
    // }, [chatContainerRef, scrollBottomRef, components])

    const scrollToBottom = () => {
        if (chatContainerRef.current) {
            chatContainerRef.current.scrollTo({
                top: chatContainerRef.current.scrollHeight,
                behavior: "smooth",
            });
        }
    }

    useEffect(() => {
        if (bottomScroller) {
            scrollToBottom();
            setBottomScroller(false);
        }
    }, [bottomScroller, components])

    const handleKeyDown = (e) => {
        // setBottomScroller(true);
        if (e.key === "Enter" && !e.ctrlKey) {
            e.preventDefault();
            generateResponse();
        }
        else if (e.key === "Enter" && e.ctrlKey) {
            e.preventDefault();
            setQuery(query + "\n");
        }
    }

    const handleClipboard = async (index) => {
        if (components[index].content[0].response.res_content === "") {
            alert("nothing to copy")
        }
        else {
            await clipboard(removeHtmlTags(components[index].content[0].response.res_content)).then(() => {
                alert("text copied to the clipboard")
            }).catch(err => {
                alert(`failed to copy text ${err}`)
            })
        }
    }


    const loadThread = () => {
        // console.log("load thread is triggered")
        // setComponents(localStorage.getItem("components") ? JSON.parse(localStorage.getItem("components")) : []);
        // setModelConfig(localStorage.getItem("modelConfig") ? JSON.parse(localStorage.getItem("modelConfig")) : modelConfig);
        // setFontSize(localStorage.getItem("fontsize") ? localStorage.getItem("fontsize") : "14");

        // If we have a currentSessionId, fetch the session event state
        if (currentSessionId && isQuaryParamLoaded) {
            console.log('Fetching session event state for session:', currentSessionId);
            getSessionEventState(currentSessionId)
                .then(eventState => {
                    console.log('Retrieved session event state:', eventState);
                    
                    // Check if eventState is an array and has elements
                    if (Array.isArray(eventState) && eventState.length > 0) {
                        // Group events by query to create conversation pairs
                        const conversations = [];
                        // Process events and organize them into conversation components
                        for (let i = 0; i < eventState.length; i++) {
                            const event = eventState[i];
                            console.log("event >>", event)
                            if (event) {
                                const cUid = uuidv4();
                                conversations.push({
                                    "content": [event],
                                    "uid": cUid
                                });
                            } 
                        }
                        setComponents(conversations);
                        console.log('Created conversation components from event state:', conversations);
                    } else {
                        console.warn('Event state is empty or not an array:', eventState);
                    }
                })
                .catch(error => {
                    console.error('Error fetching session event state:', error);
                });
        } else {
            console.log('No current session ID available, skipping event state fetch');
        }

        // Session ID is now only maintained in React state
        // No need to load from localStorage

        // await axios.post(API_BASE + "/chat/load_thread/", {
        //     "user_id": userName
        // }, config)
        //     .then(res => {
        // console.log("load_response", res)
        // if (res?.data?.response?.type == "success") {
        //     let comp_objs = res?.data?.response?.data?.components;
        //     let summary = res?.data?.response?.data?.summary;
        //     let config = res?.data?.response?.data?.thread_config;
        //     if ("input_lang" in config) {
        //         console.log("config", config)
        //         const mod_config = {
        //             input_lang: config["input_lang"],
        //             output_lang: config["output_lang"],
        //             model_type: config["model_type"],
        //         }
        //         setModelConfig(mod_config);
        //     }
        //     setComponents(comp_objs);
        //     setBottomScroller(true);
        //     setChatSummary(summary);
        //     setIsContentLoading(false);
        // }
        // }).catch(err => {
        //     console.log("error", err)
        // })
        setIsLoaed(true);
    }

    useEffect(() => {
        loadThread();
    }, [currentSessionId]) // Added currentSessionId as a dependency to reload when it changes

    // useEffect(() => {
    //     console.log("components >>", components)
    // }, [components])


    const updateComponentResponse = (uid, newResponse) => {
        setComponents((prev) => {
            const new_value = prev.map((component) =>
                component.uid === uid ?
                    {
                        ...component,
                        content: [{
                            ...component.content[0], response: [
                                ...component.content[0].response,
                                newResponse
                            ]
                        }]
                    } : component
            );
            return new_value;
        });
    }

    const delHandler = (itemIndex) => {
        // console.log("item index", itemIndex)
        // console.log("del elem", components[itemIndex])
        const delItem = async (index) => {
            // await axios.post(`${API_BASE}` + "/chat/del_component", {
            //     "user_id": auth?.userName,
            //     "c_index": index
            // }, privateConfig).then((res) => {
            // console.log("del result", res)
            // if (res.data.type == "success") {
            const comp_copy = [...components];
            comp_copy.splice(index, 1);
            setComponents(comp_copy);
            // }
            // // }).catch(err => {
            // //     console.log("error in item deletion")
            // })
        }
        delItem(itemIndex);
    }

    useEffect(() => {
        console.log("components >>", components)
    }, [components])

    const generateResponse = async () => {
        const currentQuery = query;
        if (currentQuery.trim().length === 0) {
            alert("Enter some query first to generate response");
            return;
        }

        setGenerationState("generating");
        setQuery(""); // Clear input field

        const cUid = uuidv4();
        const new_comp = {
            "content": [{
                "query": currentQuery,
                "response": []
            }],
            "uid": cUid,
        };

        setComponents(prevComponents => [...prevComponents, new_comp]);
        setActiveUid(cUid);
        setBottomScroller(true); // Scroll after adding user query

        try {
            // Define common callbacks for both methods
            const callbacks = {
                onAgentMessage: (agentMessageData) => {
                    // Process the message received from the agent
                    console.log('agent messages received >>', agentMessageData);
                    // Make sure to properly format the data before passing it to updateComponentResponse
                    // The FlexRenderer component expects certain data structures based on type
                    // if (typeof agentMessageData === 'string') {
                    //     // If it's a plain string, wrap it as a text type
                    //     updateComponentResponse(cUid, { type: 'text', content: agentMessageData });
                    // } else if (agentMessageData.type === 'markdown_text') {
                    //     // Handle markdown_text type from _processSseStream
                    //     updateComponentResponse(cUid, { type: 'text', content: agentMessageData.content });
                    // } else {
                    //     // Pass other structured data as-is
                    // }
                    updateComponentResponse(cUid, agentMessageData);
                    setBottomScroller(true); // Scroll as new content arrives
                },
                onStreamEnd: (streamEndData) => {
                    console.log('Stream Ended:', streamEndData);
                    setGenerationState("default");
                    setActiveUid("");
                    const now = new Date();
                    updateComponentResponse(cUid, { type: "complete", data: now.toISOString() }); // Pass ISO string
                    console.log('Stream completed at:', now);
                },
                onProcessingError: (error) => {
                    console.error('Session Processing Error:', error);
                    const errorMessage = `Error: ${error.message || 'Processing failed'}`;
                    updateComponentResponse(cUid, { type: "error", data: errorMessage, timestamp: new Date().toISOString() });
                    setGenerationState("default");
                    setActiveUid("");
                },
                onSetupError: (error) => {
                    console.error('Session Setup Error:', error);
                    const setupErrorMessage = `Setup Error: ${error.message || 'Setup failed'}`;
                    updateComponentResponse(cUid, { type: "error", data: setupErrorMessage, timestamp: new Date().toISOString() });
                    setGenerationState("default");
                    setActiveUid("");
                }
            };

            // If we already have a session, use sendQueryToSession
            if (currentSessionId) {
                console.log('Using existing session:', currentSessionId);
                await sendQueryToSession(
                    currentSessionId,
                    currentQuery,
                    callbacks
                );
            } else {
                // For the first query, create a new session
                console.log('Creating new session for first query');
                await createUserSession(
                    { query: currentQuery },
                    {
                        ...callbacks,
                        onSessionCreated: (data) => {
                            console.log('Session Created:', data);
                            if (data && data.session_id) {
                                setCurrentSessionId(data.session_id);
                                // updateUrlWithSessionId(data.session_id);
                            }
                        }
                    }
                );
            }
        } catch (error) {
            console.error('Error during query processing:', error);
            const catchErrorMessage = `Failed to process query: ${error.message || 'Unknown error'}`;
            updateComponentResponse(cUid, { type: "error", data: catchErrorMessage, timestamp: new Date().toISOString() });
            setGenerationState("default");
            setActiveUid("");
        }
    };


    // REMOVED: useEffect hook that called init_response based on activeUid
    // REMOVED: init_response function

    useEffect(() => {
        console.log('modelConfig', modelConfig)
    }, [modelConfig])

    useEffect(() => {
        if (generationState === 'default' && isLoaded) {
            localStorage.setItem("components", JSON.stringify(components));
            // Session ID is now only maintained in React state
        }
    }, [components, isLoaded, generationState]);

    // REMOVED: useEffect hook that called generateResponseStream
    // REMOVED: generateResponseStream function

    // This useEffect is for debugging and can be kept or removed
    useEffect(() => {
        console.log("generation state", generationState);
    }, [generationState])

    // Function to update URL with session ID without triggering a page reload
    const updateUrlWithSessionId = (sessionId) => {
        if (sessionId) {
            const url = new URL(window.location);
            url.searchParams.set('session_id', sessionId);
            window.history.pushState({}, '', url);
            console.log('Updated URL with session ID:', sessionId);
        }
    };

    const repeatedChatElements =
        components.map((item, index) => {
            return (<>
                {/* main chat element   */}
                <div className={`chat-elem ${index === (components.length - 1) ? "last-load" : ""}`}>
                    {/* query container */}
                    <div
                    // variants={sideIconVariant}
                    // initial="initial"
                    // animate="animate"
                    >
                        <div className="query-container">
                            {/* <div className={`query-space ${isPanelOpen ? "close" : ""}`}>
                        </div> */}
                            <div className={`chat-query ${isContentLoading ? "loading" : ""}`}>
                                {/* upper buttons */}

                                {/* <div className={`btn-container ${isContentLoading ? "loading" : ""}`}>
                                    <button className="del-btn"
                                        onClick={() => { delHandler(index) }}><FaTrash className="chat-delete-icon" /></button>
                                </div> */}
                                <div className="query-text-div-container">
                                    <p>
                                        {item.content[0].query}
                                    </p>
                                </div>
                            </div>
                            <div className="query-profile">
                                <FaUserCircle className="chat-user-icon" alt="user profile" /> {/* Replaced img with FaUserCircle icon */}
                            </div>
                        </div>
                    </div>
                    <div
                    // variants={sideIconVariant}
                    // initial="initial"
                    // animate="animate"
                    >
                        <div className="response-container">

                            <div className="response-profile">
                                {index === (components.length - 1) ?
                                    <img src="/images/logo_rendered_main.png"
                                        alt="logo"
                                        className={`ai-profile-image ${generationState === "init" || generationState === "generating"
                                            ? "loading" : ""}`}
                                    />
                                    :
                                    <img src="/images/logo_rendered_main.png"
                                        alt="logo"
                                        className={`ai-profile-image`}

                                    />
                                }
                            </div>
                            <div className="indiv-response-div">
                                {index === (components.length - 1) ?
                                    <div className={`chat-response ${generationState === "init" ? "loading" : ""}`}>
                                        {/* <CustomMarkdownRenderer content={item.content[0].response} /> */}
                                        <FlexRenderer items={item.content[0].response} />
                                    </div>
                                    :
                                    <div className={`chat-response`}>
                                        {/* <CustomMarkdownRenderer content={item.content[0].response} /> */}
                                        <FlexRenderer items={item.content[0].response} />
                                    </div>
                                }
                                {index === (components.length - 1) && generationState !== "default" ?
                                    ""
                                    :
                                    <div
                                    // variants={sideIconVariant}
                                    // initial="initial"
                                    // animate="animate"
                                    >
                                        {/* <div className={`res-btn-container ${isContentLoading ? "loading" : ""}`}>
                                            <button className="copy-btn" onClick={() => { handleClipboard(index) }}
                                            // onMouseEnter={(e) => { handleShowTag(e, "Copy Response") }}
                                            // onMouseLeave={handleHideTag}
                                            ><img src="/icons/copy to clipboard.png" /></button>
                                        </div> */}
                                    </div>
                                }
                            </div>
                        </div>
                    </div>
                </div>
            </>)
        })

    // Simple UI to display messages
    // TODO: Enhance this UI significantly
    return (
        <>
            <div className="inter-main-container">
                <div className={`inter-inner-container`}
                    style={{ fontSize: `${fontSize}px` }}
                >
                    <div
                        ref={chatContainerRef}
                        className={`chat-container`} id="blur-div">
                        {components.length === 0 ? (
                            <motion.div className="suggestion-queries-container"
                                initial={{ opacity: 0, y: 20 }}
                                animate={{ opacity: 1, y: 0 }}
                                transition={{ duration: 0.5 }}
                            >
                                {username && (
                                    <motion.div className="welcome-message"
                                        initial={{ opacity: 0, y: -20 }}
                                        animate={{ opacity: 1, y: 0 }}
                                        transition={{ duration: 0.5, delay: 0.2 }}
                                    >
                                        Hello <span className="username-gradient">{username}</span>, Welcome to Valency!
                                    </motion.div>
                                )}
                                <motion.div className="suggestion-queries-list-wrapper"
                                    initial={{ opacity: 0 }}
                                    animate={{ opacity: 1 }}
                                    transition={{ duration: 0.5, delay: 0.4 }}
                                >
                                    <div className="suggestion-queries-list">
                                        {suggestionQueries.map((suggestion, index) => (
                                            <button
                                                key={index}
                                                className="suggestion-query-button"
                                                onClick={() => {
                                                    setQuery(suggestion);
                                                    // Optionally, you can also immediately send the query
                                                    // generateResponse();
                                                }}
                                            >
                                                {suggestion}
                                            </button>
                                        ))}
                                    </div>
                                </motion.div>
                            </motion.div>
                        ) : (
                            repeatedChatElements
                        )}
                        {/* scroll to bottom button */}
                        {/* <div ref={bottomRef} id="bottomRef"></div> */}
                    </div>

                    <div className="bottom-container" style={{ marginBottom: `${components.length === 0 ? "10vh" : "5vh"}` }}>
                        <div className="query-gen-container">
                            <div className={`query-text-container ${isContentLoading ? "loading" : ""}`}>
                                <textarea placeholder="Enter your query here"
                                    style={{ fontSize: `${fontSize}px` }}
                                    disabled={generationState === "generating"}
                                    value={query}
                                    onChange={(e) => setQuery(e.target.value)}
                                    className="query-text-object"
                                    onKeyDown={handleKeyDown}
                                />
                                <div className={`query-button-container ${isContentLoading ? "loading" : ""}`}>
                                    <button className="chat-query-btn" onClick={() => generateResponse()} disabled={generationState === "generating"}>
                                        <IoSend className="" />
                                    </button>
                                </div>
                            </div>
                        </div>
                        {currentSessionId && (
                            <div className="session-info" style={{marginTop: "20px", color: "gray"}}>
                                <span>Current Session ID: </span>
                                <span className="session-id">{currentSessionId}</span>
                            </div>
                        )}
                    </div>
                </div>

                {/* <div
            ref={scrollBottomRef}
            className={`bottom-scroller ${isDark ? "dark" : ""}`} // KEPT isDark here if it's a separate component not affected by html class
            onClick={scrollToBottom}
        >
            <img src="/icons/generate_dark.svg" alt="logo" />
        </div> */}

                {
                    showTag ?
                        <div className={`info-tag`}
                            style={{ transform: `translate(${tagPosition.x}px, ${tagPosition.y}px)` }}
                        >
                            <h1>{tagText}</h1>
                        </div>
                        :
                        ""
                }
            </div >
            {/* </GlassyContainer> */}
        </>
    );
}
export default ChatInterface;