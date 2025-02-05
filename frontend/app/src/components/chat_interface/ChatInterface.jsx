import React, { useEffect, useState, useRef, use } from "react";
import { useNavigate } from "react-router-dom";
import axios from "axios";
import { v4 as uuidv4 } from "uuid";
import { EventSourcePolyfill } from "event-source-polyfill";
import clipboard from "clipboard-copy";
import { removeHtmlTags, responseParser } from "../../components/parsers/Basic";
import { motion } from "framer-motion";
import "./ChatInterface.css";
import { IoSend } from "react-icons/io5";
import { chat_endpoint } from "../../endpoints/endpoints";
import { baseURL } from "../../endpoints/base";
import { call_endpoint_async, call_eventsource } from "../../endpoints/caller";

import { MODEL_OPTIONS_TYPES } from "./extendedVariables";
import masterParser from "../parsers/MasterParser";
import SpinningLoader from "../UI/SpinLoader/SpinLoader";
import CustomToggler from "../UI/CustomToggler/CustomToggler";
import CustomDropdown from "../UI/CustomDropdown/CustomDropdown";
import GlassyContainer from "../glassy_container/gc";

const ChatInterface = () => {
    const [isPanelOpen, setIsPanelOpen] = useState(false);
    const [profileExpanded, setProfileExpanded] = useState(false);
    const [chatStream, setChatStream] = useState("");
    const [chatSummary, setChatSummary] = useState("");
    const [isStreaming, setIsStreaming] = useState(false);
    const [isCompleted, setIsCompleted] = useState(false);
    const [generationState, setGenerationState] = useState("default"); // init, complete, default, generating (four value can be given)
    const [isContentLoading, setIsContentLoading] = useState(true);
    const [query, setQuery] = useState("");
    const navigate = useNavigate();
    const bottomRef = useRef(null);
    const [isRecording, setIsRecording] = useState(false);
    const [activeSpeaker, setActiveSpeaker] = useState(null);
    const [sideContent, setSideContent] = useState(null);
    const [serverStreamId, setServerStreamId] = useState(null);
    const [isLoaded, setIsLoaed] = useState(false);

    //ui changes 
    const [isDark, setIsDark] = useState(false);
    const [fontSize, setFontSize] = useState("14");


    //info-tag setup
    const [showTag, setShowTag] = useState(false);
    const [tagText, setTagText] = useState("test");
    const [tagPosition, setTagPosition] = useState({ x: 0, y: 0 });

    // input config for model
    const [modelConfig, setModelConfig] = useState({
        "model_name": "Gemini 2 Flash",
        "model_code": "gemini-2.0-flash-exp",
    })

    // change font size 
    const changeFontSize = (e) => {
        setFontSize(e.target.value);
    }

    useEffect(() =>{
        if (isLoaded){
            localStorage.setItem("fontsize", fontSize)
        }
    }, [fontSize])

    const sideIconVariant = {
        initial: {
            opacity: 0,
            y: "2vh", // Icon starts outside the viewport at the bottom
        },
        animate: {
            opacity: 1,
            y: 0, // Icon slides up to its final position
            transition: {
                duration: 0.5,
                ease: "easeOut",
            },
        },
    };

    const sidePanelVariant = {
        collapsed: {
            width: "0%",
        },
        expanded: {
            width: "100%",
            transition: {
                duration: 0.5,
                ease: "easeOut",
            },
        },
    };

    // useEffect(() => {
    //     console.log("isDark", isDark)
    // }, [isDark])

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

    //advance setting show hide
    const [showAdvanceSetting, setShowAdvanceSetting] = useState(false);


    // component states 
    const [components, setComponents] = useState([]);
    const [bottomScroller, setBottomScroller] = useState(false);
    const [activeUid, setActiveUid] = useState("");


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

    const updateEnResponse = (uid, newResponse) => {
        setComponents((prev) => {
            const new_value = prev.map((component) =>
                component.uid === uid ?
                    {
                        ...component,
                        content: [{
                            ...component.content[0],
                            response: {
                                ...component.content[0]["response"],
                                en_response: newResponse
                            }
                        }]
                    } : component
            )
            return new_value;
        }
        );
    }

    const updateEnQuery = (uid, newResponse) => {
        setComponents((prev) => {
            const new_value = prev.map((component) =>
                component.uid === uid ?
                    {
                        ...component,
                        content: [{
                            ...component.content[0], en_query: newResponse
                        }]
                    } : component
            )
            return new_value;
        }
        );
    }

    const updateConfig = async () => {
        localStorage.setItem("modelConfig", JSON.stringify(modelConfig));
    }

    const handlePanelToggle = () => {
        setIsPanelOpen(!isPanelOpen);
        setShowTag(false);
    }

    useEffect(() => {
        setQuery("");
    }, [modelConfig])

    //show config content
    const showConfig = () => {
        setIsPanelOpen(true);
        setSideContent("show_config");
    }

    // show profile content
    const showProfile = () => {
        setIsPanelOpen(true);
        setSideContent("show_profile");
    }

    const changeQuery = (text_value) => {
        setQuery(query => {
            return query + " " + text_value
        }
        )
    }

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


    const handleLike = async (index) => {
        if (components[index].content[0].response.like === 0 || components[index].content[0].response.like === 1) {
            // await axios.post(API_BASE + "/chat/handle_like", {
            //     "user_id": auth?.userName,
            //     "c_index": index,
            //     "like_value": 2
            // }, privateConfig).then((res) => {
            //     console.log("like_res", res)
            setComponents((prev) => {
                const new_value = prev.map((component, id) =>
                    id === index ?
                        {
                            ...component,
                            content: [{
                                ...component.content[0], response: {
                                    ...component.content[0].response,
                                    like: 2,
                                }
                            }]
                        } : component
                )
                return new_value;
            });
            // }).catch(err => {
            //     console.log("error in like data transfer", err)
            // })
        }
        else if (components[index].content[0].response.like === 2) {
            // await axios.post(API_BASE + "/chat/handle_like", {
            //     "user_id": auth?.userName,
            //     "c_index": index,
            //     "like_value": 0
            // }, privateConfig).then((res) => {
            //     console.log("like_res", res)
            setComponents((prev) => {
                const new_value = prev.map((component, id) =>
                    id === index ?
                        {
                            ...component,
                            content: [{
                                ...component.content[0], response: {
                                    ...component.content[0].response,
                                }
                            }]
                        } : component
                )
                return new_value;
            });
            // }).catch(err => {
            //     console.log("error in like data transfer", err)
            // })

        }
    }


    const handleDislike = async (index) => {
        if (components[index].content[0].response.like === 0 || components[index].content[0].response.like === 2) {
            // await axios.post(API_BASE + "/chat/handle_like", {
            //     "user_id": auth?.userName,
            //     "c_index": index,
            //     "like_value": 1
            // }, privateConfig).then((res) => {
            setComponents((prev) => {
                const new_value = prev.map((component, id) =>
                    id === index ?
                        {
                            ...component,
                            content: [{
                                ...component.content[0], response: {
                                    ...component.content[0].response,
                                    like: 1,
                                }
                            }]
                        } : component
                )
                return new_value;
            });
            // }).catch(err => {
            //     console.log("error in dislike data transfer", err)
            // });
        }

        else if (components[index].content[0].response.like === 1) {
            // await axios.post(API_BASE + "/chat/handle_like", {
            //     "user_id": auth?.userName,
            //     "c_index": index,
            //     "like_value": 0
            // }, privateConfig).then((res) => {
            setComponents((prev) => {
                const new_value = prev.map((component, id) =>
                    id === index ?
                        {
                            ...component,
                            content: [{
                                ...component.content[0], response: {
                                    ...component.content[0].response,
                                    like: 0,
                                }
                            }]
                        } : component
                )
                return new_value;
            });
            // }).catch(err => {
            //     console.log("error in dislike data transfer", err)
            // });
        }
    }

    const loadThread = () => {
        // console.log("load thread is triggered")
        setComponents(localStorage.getItem("components") ? JSON.parse(localStorage.getItem("components")) : []);
        setModelConfig(localStorage.getItem("modelConfig") ? JSON.parse(localStorage.getItem("modelConfig")) : modelConfig);
        setFontSize(localStorage.getItem("fontsize") ? localStorage.getItem("fontsize") : "14");
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
    }, [])

    // useEffect(() => {
    //     console.log("components >>", components)
    // }, [components])


    const updateComponentResponse = (uid, newResponse) => {
        setComponents((prev) => {
            // console.log('new chat stream >>', newResponse)
            // console.log("incoming uid >>", uid)
            const new_value = prev.map((component) =>
                component.uid === uid ?
                    {
                        ...component,
                        content: [{
                            ...component.content[0], response: {
                                ...component.content[0].response,
                                res_content: newResponse,
                                like: 0,
                            }
                        }]
                    } : component
            )
            // console.log('new value >>', newResponse)
            return new_value;
        }
        );
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

    const generateResponse = async () => {
        var requestable = "";
        requestable = query;
        if (requestable.length == 0) {
            alert("Enter some query first to generate response");
        }
        else {
            setChatStream("");
            let cUid = uuidv4();
            let new_comp = {
                "content": [{
                    "query": requestable,
                    "response": {
                        "res_content": "",
                        "like": 0,
                    },
                }],
                "uid": cUid,
            }
            setComponents([...components, new_comp])
            setActiveUid(cUid);
        }
    }


    useEffect(() => {
        // console.log("active", activeUid);
        // console.log("updated_component", components)

        const init_response = async () => {
            const response = await call_endpoint_async(chat_endpoint.init, {
                query: query,
                config: {
                    "model_name": modelConfig["model_code"]
                }
            })
            if (response.data.status === "success") {
                setServerStreamId(response.data.id);
                setGenerationState("init");
            }
        }
        if (activeUid !== "") {
            init_response();
        }
        scrollToBottom();
    }, [activeUid])

    useEffect(() => {
        console.log('modelConfig', modelConfig)
    }, [modelConfig])

    useEffect(() => {
        if (generationState === 'default' && isLoaded) {
            localStorage.setItem("components", JSON.stringify(components));
        }
    }, [components])

    useEffect(() => {
        if (generationState === "init") {
            if (serverStreamId) {
                generateResponseStream(serverStreamId).then(() => {
                    setGenerationState("complete");
                }).catch(err => {
                    console.log("error in generating streams", err)
                })
            }
        }
        else if (generationState === "complete") {
            localStorage.setItem("components", JSON.stringify(components));
            setGenerationState("default")
            // const save_response = async () => {
            //     var save_dict = {
            //         "user_id": auth?.userName,
            //         "component": {
            //             "content": [{
            //                 "query": components[components.length - 1]?.content[0].query,
            //                 "audio_query": null,
            //                 "en_query": components[components.length - 1]?.content[0].en_query,
            //                 "response": components[components.length - 1]?.content[0].response
            //             }],
            //             "location": [20, 3, 25.3],
            //             "uid": activeUid,
            //             "config": modelConfig
            //         },
            //         "summary": chatSummary,
            //         "component_length": components.length,
            //         "uid": activeUid
            //     }
            //     await axios.post(`${API_BASE}/chat/save_response`, save_dict, privateConfig).then((res) => {
            //         console.log("save response", res)
            //         setGenerationState("default")
            //     }, privateConfig).catch((err) => {
            //         console.log("error >>", err);
            //     })
            // }
            // save_response();
        }
        else if (generationState === "default") {
            setActiveUid("");
            setServerStreamId(null);
        }
    }
        , [generationState])

    // useEffect(() => {
    // if (auth.userName && auth?.localAccessToken) {
    //     setSideContent("show_profile");
    //     console.log(auth.user);
    //     console.log("auth is triggered")
    //     if (auth?.localAccessToken) {
    //         console.log("local access is triggered")
    //         setPrivateConfig(() => {
    //             console.log("access token", auth?.localAccessToken)
    //             return {
    //                 headers: {
    //                     "Content-Type": "application/json",
    //                     "Authorization": `Bearer ${auth?.localAccessToken}`,
    //                 },
    //             }
    //         });
    //     }
    //     loadThread(auth?.userName, auth?.localAccessToken)
    // }
    // }, [auth])

    const generateResponseStream = async (cUid) => {
        // console.log("gen res is called");
        return new Promise((resolve) => {
            if (!isStreaming) {
                const eventSource = call_eventsource(chat_endpoint.stream, cUid);
                eventSource.onopen = () => {
                    setIsStreaming(true);
                    setIsCompleted(false);
                    setQuery("");
                    setGenerationState("generating")
                }
                eventSource.onmessage = async (event) => {
                    // console.log("msg data >>", event.data);
                    if (event.data === "<|end|>") {
                        console.log("end token is received")
                        setIsStreaming(false);
                        eventSource.close();
                        resolve();
                    }
                    else {
                        setChatStream((prev) => {
                            const newChatStream = prev + event.data;
                            updateComponentResponse(activeUid, newChatStream);
                            scrollToBottom();
                            return newChatStream;
                        });
                    }
                    // responseParser(
                    //     event.data,
                    //     setIsStreaming,
                    //     setIsCompleted,
                    //     updateEnResponse,
                    //     updateEnQuery,
                    //     setChatStream,
                    //     updateComponentResponse,
                    //     cUid,
                    //     resolve,
                    //     isDark,
                    //     SpinningLoader,
                    //     scrollToBottom,
                    // )
                }

                eventSource.onerror = () => {
                    setIsStreaming(false);
                    eventSource.close();
                    resolve();
                }
                // setIsStreaming(true);
            }
            else {
                setIsStreaming(false);
                resolve();
            }
        })
    }

    useState(() => {
        console.log("generation state", generationState);
    }, [generationState])

    const profileExpandContent = (
        <>
            <div
            // variants={sideIconVariant}
            // initial="initial"
            // animate="animate"
            >
                <div className={`profile-expand-content ${isDark ? "dark" : ""} ${profileExpanded ? "open" : ""}`}>
                    <div className="profile-content">
                        <h1>Adjust Font Size</h1>
                        <div className={`slider-container ${isDark ? "dark" : ""} ${profileExpanded ? "open" : ""}`}>
                            <input type="range" min="12" max="24" value={fontSize} step="4" onChange={changeFontSize} className={`slider ${isDark ? "dark" : ""} ${profileExpanded ? "open" : ""}`} />
                            <div className={`slider-labels ${isDark ? "dark" : ""} ${profileExpanded ? "open" : ""}`}>
                                <span className={fontSize === "12" ? 'slider-active-label' : ''}>Small</span>
                                <span className={fontSize === "16" ? 'slider-active-label' : ''}>Medium</span>
                                <span className={fontSize === "20" ? 'slider-active-label' : ''}>Large</span>
                                <span className={fontSize === "24" ? 'slider-active-label' : ''}>Extra <br />Large</span>
                            </div>
                        </div>

                        <br />
                        <br />
                        <p style={{ fontSize: `${fontSize}px` }}>This text's font size will adjust.</p>
                        <br />
                        <br />
                        {/* <button onClick={logoutHandler}>Log Out</button> */}

                    </div>
                </div>
            </div>
        </>
    )


    const configContainer = (
        <>
            <div
            // variants={sideIconVariant}
            // initial="initial"
            // animate="animate"
            >
                <div className={`config-container ${isDark ? "dark" : ""}`}>
                    {/* <h1>Choose Language</h1>
                <CustomDropdown
                    options={INPUT_OPTIONS_TYPES}
                    configState={modelConfig}
                    configNames={["input_lang", "output_lang"]}
                    setConfig={setModelConfig}
                    isDark={isDark} /> */}
                    {/* <h1>Choose Language</h1>
            <CustomDropdown
                options={INPUT_OPTIONS_TYPES}
                configState={modelConfig}
                configName={"output_lang"}
                setConfig={setModelConfig} /> */}
                    <br />
                    <div>
                        <h1>Choose Model</h1>
                        <CustomDropdown
                            options={MODEL_OPTIONS_TYPES}
                            configState={modelConfig}
                            configNames={["model_name", "model_code"]}
                            setConfig={setModelConfig}
                            isDark={isDark}
                        />
                    </div>
                    <br />
                    <button className={`config-save-button ${isDark ? "dark" : ""}`}
                        onClick={updateConfig}>SAVE</button>
                </div>
            </div>
        </>)

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
                            <div className={`chat-query ${isContentLoading ? "loading" : ""} ${isDark ? "dark" : ""}`}>
                                {/* upper buttons */}

                                <div className={`btn-container ${isContentLoading ? "loading" : ""} ${isDark ? "dark" : ""}`}>
                                    <button className="del-btn"
                                        // onMouseEnter={(e) => { handleShowTag(e, "Delete Query and Response") }}
                                        // onMouseLeave={handleHideTag}
                                        onClick={() => { delHandler(index) }}><img src="icons/delete-icon-2.png" /></button>
                                </div>
                                <div>
                                    <p>
                                        {item.content[0].en_show ? item.content[0].en_query : item.content[0].query}
                                    </p>
                                </div>
                            </div>
                            <div className="query-profile">
                                <img src="icons/user-picture_dark.svg"
                                    alt="logo"
                                />
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
                                    <img src="/images/valency_logo_light_600x600.png"
                                        alt="logo"
                                        className={`ai-profile-image ${generationState === "init" || generationState === "generating"
                                            ? "loading" : ""} ${isDark ? "dark" : ""}`}
                                    />
                                    :
                                    <img src="/images/valency_logo_light_600x600.png"
                                        alt="logo"
                                        className={`ai-profile-image ${isDark ? "dark" : ""}`}

                                    />
                                }
                            </div>
                            <div className="indiv-response-div">
                                {index === (components.length - 1) ?
                                    <div className={`chat-response ${generationState === "init" ? "loading" : ""} ${isDark ? "dark" : ""}`}>
                                        <p style={{ fontSize: `${fontSize}px` }}>
                                            {masterParser(item.content[0].response.res_content, isDark)}
                                        </p>
                                    </div>
                                    :
                                    <div className={`chat-response ${isDark ? "dark" : ""}`}>
                                        <p>
                                            {item.content[0].en_show ? masterParser(item.content[0].response.en_response, isDark) : masterParser(item.content[0].response.res_content, isDark)}
                                        </p>
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
                                        <div className={`res-btn-container ${isContentLoading ? "loading" : ""} ${isDark ? "dark" : ""}`}>
                                            <button className="copy-btn" onClick={() => { handleClipboard(index) }}
                                            // onMouseEnter={(e) => { handleShowTag(e, "Copy Response") }}
                                            // onMouseLeave={handleHideTag}
                                            ><img src="icons/copy to clipboard.png" /></button>
                                            {item.content[0].response.like === 0 || item.content[0].response.like === 1 ?
                                                <button
                                                    className="like-btn"
                                                    onClick={() => { handleLike(index) }}><img src="icons/like-unfilled.png"
                                                    // onMouseEnter={(e) => { handleShowTag(e, "Like Response") }}
                                                    // onMouseLeave={handleHideTag}
                                                    /></button>
                                                :
                                                ""
                                            }
                                            {item.content[0].response.like === 2 ?
                                                <button className="like-btn"
                                                    // onMouseEnter={(e) => { handleShowTag(e, "Like Response") }}
                                                    // onMouseLeave={handleHideTag}
                                                    onClick={() => { handleLike(index) }}
                                                ><img src="icons/like-filled.png" /></button>
                                                :
                                                ""
                                            }

                                            {item.content[0].response.like === 0 || item.content[0].response.like === 2 ?
                                                <button className="dislike-btn" onClick={() => { handleDislike(index) }}
                                                // onMouseEnter={(e) => { handleShowTag(e, "Dislike Response") }}
                                                // onMouseLeave={handleHideTag}
                                                ><img src="icons/dislike-unfilled.png" /></button>
                                                :
                                                ""
                                            }
                                            {item.content[0].response.like === 1 ?
                                                <button className="dislike-btn" onClick={() => { handleDislike(index) }}
                                                // onMouseEnter={(e) => { handleShowTag(e, "Dislike Response") }}
                                                // onMouseLeave={handleHideTag}
                                                ><img src="icons/dislike-filled.png" /></button>
                                                :
                                                ""
                                            }
                                        </div>
                                    </div>
                                }
                            </div>
                        </div>
                    </div>
                </div>
            </>)
        })

    return (
        <>
            {/* <GlassyContainer> */}
            <div className="inter-main-container">
                <div className={`side-pane ${isPanelOpen ? "open" : ""} ${isDark ? "dark" : ""} custom-scrollbar`}>
                    <div className={`side-option ${isPanelOpen ? "open" : ""} custom-scrollbar`}>
                        <div className={`top-side-container ${isPanelOpen ? "open" : ""} ${isDark ? "dark" : ""} custom-scrollbar`}>
                            {isPanelOpen ? <>
                                <div
                                // variants={sideIconVariant}
                                // initial="initial"
                                // animate="animate"
                                >
                                    <img src="icons/anchor-left-arrow-dark.svg"
                                        alt="logo"
                                        // onMouseEnter={(e) => { handleShowTag(e, "Hide Pane") }}
                                        // onMouseLeave={handleHideTag}
                                        className={`side-collapser ${isDark ? "dark" : ""}`}
                                        onClick={handlePanelToggle}
                                    />
                                </div>
                            </> : ""}
                        </div>
                        <div className={`side-expand-container ${isPanelOpen ? "open" : ""} custom-scrollbar`}>
                            <div className="side-icon-container">
                                <div
                                // variants={sideIconVariant}
                                // initial="initial"
                                // animate="animate"
                                >
                                    <img src="icons/user-picture_dark.svg"
                                            alt="logo"
                                            id="side-profile-image"
                                            // onMouseEnter={(e) => { handleShowTag(e, "Profile", 12) }}
                                            // onMouseLeave={handleHideTag}
                                        onClick={showProfile}
                                        />
                                </div>
                                <div
                                // variants={sideIconVariant}
                                // initial="initial"
                                // animate="animate"
                                >
                                    <img src="icons/configuration-with-gear_dark.svg"
                                        alt="logo"
                                        className={`side-logo ${isDark ? "dark" : ""}`}
                                        onClick={showConfig}
                                    // onMouseEnter={(e) => { handleShowTag(e, "Config", 12) }}
                                    // onMouseLeave={handleHideTag}
                                    />
                                </div>
                            </div>
                            <div className={`side-content-container ${isPanelOpen ? "open" : ""} custom-scrollbar`}>
                                <div className={`side-panel ${isPanelOpen ? "open" : ""} custom-scrollbar`}>
                                    <div className="panel-content">
                                        {isPanelOpen && sideContent === "show_config" ? configContainer : ""}
                                        {isPanelOpen && sideContent === "show_profile" ? profileExpandContent : ""}
                                    </div>
                                </div>
                            </div>
                        </div>
                    </div>
                </div>
                <div className={`inter-inner-container ${isPanelOpen ? "panel-open" : ""} ${profileExpanded ? "profile-open" : ""} ${isDark ? "dark" : ""}`}
                    style={{ fontSize: `${fontSize}px` }}
                >

                    <div
                        ref={chatContainerRef}
                        className={`chat-container ${isPanelOpen ? "panel-open" : ""} ${isDark ? "dark" : ""}`} id="blur-div">
                        {repeatedChatElements}
                        {/* <div ref={bottomRef} id="bottomRef"></div> */}
                    </div>
                    <div
                    // variants={sideIconVariant}
                    // initial="initial"
                    // animate="animate"
                    >
                        <div className="bottom-container"


                        >
                            <div className="query-gen-container">
                                <div className={`query-text-container ${isContentLoading ? "loading" : ""} ${isDark ? 'dark' : ""}`}>
                                    <textarea placeholder="Enter your query here"
                                        style={{ fontSize: fontSize }}
                                        disabled={isRecording || generationState === "generating"}
                                        value={query}
                                        onChange={(e) => setQuery(e.target.value)}
                                        className="query-text-object"
                                        onKeyDown={handleKeyDown}
                                    />
                                </div>
                                <div className={`query-button-container ${isContentLoading ? "loading" : ""} ${isDark ? 'dark' : ""}`}>
                                    <button onClick={() => generateResponse()} disabled={isRecording && generationState === "generating"}>
                                        {/* <img src="/icon_second/send-filled.svg" alt="logo"
                                                onMouseEnter={(e) => { handleShowTag(e, "Generate Resonse") }}
                                                onMouseLeave={handleHideTag}
                                            /> */}
                                        <IoSend className="chat-query-btn" />
                                    </button>
                                </div>
                            </div>
                        </div>
                    </div>
                </div>

                {/* <div
                    ref={scrollBottomRef}
                    className={`bottom-scroller ${isDark ? "dark" : ""}`}
                    onClick={scrollToBottom}
                >
                    <img src="/icons/generate_dark.svg" alt="logo" />
                </div> */}

                {
                    showTag ?
                        <div className={`info-tag ${isDark ? "dark " : ""}`}
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
    )
}
export default ChatInterface;