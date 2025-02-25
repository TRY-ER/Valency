import React from "react"
import ReactDOMServer from "react-dom/server";

const END_TOKENS = ["\nUser:", "<|endoftext|>", "</s>", "|e|"]

const getHTMLString = (element) => {
    return ReactDOMServer.renderToString(element)
}

const responseParser = (
    parsable,
    setIsStreaming,
    setIsCompleted,
    updateEnResponse,
    updateEnQuery,
    setChatStream,
    updateComponentResponse,
    cUid,
    resolve,
    isDark,
    SpinningLoader,
    scrollToBottom
) => {
    if (END_TOKENS.includes(parsable)) {
        setIsStreaming(false);
        resolve();
        console.log("promise is complete")
        setIsCompleted(true);
    }
    else if (parsable.startsWith("[final-res]")) {
        const final_res = parsable.split('[final-res]')[1];
        updateEnResponse(cUid, final_res);
    }
    else if (parsable.startsWith("[final-query]")) {
        const final_query = parsable.split('[final-query]')[1];
        updateEnQuery(cUid, final_query);
    }
    else if (parsable.startsWith("<|tag|>")) {
        var tag_value = parsable.split('<|tag|>')[1];
        // replace "\" from tag_value
        tag_value = tag_value.replace(/\\/g, "");
        // read json object from final_query
        const tag_json = JSON.parse(tag_value);
        const loadingHTMLString = getHTMLString(<SpinningLoader isDark={isDark} />);
        if (tag_json.status === "loading") {
            setChatStream((prev) => {
                const newChatStream = prev + `<div className="int-tag pending ${isDark ? "dark" : ""}" ><h1>${tag_json.title}</h1>${loadingHTMLString}</div>`;
                updateComponentResponse(cUid, newChatStream);
                return newChatStream;
            });
        }
        else {
            setChatStream((prev) => {
                console.log("prev chat stream", prev)
                // now replacing the <p>Loading...</p> element from previous with ""
                const newChatStream = prev.replace(
                    `<div className="int-tag pending ${isDark ? "dark" : ""}" ><h1>${tag_json.title}</h1>${loadingHTMLString}</div>`,
                    `<div className="int-tag ${isDark ? "dark" : ""}" ><h1>${tag_json.title}</h1></div>`)

                console.log("mod chat stream", newChatStream)
                updateComponentResponse(cUid, newChatStream);
                return newChatStream;
            });
        }
    }
    else {
        setChatStream((prev) => {
            const newChatStream = prev + parsable;
            updateComponentResponse(cUid, newChatStream);
            scrollToBottom();
            return newChatStream;
        });
    }
}


const removeHtmlTags = (content) => {
    return content.replace(/<[^>]*>/g, '');
}

export { removeHtmlTags, responseParser } 