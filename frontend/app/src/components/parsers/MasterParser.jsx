import parse from "html-react-parser";
import React, { useEffect, useState } from "react";
import { getHTMLStart, getHTMLEnd, getJSONStart, getJSONEnd } from "./utils.js";
import "./MasterParser.css";
import ConvertMapJSON from "./BaseJSONParser";

const HTMLContainer = (props) => {
    const content = props.content;
    return <>
        {parse(content)}
    </>
}


const masterParser = (rawResponse, isDark, isDarkboolean = false) => {
    var html_parsed_values = []
    var HTMLStarts = getHTMLStart(rawResponse);
    var JSONStarts = getJSONStart(rawResponse);
    var JSONEnds = getJSONEnd(rawResponse);
    var HTMLEnds = getHTMLEnd(rawResponse);
    if (HTMLStarts.length === 0) {
        if (JSONStarts.length > 0) {
            html_parsed_values.push(<HTMLContainer content={rawResponse.substring(0, JSONStarts[0]["startIndex"])} />);
            for (var i = 0; i < JSONStarts.length; i++) {
                html_parsed_values.push(ConvertMapJSON(rawResponse.substring(JSONStarts[i]["endIndex"], JSONEnds[i]["startIndex"])));
            }
        }
        else {
            html_parsed_values.push(<HTMLContainer content={rawResponse} />);
        }
    }
    else {
        html_parsed_values.push(<HTMLContainer content={rawResponse.substring(0, HTMLStarts[0]["startIndex"])} />);
    }
    if (HTMLStarts.length === HTMLEnds.length) {
        for (var i = 0; i < HTMLStarts.length; i++) {
            html_parsed_values.push(<HTMLContainer content={rawResponse.substring(HTMLStarts[i]["endIndex"], HTMLEnds[i]["startIndex"])} />);
        }
    }
    else if (HTMLStarts.length > HTMLEnds.length) {
        for (var i = 0; i < HTMLEnds.length; i++) {
            html_parsed_values.push(<HTMLContainer content={rawResponse.substring(HTMLStarts[i]["endIndex"], HTMLEnds[i]["startIndex"])} />);
        }
        html_parsed_values.push(<HTMLContainer content={rawResponse.substring(HTMLStarts[i]["endIndex"])} />);
    }
    if (JSONStarts.length > 0 && HTMLStarts.length !== 0) {
        for (var i = 0; i < JSONStarts.length; i++) {
            html_parsed_values.push(ConvertMapJSON(rawResponse.substring(JSONStarts[i]["endIndex"], JSONEnds[i]["startIndex"])));
        }
    }
    return <>
        <div className={`parsable-div ${isDark ? "dark" : ""}`}>
            {html_parsed_values.map((value, index) => {
                return <>{value}</>
            })}
        </div>
    </>
        ;
}


export default masterParser;    