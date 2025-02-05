

function getHTMLStart(text) {
    const tags = [];
    let startIndex = 0;
    let match;

    while ((match = text.substring(startIndex).match(/<#HTML: start#>/))) {
        const tagIndex = startIndex + match.index;
        tags.push({ startIndex: tagIndex, endIndex: tagIndex + match[0].length });
        startIndex = tagIndex + match[0].length;
    }

    return tags;
}

function getHTMLEnd(text) {
    const tags = [];
    let startIndex = 0;
    let match;

    while ((match = text.substring(startIndex).match(/<#HTML: end#>/))) {
        const tagIndex = startIndex + match.index;
        tags.push({ startIndex: tagIndex, endIndex: tagIndex + match[0].length });
        startIndex = tagIndex + match[0].length;
    }

    return tags;
}

function getJSONStart(text) {
    const tags = [];
    let startIndex = 0;
    let match;

    while ((match = text.substring(startIndex).match(/<#JSON: start#>/))) {
        const tagIndex = startIndex + match.index;
        tags.push({ startIndex: tagIndex, endIndex: tagIndex + match[0].length });
        startIndex = tagIndex + match[0].length;
    }

    return tags;
}

function getJSONEnd(text) {
    const tags = [];
    let startIndex = 0;
    let match;

    while ((match = text.substring(startIndex).match(/<#JSON: end#>/))) {
        const tagIndex = startIndex + match.index;
        tags.push({ startIndex: tagIndex, endIndex: tagIndex + match[0].length });
        startIndex = tagIndex + match[0].length;
    }

    return tags;
}

export { getHTMLStart, getHTMLEnd, getJSONStart, getJSONEnd };
