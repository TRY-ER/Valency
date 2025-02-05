import { jsonTypeMap } from "./JSONMapper";


const ConvertMapJSON = (jsonObj) => {
    const json = JSON.parse(jsonObj);
    if ("type" in json) {
        if (json["type"] in jsonTypeMap) {
            return jsonTypeMap[json["type"]](json);
        }
    }
}

export default ConvertMapJSON;