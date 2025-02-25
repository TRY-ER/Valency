import axios from 'axios';
import { EventSourcePolyfill } from 'event-source-polyfill';

const call_endpoint_async = async (data, payload) => {
    if (data.method === 'post') {
        return await axios.post(data.url, payload)
    }
    else if(data.method === 'get') {
        return await axios.get(data.url)
    }
    // else {
    //     throw new Error('Invalid method')
    // }
}

const call_eventsource = (data, url_payload) =>{
    if (data.method === 'get@sse') {
        return new EventSourcePolyfill(data.url + "/" + url_payload)
    }
}

export {
    call_endpoint_async,
    call_eventsource
}