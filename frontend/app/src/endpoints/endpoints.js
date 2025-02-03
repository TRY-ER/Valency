import { baseURL } from './base';

const endpoints = {
    validate: {
        method: 'post',
        url: `${baseURL}/validate`,
        params: "type: for validation type. Should be one of [MOL, PROT, POLY], value: string representation of the molecule, protein, or polymer", 
        description: "This endpoint is used to validate the molecules, protiens, and polymers from their string representation.", 
        return_types: "status: string, valid: boolean"
    },
    visualize: {
        method: 'post',
        url: `${baseURL}/visualize`,
        params: "type: for visualization type. Should be one of [MOL, PROT, POLY], value: string representation of the molecule, protein, or polymer", 
        description: "This endpoint is used to visualize the molecules, protiens, and polymers from their string representation.", 
        return_types: "status: string, image_data: base64 encoded image data"
    },
    inform: {
        method: 'post',
        url: `${baseURL}/inform`,
        params: "type: for information type. Should be one of [MOL, PROT, POLY], value: string representation of the molecule, protein, or polymer", 
        description: "This endpoint is used to get information of the molecules, protiens, and polymers from their string representation.", 
        return_types: "status: string, info: dict of informations"
    }, 
} 


const generate_endpoints = {
    brics: {
        smiles: {
            method: 'post',
            url: `${baseURL}/generate/brics/smiles/upload`,
            params: "file: file containing the SMILES strings",
            description: "This endpoint is used to upload text files with SMILES strings",
            return_types: "status: string, file_id: uuid of the uploaded file, line_count: number of lines in the file", 
        },
        psmiles: {
            method: 'post',
            url: `${baseURL}/generate/brics/psmiles/upload`,
            params: "file: file containing the SMILES strings",
            description: "This endpoint is used to upload text files with SMILES strings",
            return_types: "status: string, file_id: uuid of the uploaded file, line_count: number of lines in the file", 
        },
        smiles_stream: {
            method: 'get@sse',
            url: `${baseURL}/generate/brics/smiles/stream`,
            params: "uuid: uuid of the uploaded file (passed in the url)",
            description: "This endpoint is used to stream the steps content and generated SMILES string files", 
            return_types: "server side event streaming", 
        },
        psmiles_stream: {
            method: 'get@sse',
            url: `${baseURL}/generate/brics/psmiles/stream`,
            params: "uuid: uuid of the uploaded file (passed in the url)",
            description: "This endpoint is used to stream the steps content and generated SMILES string files", 
            return_types: "server side event streaming", 
        },
    },
    lstm: {
        set: {
            method: 'post',
            url: `${baseURL}/generate/lstm/set`,
            params: "input_type: string, input_type should be one of [psmiles, wdg], num_gen: number of candidates to generate",
            description: "This endpoint is used to set the input type and number of candidates to generate",
            return_types: "status: string, id: uuid of the set up request"
        },
        stream: {
            method: 'get@sse',
            url: `${baseURL}/generate/lstm/stream`,
            params: "gen_id: uuid of the set up request (passed in the url)",
            description: "This endpoint is used to stream the generated content",
            return_types: "server side event streaming"
        }
    } 
}

const discriminate_endpoints = {
    ssearch: {
        method: 'post',
        url: `${baseURL}/discriminator/ssearch/local`,
        params: "input_type: string, input_type should be one of ['MOL', 'POLY', 'PROT'], data: query to do the similarity search with, k: number of candidates to generate",    
        description: "This endpoint is used to conduct structural similarity search with the query string",
        return_types: "status: string, results: list of results, message: string (retured if status is failed)"
    }
}

export {
    endpoints,
    generate_endpoints,
    discriminate_endpoints
}