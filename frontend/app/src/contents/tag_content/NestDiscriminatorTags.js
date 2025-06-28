import SSContent from "../../pages/Discriminator/SimilaritySearch/SSContent";
import DocRenderer from "../doc_content/DocRenderer";
import ChemBLGetter from "../../pages/Explorer/ChemBLExplorer/ChemBLGetter";
import ChemBLSimilarityGetter from "../../pages/Explorer/ChemBLExplorer/ChemBLSimilarityGetter";
import ApprovedDrugsViewer from "../../pages/Explorer/ChemBLExplorer/ApprovedDrugsViewer";
import ChemBLActivityFetcher from "../../pages/Explorer/ChemBLExplorer/ChemBLActivityFetcher";
import ChemBLUtilities from "../../pages/Explorer/ChemBLExplorer/ChemBLUtilities";
import PubChemGetter from "../../pages/Explorer/PubChemExplorer/PubChemGetter";
import PubChemSimilarityGetter from "../../pages/Explorer/PubChemExplorer/PubChemSimilarityGetter";
import PubChemUtilities from "../../pages/Explorer/PubChemExplorer/PubChemUtilities";

const SimilaritySearchTabContents = [
    {
        id: 1,
        title: 'SMILES',
        link: '',
        description: `This tool uses SMILES string to retrieve similar molecules from the chemical space.
                      This tool takes the SMILES string with a the number of candidates to retrieve and does
                      the operation and returns a section that shows the image of the molecules with it's SMILES and 
                      the similarity distance from the query SMILES. The results can be downloadable in a text file
                      as well as can be copied to the clipboard.`,
        component: <SSContent key="smiles" inputType={"MOL"} />, 
        docs: <DocRenderer filePath="/markdown_repo/MoleculeSimilaritySearch.md"/>,
        includeDocs: true
    },
    {
        id: 2,
        title: 'PSMILES',
        link: 'psmiles',
        description: `This tool uses PSMILES string to retrieve similar polymers from the polymer space.
                      This tool takes the PSMILES string with a the number of candidates to retrieve and does
                      the operation and returns a section that shows the image of the polymers with it's PSMILES and 
                      the similarity distance from the query PSMILES.The results can be downloadable in a text file
                      as well as can be copied to the clipboard.`,
        component:  <SSContent key="psmiles" inputType={"POLY"} />,
        docs: <DocRenderer filePath="/markdown_repo/PolymerSimilaritySearch.md"/>,
        includeDocs: true
    },
    {
        id: 3,
        title: 'Protein',
        link: 'protein',
        description: `This tool uses PDB ID string to retrieve similar proteins from the PDB space.
                      This tool takes the PDB ID string with a the number of candidates to retrieve and does
                      the operation and returns a section that shows the image of the proteins with it's PDB Id 
                      and similarity score. The results can be downloadable in a text file as well as can be 
                      copied to the clipboard.`,
        component:  <SSContent key="pdb" inputType={"PROT"} />,
        docs: <DocRenderer filePath="/markdown_repo/ProteinSimilaritySearch.md"/>,
        includeDocs: true
    }
]; 

const ChemBLTabContents = [
    {
        id: 1,
        title: 'Getter',
        link: '',
        description: `This tool uses ChemBL database to retrieve molecule information by ChEMBL ID or preferred name.
                      It allows you to search for molecules using either their unique ChEMBL identifier or their 
                      common/preferred name and displays comprehensive molecular data.`,
        component: <ChemBLGetter />,
        docs: <DocRenderer filePath="/markdown_repo/ChEMBLGetter.md"/>,
        includeDocs: true
    },
    {
        id: 2,
        title: 'Similarity Search',
        link: 'similarity-search',
        description: `This tool uses ChemBL database to conduct a similarity search for molecules.`,
        component: <ChemBLSimilarityGetter />,
        docs: <DocRenderer filePath="/markdown_repo/ChEMBLSimilarityGetter.md"/>,
        includeDocs: true
    },
    {
        id: 3,
        title: 'Activity Fetcher',
        link: 'activity-fetcher',
        description: `This tool fetches bioactivity data for biological targets and molecules using ChEMBL database.
                      You can search for target activities using target ChEMBL ID or molecule activities using molecule ChEMBL ID.`,
        component: <ChemBLActivityFetcher />,
        docs: <DocRenderer filePath="/markdown_repo/ChEMBLActivityFetcher.md"/>,
        includeDocs: true
    },
    // {
    //     id: 4,
    //     title: 'Approved Drugs',
    //     link: 'app-drugs',
    //     description: `This tool uses ChemBL database to retrieve information about approved drugs.`,
    //     component: <ApprovedDrugsViewer />,
    //     docs: null,
    // },
    {
        id: 4,
        title: 'Utilities',
        link: 'utils',
        description: `This tool provides various utilities for working with ChemBL data including target searches by gene name, SMILES to CTAB conversion, molecular descriptor calculations, structural alerts analysis, molecule standardization, and parent molecule extraction.`,
        component: <ChemBLUtilities />,
        docs: <DocRenderer filePath="/markdown_repo/ChEMBLUtilities.md"/>,
        includeDocs: true
    },
];

const PubChemTabContents = [
    {
        id: 1,
        title: 'Getter',
        link: '',
        description: `This tool uses PubChem database to retrieve compound information by CID, name, SMILES, or InChI Key.
                      It allows you to search for compounds using various identifiers and displays comprehensive molecular data
                      including structure, properties, and synonyms.`,
        component: <PubChemGetter />,
        docs: <DocRenderer filePath="/markdown_repo/PubChemGetter.md"/>,
        includeDocs: true
    },
    {
        id: 2,
        title: 'Similarity Search',
        link: 'similarity-search',
        description: `This tool uses PubChem database to conduct similarity searches for compounds based on CID or SMILES.
                      Find structurally similar compounds with adjustable similarity thresholds and comprehensive results display.`,
        component: <PubChemSimilarityGetter />,
        docs: <DocRenderer filePath="/markdown_repo/PubChemSimilarityGetter.md"/>,
        includeDocs: true,
    },
    {
        id: 3,
        title: 'Utilities',
        link: 'utils',
        description: `This tool provides various utilities for working with PubChem data including compound search by various identifiers,
                      substructure search, mass-based search, cross-reference lookup, synonym retrieval, and property calculation.`,
        component: <PubChemUtilities />,
        docs: <DocRenderer filePath="/markdown_repo/PubChemUtilities.md"/>,
        includeDocs: true,
    },
];

export {
    SimilaritySearchTabContents,
    ChemBLTabContents,
    PubChemTabContents
}