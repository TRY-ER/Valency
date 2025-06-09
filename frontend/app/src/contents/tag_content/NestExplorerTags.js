import BasicProtViewer from "../../pages/Explorer/ProtExplorer/BasicProtViewer"
import MolStarViewer from "../../pages/Explorer/ProtExplorer/MolStarViewer";
import UniProtViewer from "../../pages/Explorer/ProtExplorer/UniProtViewer"; // Import the new component
import DocRenderer from "../doc_content/DocRenderer";
import UniProtSummaryViewer from "../../pages/Explorer/ProtExplorer/UniProtSummaryViewer";
import AlphafoldAnnotationsViewer from "../../pages/Explorer/ProtExplorer/AlphafoldAnnotationsViewer";
import RCSBPDBExplorer from "../../pages/Explorer/RCSBPDBExplorer/RCSBPDBExplorer";
import RCSBStructureSimilaritySearch from "../../pages/Explorer/RCSBPDBExplorer/RCSBStructureSimilaritySearch";
import UniProtExplorer from "../../pages/Explorer/UniProtExplorer/UniProtExplorer";

const ProtExploreTabContents = [
    {
        id: 1,
        title: 'Basic',
        link: '',
        description: `This tool uses SMILES string to generate molecules using BRICS algorithm. The algorithm generates hypothetical molecules based on the given input.
                      The input can be a SMILES string in a text file each represented in a new line. The output will be a set of SMILES that are generated 
                      based on the input. The tool has three sections, the first section shows a progress bar showing the completed steps for the job.
                      The second section contains the file upload setup with a panel to show first few files of file to ensure the right text file has been selected.
                      The third section is the output section that consists of a output panel, download button and the reset button.
                      the output panel shows the steps details as well as completed steps and the molecules generated. This tool is specifically useful when you have
                      a set of molecules and you want to generate more molecules based on the given set.`,
        component: <BasicProtViewer />, 
        docs: <DocRenderer filePath="/markdown_repo/BRICSSMILES.md"/>,
    },
    {
        id: 2,
        title: 'Advanced',
        link: 'adv',
        description: `This tool uses SMILES string to generate molecules using BRICS algorithm. The algorithm generates hypothetical molecules based on the given input.
                      The input can be a SMILES string in a text file each represented in a new line. The output will be a set of SMILES that are generated 
                      based on the input. The tool has three sections, the first section shows a progress bar showing the completed steps for the job.
                      The second section contains the file upload setup with a panel to show first few files of file to ensure the right text file has been selected.
                      The third section is the output section that consists of a output panel, download button and the reset button.
                      the output panel shows the steps details as well as completed steps and the molecules generated. This tool is specifically useful when you have
                      a set of molecules and you want to generate more molecules based on the given set.`,
        component: <MolStarViewer />, 
        // docs: <DocRenderer filePath="/markdown_repo/BRICSSMILES.md"/>,
    }
]

const UniProtExploreTabContents = [
    {
        id: 1, // New ID for the UniProt tab
        title: 'Basic',
        link: '', // New link for the UniProt tab
        description: `Explore protein structures and information from AlphaFold using a UniProt accession key. Fetches data from the AlphaFold DB API and displays the 3D model, sequence, and other relevant details.`,
        component: <UniProtViewer />,
        docs: null 
    },
    {
        id: 2, // New ID for the UniProt tab
        title: 'Summary',
        link: 'summary', // New link for the UniProt tab
        description: `Explore protein structures and information from AlphaFold using a UniProt accession key. Fetches data from the AlphaFold DB API and displays the 3D model, sequence, and other relevant details.`,
        component: <UniProtSummaryViewer />, 
        docs: null 
    },
    {
        id: 3, // New ID for the AlphaFold Annotations tab
        title: 'Annotations',
        link: 'annotations', // New link for the AlphaFold Annotations tab
        description: `View AlphaFold annotations for a specific UniProt accession. Fetches annotation data from the AlphaFold API and displays it in an interactive visualization with detailed annotation tracks.`,
        component: <AlphafoldAnnotationsViewer />, 
        docs: <DocRenderer filePath="/markdown_repo/AlphaFoldAnnotationsViewer.md"/>
    }
]

const RCSBExplorerTabContent= [
    {
        id: 1, // New ID for the UniProt tab
        title: 'Getter',
        link: '', // New link for the UniProt tab
        description: `RCSB PDB Explorer is a tool that allows users to search and retrieve protein structures from the RCSB Protein Data Bank (PDB) using various identifiers such as PDB ID, UniProt ID, or Gene Name. It provides a user-friendly interface for exploring protein structures and their associated data.`,
        component: <RCSBPDBExplorer />,
        docs: null 
    },
    {
        id: 2, // New ID for the UniProt tab
        title: 'Similarity Search',
        link: 'ssearch', // New link for the UniProt tab
        description: `Similarity search in RCSB PDB Explorer allows users to find protein structures that are similar to a given structure based on sequence or structural similarity. This feature helps in identifying homologous proteins and understanding structural relationships.`,
        component: <RCSBStructureSimilaritySearch />,
        docs: null 
    },
]

const UniprotTabContent= [
    {
        id: 1, // New ID for the UniProt tab
        title: 'Search',
        link: '', // New link for the UniProt tab
        description: `Uniprot Explorer is a tool that allows users to search and retrieve protein information from the Uniprot database using various identifiers such as UniProt ID, Gene Name, or Protein Name. It provides a user-friendly interface for exploring protein sequences, functions, and annotations. with query and UniProt ID.`,
        component: <UniProtExplorer />,
        docs: <DocRenderer filePath="/markdown_repo/UniProtExplorer.md"/>
    },
]


export { ProtExploreTabContents, UniProtExploreTabContents, RCSBExplorerTabContent, UniprotTabContent};