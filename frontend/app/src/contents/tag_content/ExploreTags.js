import DocRenderer from "../doc_content/DocRenderer";
import MolEComponent from "../../pages/Explorer/MolExplorer/MolExplorer";
import ProtEComponent from "../../pages/Explorer/ProtExplorer/ProtExplorer";
import PolyEComponent from "../../pages/Explorer/PolyExplorer/PolyExplorer";
import RCSBComponent from "../../pages/Explorer/ProtExplorer/RCSBExplorer";
import UniProtExplorer from "../../pages/Explorer/UniProtExplorer/UniProtExplorer";
import { ProtExploreTabContents, UniProtExploreTabContents, RCSBExplorerTabContent,
    UniprotTabContent
 } from "./NestExplorerTags";
import { sub } from "framer-motion/client";

const ExplorerTabContent = [
    {
        id: 1,
        title: 'Molecule Explorer',
        link: '', // Keep this as empty for the default/index route under /explorer
        description: `The molecule explorer tool takes SMILES string of the molecule and provides interface for 
                     showcasing the 2D image of the chemical representation of the molecule and a information panel
                     for showcasing other properties. The information panel shows details as Molecular Formula,
                     Molecular Weight, Heavy Atoms Count, H Bond Doner Count, H Bond Acceptor Count, Rotatale Bonds Count,
                     Topological Polar Surface Area (TPSA), and Number of rings`,
        component: <MolEComponent />,
        docs: <DocRenderer filePath="/markdown_repo/MoleculeExplorer.md"/>,
    },
    {
        id: 2,
        title: 'Protein Explorer',
        link: 'proe', // This is correct for /explorer/proe
        description: `The protein explorer tool takes the Protein Data Bank ID (PDB ID) to show it's details and 3D structure
                      in two panels. The 3D structure is rotatable and scalable for exploration purposes. The information panel 
                      whows it's several details consisting Pdb Id, Title,Authors,Journal,Year,Volume,Pages,Doi,Pubmed Id,
                     Experiment Method,Molecular Weight (kDa),Deposited Model Count,Polymer entity count,Polymer monomer count,
                     Structural Features,Release Date,Resolution`,
        component: <ProtEComponent tabContent={ProtExploreTabContents} basePath="explorer/proe"/>, // Corrected basePath
        // docs: <DocRenderer filePath="/markdown_repo/ProteinExplorer.md"/>,
        subElements: ProtExploreTabContents,
        includeDocs: true,
    },
    {
        id: 3,
        title: 'Alphafold',
        link: 'alphafold', // This is correct for /explorer/proe
        description: `The alfafold explorer tool takes the Uniprot ID to show it's details and 3D structure`,
        component: <ProtEComponent tabContent={UniProtExploreTabContents} basePath="explorer/alphafold"/>, // Corrected basePath
        // docs: <DocRenderer filePath="/markdown_repo/ProteinExplorer.md"/>,
        subElements: UniProtExploreTabContents,
        includeDocs: true,
    },
    {
        id: 4,
        title: 'RCSB',
        link: 'rcsb', // This is correct for /explorer/proe
        description: `The RCSB explorer tool enables various mode of search for the PDB ID, Uniprot ID, or Gene Name.`,
        component: <RCSBComponent tabContent={RCSBExplorerTabContent} basePath="explorer/rcsb"/>, // Corrected basePath
        // docs: <DocRenderer filePath="/markdown_repo/ProteinExplorer.md"/>,
        subElements: RCSBExplorerTabContent,
        includeDocs: true,
    },
    {
        id: 5,
        title: 'Uniprot',
        link: 'uniprot', // This is correct for /explorer/uniprot
        description: `The Uniprot explorer tool enables various mode of search for query and UniProt ID`,
        component: <UniProtExplorer />, // Use UniProtExplorer directly
        docs: <DocRenderer filePath="/markdown_repo/UniProtExplorer.md"/>,
    },
    {
        id: 6,
        title: 'Polymer Explorer',
        link: 'polye', // This is correct for /explorer/polye
        description: `The polymer explorer tool takes PSMILES of the polymer and gives visual representation in a sidepanel,
                      as well as shows it's relevant details in the information panel containing details of molecular formula
                      of the PSMILES, Monomer Molcular Weight, Number of rings in the Monomer, and the corresponding open bond
                      indexes to indicate wher the potential bonds can be formed to connect with the next monomer.`,
        component: <PolyEComponent />,
        docs: <DocRenderer filePath="/markdown_repo/PolymerExplorer.md"/>,
    }
]

export default ExplorerTabContent;