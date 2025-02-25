import SSContent from "../../pages/Discriminator/SimilaritySearch/SSContent";
import DocRenderer from "../doc_content/DocRenderer";

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
    }
]; 

export {
    SimilaritySearchTabContents,
}