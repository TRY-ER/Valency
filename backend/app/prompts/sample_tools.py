

tool_config = """
Following are the tools that are available for use:

    *Remember to return only tools that are functional not the master tools which have sub-tools.*
    
    1.
Tool Name: Explorers
Tool Link: /
Tool Description: Explorers are used to explore molecules, proteins, and polymers. This is the main tool below which there are three more tools consisting tools specfic to molecules, proteins and polymers

2.
Tool Name: Explorers > Molecule Explorer
Tool Link: /
Tool Description: The molecule explorer tool takes SMILES string of the molecule and provides interface for showcasing the 2D image of the chemical representation of the molecule and a information panel for showcasing other properties. The information panel shows details as Molecular Formula, Molecular Weight, Heavy Atoms Count, H Bond Doner Count, H Bond Acceptor Count, Rotatale Bonds Count, Topological Polar Surface Area (TPSA), and Number of rings

3.
Tool Name: Explorers > Protein Explorer
Tool Link: /proe
Tool Description: The protein explorer tool takes the Protein Data Bank ID (PDB ID) to show it's details and 3D structure in two panels. The 3D structure is rotatable and scalable for exploration purposes. The information panel whows it's several details consisting Pdb Id, Title,Authors,Journal,Year,Volume,Pages,Doi,Pubmed Id, Experiment Method,Molecular Weight (kDa),Deposited Model Count,Polymer entity count,Polymer monomer count, Structural Features,Release Date,Resolution

4.
Tool Name: Explorers > Polymer Explorer
Tool Link: /polye
Tool Description: The polymer explorer tool takes PSMILES of the polymer and gives visual representation in a sidepanel, as well as shows it's relevant details in the information panel containing details of molecular formula of the PSMILES, Monomer Molcular Weight, Number of rings in the Monomer, and the corresponding open bond indexes to indicate wher the potential bonds can be formed to connect with the next monomer.

5.
Tool Name: Discriminators
Tool Link: /discriminators
Tool Description: Discriminators tool cosists of tools to discriminate molecules, proteins, and polymers. But as of now we have similarity search implemented as a sub tool. This sub tool contains structural similarity search for molecules, proteins, and polymers handled as sub-sub tools.

6.
Tool Name: Discriminators > Similarity Search
Tool Link: /discriminators
Tool Description: This tool uses vectors generated from SMILES, PSMILES, Protein sequences to conduct a similarity search and get relevant candidates from those. This tool contains three sub-tools for similarity search with SMILES for molecules, PSMILES for polymers, and PDB-ID for proteins.

7.
Tool Name: Discriminators > Similarity Search > SMILES
Tool Link: /discriminators
Tool Description: This tool uses SMILES string to retrieve similar molecules from the chemical space.
                      This tool takes the SMILES string with a the number of candidates to retrieve and does
                      the operation and returns a section that shows the image of the molecules with it's SMILES and 
                      the similarity distance from the query SMILES. The results can be downloadable in a text file
                      as well as can be copied to the clipboard.

8.
Tool Name: Discriminators > Similarity Search > PSMILES
Tool Link: /discriminators/psmiles
Tool Description: This tool uses PSMILES string to retrieve similar polymers from the polymer space.
                      This tool takes the PSMILES string with a the number of candidates to retrieve and does
                      the operation and returns a section that shows the image of the polymers with it's PSMILES and 
                      the similarity distance from the query PSMILES.The results can be downloadable in a text file
                      as well as can be copied to the clipboard.

9.
Tool Name: Discriminators > Similarity Search > Protein
Tool Link: /discriminators/protein
Tool Description: This tool uses PDB ID string to retrieve similar proteins from the PDB space.
                      This tool takes the PDB ID string with a the number of candidates to retrieve and does
                      the operation and returns a section that shows the image of the proteins with it's PDB Id 
                      and similarity score. The results can be downloadable in a text file as well as can be 
                      copied to the clipboard.

10.
Tool Name: Generators
Tool Link: /generators
Tool Description: Generators are used to create hypothetical molecules, polymers (protein generation is not yet implemented). It has two sub-tools name BRICS and LSTM generators. Both of these tools have sub tools for generation for SMILES and PSMILES format for molecules and polymers respectively.

11.
Tool Name: Generators > BRICS Generator
Tool Link: /generators
Tool Description: Generate molecules using the BRICS algorithm. The BRICS algorithm is a rule-based algorithm that generates molecules based on the given input. The input can be a SMILES string or a PSMILES string. The output will be a set of molecules that are generated based on the input. This tool has two sub-tools for generation for SMILES and PSMILES format for molecules and polymers respectively.

12.
Tool Name: Generators > BRICS Generator > SMILES
Tool Link: /generators
Tool Description: This tool uses SMILES string to generate molecules using BRICS algorithm. The algorithm generates hypothetical molecules based on the given input. The input can be a SMILES string in a text file each represented in a new line. The output will be a set of SMILES that are generated based on the input. The tool has three sections, the first section shows a progress bar showing the completed steps for the job. The second section contains the file upload setup with a panel to show first few files of file to ensure the right text file has been selected. The third section is the output section that consists of a output panel, download button and the reset button. The output panel shows the steps details as well as completed steps and the molecules generated. This tool is specifically useful when you have a set of molecules and you want to generate more molecules based on the given set.

13.
Tool Name: Generators > BRICS Generator > PSMILES
Tool Link: /generators/psmiles
Tool Description: This tool uses PSMILES string to generate molecules using BRICS algorithm. The algorithm generates hypothetical polymers based on the given input. The input can be a PSMILES string in a text file each represented in a new line. The output will be a set of PSMILES that are generated based on the input. The tool has three sections, the first section shows a progress bar showing the completed steps for the job. The second section contains the file upload setup with a panel to show first few files of file to ensure the right text file has been selected. The third section is the output section that consists of a output panel, download button and the reset button. The output panel shows the steps details as well as completed steps and the psmiles generated upon completed. This tool is specifically useful when you have a set of polymers and you want to generate more polymers based on the given set.

14.
Tool Name: Generators > LSTM Generator
Tool Link: /generators/lstm
Tool Description: Generate molecules using the LSTM algorithm. The LSTM algorithm is a deep learning algorithm that generates molecules based on the given input.The input can be a SMILES string or a PSMILES string. The output will be a set of molecules that are generated based on the input. This tool has two sub-tools for generation for SMILES and PSMILES format for molecules and polymers respectively.

15.
Tool Name: Generators > LSTM Generator > PSMILES
Tool Link: /generators/lstm
Tool Description: This tool uses PSMILES string to generate molecules using LSTM algorithm. The algorithm generates hypothetical polymers without any given input polymer string. It's very useful when you want to consider entirely new polymer rather refering from a specific kind. The tool takes input integer for number of generations. Then it generates the polymers based on the given number of generations. The tool has three sections, the first section shows a progress bar showing the number of genertions recieved. The second section contains the input setup with a input box to take the number of generations. The third section is the output section that consists of a output panel, download button and the reset button. the output panel shows the steps details as well as completed steps and the psmiles generated upon completed.

16.
Tool Name: Generators > LSTM Generator > WDG
Tool Link: /generators/lstm/wdg
Tool Description: This tool uses weighted directed graph string to generate molecules using LSTM algorithm. The algorithm generates hypothetical polymers without any given input polymer string. It's very useful when you want to consider entirely new polymer rather refering from a specific kind. The tool takes input integer for number of generations. Then it generates the polymers based on the given number of generations. The tool has three sections, the first section shows a progress bar showing the number of genertions recieved. The second section contains the input setup with a input box to take the number of generations. The third section is the output section that consists of a output panel, download button and the reset button. the output panel shows the steps details as well as completed steps and the psmiles generated upon completed.
"""