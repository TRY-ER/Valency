# Molecule Explorer Component

## Overview

The **Molecule Explorer** tool allows users to input a SMILES string (Simplified Molecular Input Line Entry System) representing a chemical molecule and provides an interactive interface for visualizing its 2D structure. In addition, the component displays key molecular properties in an information panel.

## Features

- **SMILES Input**: The tool accepts a SMILES string as input, which is used to generate a 2D representation of the molecule.
- **2D Molecular Visualization**: A 2D image of the molecule is rendered, showcasing its chemical structure.
- **Information Panel**: Displays additional molecular properties and metrics.

### Molecular Properties Displayed

1. **Molecular Formula**: The chemical formula of the molecule.
2. **Molecular Weight**: The total molecular weight of the compound.
3. **Heavy Atoms Count**: The number of non-hydrogen atoms in the molecule.
4. **H Bond Donor Count**: The count of hydrogen bond donors.
5. **H Bond Acceptor Count**: The count of hydrogen bond acceptors.
6. **Rotatable Bonds Count**: The number of rotatable bonds within the molecule.
7. **Topological Polar Surface Area (TPSA)**: A property that predicts the ability of the molecule to pass through cell membranes.
8. **Number of Rings**: The total number of ring structures present in the molecule.

## Usage

The `Molecule Explorer` component can be used by passing a SMILES string to the interface. Once the SMILES is inputted, the tool will generate a 2D chemical representation and populate the information panel with the molecular properties listed above.

