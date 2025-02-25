# LSTM WDG Generator Documentation
---

## Overview
---

The **LSTM WDG Generator** tool implements a Long Short-Term Memory (LSTM) neural network architecture to generate novel polymer structures using Weighted Directed Graph (WDG) representation. Based on approaches similar to those described in our research paper "Open-source Polymer Generative Pipeline" `[1]`, this tool generates polymer structures represented as weighted directed graphs `[2]` without requiring input polymers.
For better understanding of weighted directed graphs follow [this tutorial](https://deepchem.io/tutorials/understanding-weighted-directed-graphs-for-polymer-implimentations/) `[3]`

*Note - This tool can be significantly slower if hosted on a CPU server (Given it runs a model compatible with pytorch CUDA)*

## Features
---
- **Parameter Input**: Specify the number of polymer structures to generate
- **Progress Tracking**: Visual progress bar showing generation status
- **Results Panel**: Interactive output display showing generated WDG structures
- **Export Options**: Download generated structures or reset the process

### Tool Sections
---
1. **Progress Section**
   - Progress bar showing number of generations completed
   - Current generation count
   - Overall process status

2. **Input Section**
   - Number of generations input box
   - Generation parameter controls
   - Start/Stop generation buttons

3. **Output Section**
   - Generated WDG structures display
   - Download results button
   - Reset operation button
   - Generation log details

## Usage
---
1. Enter the desired number of polymer structures to generate
2. Initiate the generation process
3. Monitor generation progress
4. Review generated structures
5. Export results as needed

## Example
---
To generate 10 novel polymer structures:

1. Enter number of generations: `10`
2. Click "Generate"
3. Results will appear as weighted directed graphs as strings:
```
1. [1*]c1ccc2c(c1)S(=O)(=O)c1cc([2*])ccc1-2.[3*]c1ccc2c(-c3c(O)ccc4cc([4*])ccc34)c(O)ccc2c1|0.25|0.75|<1-2:0.375:0.375<1-1:0.375:0.375<2-2:0.375:0.375<3-4:0.375:0.375<3-3:0.375:0.375<4-4:0.125:0.125<1-3:0.125:0.125<1-4:0.125:0.125<2-3:0.125:0.125<2-4:0.125:0.125

2. [1*]c1ccc2c(c1)S(=O)(=O)c1cc([2*])ccc1-2.[3*]c1cc([4*])cc(Oc2ccccc2)c1|0.25|0.75|<1-3:0.25:0.25<1-4:0.25:0.25<2-3:0.25:0.25<2-4:0.25:0.25<1-2:0.25:0.25<3-4:0.25:0.25<1-1:0.25:0.25<2-2:0.25:0.25<3-3:0.25:0.25<4-4:0.25:0.25

3. [1*]c1ccc([2*])c2nsnc12.[3*]c1cc([4*])cc2[nH]nc(N)c12|0.25|0.75|<1-3:0.25:0.25<1-4:0.25:0.25<2-3:0.25:0.25<2-4:0.25:0.25<1-2:0.25:0.25<3-4:0.25:0.25<1-1:0.25:0.25<2-2:0.25:0.25<3-3:0.25:0.25<4-4:0.25:0.25
```

### Use Cases
---
This tool is particularly useful for:
- Generating complex polymer architectures
- Exploring branched polymer structures
- Designing polymer networks
- Studying polymer topology
- Analyzing connectivity patterns

## References
---
[1] Mohanty, Debasish, et al. "Open-source Polymer Generative Pipeline." arXiv preprint arXiv:2412.08658 (2024).

[2] Aldeghi, Matteo, and Connor W. Coley. "A graph representation of molecular ensembles for polymer property prediction." Chemical Science 13.35 (2022): 10486-10498.

[3] Ramsundar, Bharath, et al. Deep learning for the life sciences: applying deep learning to genomics, microscopy, drug discovery, and more. O'Reilly Media, 2019.

### Notes
---
- Generated structures are hypothetical and require experimental validation
- Results may vary between generation runs
- WDG representation allows for complex polymer architectures beyond linear chains
- Graph structures maintain chemical validity constraints