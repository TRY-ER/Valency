# LSTM PSMILES Generator Documentation
---

## Overview
---

The **LSTM PSMILES Generator** tool implements a Long Short-Term Memory (LSTM) neural network architecture to generate novel polymer structures. Based on approaches similar to those described in our research paper "Open-source Polymer Generative Pipeline [1], this tool generates PSMILES (Polymer SMILES) strings representing new polymer structures without requiring input polymers.

*Note - This tool can be significantly slower if hosted on a CPU server (Given it runs a model compatible with pytorch CUDA)*

## Features
---
- **Parameter Input**: Specify the number of polymer structures to generate
- **Progress Tracking**: Visual progress bar showing generation status
- **Results Panel**: Interactive output display showing generated polymers
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
   - Generated polymers display
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
3. Results will appear as:
```
[*]CC(C)C[*]
[*]CC(=O)O[*]
[*]CCC(C)(C)CC[*]
```

### Use Cases
---
This tool is particularly useful for:
- Discovering entirely new polymer structures
- Exploring uncharted regions of polymer space
- Generating diverse polymer libraries
- Initial screening of potential new materials
- Polymer drug discovery applications

## References
---
[1] Mohanty, Debasish, et al. "Open-source Polymer Generative Pipeline." arXiv preprint arXiv:2412.08658 (2024). 

### Notes
---
- Generated structures are hypothetical and require experimental validation
- Results may vary between generation runs
- All structures follow PSMILES notation conventions