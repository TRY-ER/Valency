# ADMET Predictor

## Overview
The ADMET Predictor is a comprehensive molecular property prediction tool that uses advanced machine learning models to predict Absorption, Distribution, Metabolism, Excretion, and Toxicity (ADMET) properties of chemical compounds. By inputting a SMILES string, users can obtain detailed predictions about how a molecule behaves in biological systems, making it invaluable for drug discovery and safety assessment.

## Features

### Core ADMET Predictions
- **Absorption Properties**: Bioavailability, human intestinal absorption, Caco-2 permeability
- **Distribution Properties**: Blood-brain barrier penetration, plasma protein binding, volume of distribution
- **Metabolism Properties**: CYP enzyme interactions, clearance rates, half-life predictions
- **Excretion Properties**: Clearance mechanisms and elimination pathways
- **Toxicity Assessment**: AMES mutagenicity, hERG cardiotoxicity, drug-induced liver injury

### Physicochemical Properties
- **Basic Descriptors**: Molecular weight, LogP, hydrogen bond donors/acceptors
- **Drug-Likeness**: Lipinski's Rule of Five compliance, QED score
- **Structural Features**: Stereo centers, topological polar surface area, molecular complexity
- **Solubility Predictions**: Aqueous solubility and lipophilicity measurements

### Advanced Toxicity Screening
- **Mutagenicity Testing**: AMES test predictions for genetic toxicity
- **Carcinogenicity Assessment**: Cancer-causing potential evaluation
- **Nuclear Receptor Interactions**: Endocrine disruption potential
- **Stress Response**: Cellular stress pathway activation predictions

### Benchmarking Analysis
- **DrugBank Percentiles**: Compare your molecule against approved drugs
- **Relative Positioning**: Understand how your compound ranks in pharmaceutical space
- **Property Distribution**: See where your molecule fits in the drug landscape

## Data Source
**Machine Learning Models** - Trained on pharmaceutical and toxicological databases
- Models trained on thousands of experimental measurements
- Incorporates data from DrugBank, ChEMBL, and toxicology databases
- Continuous model updates with latest research findings
- Validated against experimental data for accuracy

## Usage

### Step 1: Enter Molecular Structure
1. **Input SMILES**: Enter the SMILES string of your molecule in the input field
   - Example demo molecule: `CC(C)CC1=CC=C(C=C1)C(C)C(=O)O` (Ibuprofen)
   - The system validates SMILES format in real-time
2. **Structure Validation**: Green checkmark appears for valid SMILES
3. **Molecular Visualization**: 2D structure display shows your input molecule

### Step 2: Generate Predictions
1. **Click "Predict ADMET Properties"**: Initiates the prediction process
2. **Processing Indicator**: Loading spinner shows prediction is in progress
3. **Wait for Results**: Typical prediction time is 10-30 seconds

### Step 3: Analyze Results
1. **Comprehensive Report**: Detailed ADMET profile appears below
2. **Data Viewer**: Expandable sections for different property categories
3. **Reset Option**: Clear results and start with a new molecule

## What You'll See

### Input Interface
- **SMILES Entry Field**: Text box with placeholder showing demo molecule
- **Header**: "Enter Molecular SMILES for ADMET Prediction"
- **Validation Status**: Real-time feedback on SMILES validity
- **2D Structure**: Visual representation of your input molecule

### Prediction Process
- **Loading Animation**: Spinner with "Predicting ADMET properties..." message
- **Button States**: 
  - "Enter Valid Molecule to Predict" (disabled)
  - "Predict ADMET Properties" (active)
  - "Reset" (after prediction)

### Results Display
When predictions are complete, you'll see a comprehensive data viewer with:

#### Basic Molecular Properties
```
Molecular Weight: 46.069 g/mol
LogP: -0.001 (lipophilicity)
Hydrogen Bond Acceptors: 1
Hydrogen Bond Donors: 1
Lipinski Rule Compliance: 4/4 (drug-like)
QED Score: 0.407 (drug-likeness)
Stereo Centers: 0
TPSA: 20.23 Ų (polar surface area)
```

#### Absorption Predictions
```
Human Intestinal Absorption: 99.7% (excellent)
Bioavailability: 87.8% (good)
Caco-2 Permeability: -4.07 (moderate)
PAMPA Permeability: 63.4% (good)
```

#### Distribution Predictions
```
Blood-Brain Barrier: 97.2% (crosses easily)
Plasma Protein Binding: 4.0% (low binding)
Volume of Distribution: 5.0 L/kg
```

#### Metabolism Predictions
```
CYP1A2 Inhibition: 0.07% (unlikely)
CYP2C19 Inhibition: 0.11% (unlikely)
CYP2C9 Substrate: 9.1% (unlikely)
CYP2D6 Substrate: 9.4% (unlikely)
CYP3A4 Substrate: 20.0% (possible)
Hepatic Clearance: 24.5 mL/min/kg
Half-Life: 19.5 hours
```

#### Toxicity Predictions
```
AMES Mutagenicity: 4.4% (low risk)
hERG Cardiotoxicity: 1.1% (low risk)
Drug-Induced Liver Injury: 9.8% (low risk)
Carcinogenicity: 44.6% (moderate concern)
Clinical Toxicity: 0.05% (very low risk)
Skin Sensitization: 21.6% (moderate risk)
```

#### DrugBank Percentiles
```
Molecular Weight Percentile: 1.1% (very small)
LogP Percentile: 23.1% (hydrophilic)
Bioavailability Percentile: 67.0% (good)
BBB Penetration Percentile: 82.6% (high)
[Additional comparative metrics...]
```

## Examples

### Example 1: Small Drug Molecule (Ethanol - CCO)
**Properties**: Very small, highly soluble, crosses BBB easily
**ADMET Profile**: Good absorption, high BBB penetration, low toxicity
**Drug-Likeness**: Below typical drug molecular weight range

### Example 2: Anti-inflammatory Drug (Ibuprofen)
**SMILES**: `CC(C)CC1=CC=C(C=C1)C(C)C(=O)O`
**Properties**: Moderate size, good lipophilicity, drug-like
**ADMET Profile**: Good oral bioavailability, low toxicity, moderate metabolism

### Example 3: Large Complex Molecule
**Properties**: High molecular weight, multiple functional groups
**ADMET Profile**: May show poor absorption, potential toxicity concerns

## Use Cases

### Drug Discovery
- **Lead Optimization**: Improve ADMET properties of drug candidates
- **Compound Prioritization**: Select compounds with favorable profiles for testing
- **Safety Assessment**: Early identification of potential toxicity issues
- **Formulation Planning**: Understand absorption and distribution challenges

### Chemical Safety
- **Toxicity Screening**: Assess potential health risks of chemicals
- **Environmental Impact**: Predict bioaccumulation and persistence
- **Regulatory Compliance**: Support safety documentation requirements
- **Risk Assessment**: Quantify potential exposure risks

### Academic Research
- **Property Prediction**: Understand structure-activity relationships
- **Method Validation**: Compare predictions with experimental data
- **Chemical Space Exploration**: Map ADMET landscapes of compound libraries
- **Teaching Tool**: Demonstrate drug discovery principles

### Pharmaceutical Development
- **Candidate Selection**: Choose compounds for clinical development
- **Dose Prediction**: Estimate required dosing regimens
- **Formulation Strategy**: Design appropriate delivery systems
- **Clinical Planning**: Anticipate potential side effects and monitoring needs

## Understanding the Results

### Interpreting Percentiles
- **Low Percentiles (0-25%)**: Below typical drug ranges
- **Medium Percentiles (25-75%)**: Within normal drug space
- **High Percentiles (75-100%)**: Above typical drug ranges

### Risk Assessment Guidelines
- **Toxicity Scores < 10%**: Generally considered low risk
- **Toxicity Scores 10-50%**: Moderate concern, requires evaluation
- **Toxicity Scores > 50%**: High concern, likely problematic

### Drug-Likeness Indicators
- **Lipinski Rule**: 4/4 = fully compliant, drug-like
- **QED Score**: 0.0-1.0 scale, >0.5 considered drug-like
- **Bioavailability**: >30% typically required for oral drugs

## Technical Notes
- **Prediction Speed**: Results typically available within 30 seconds
- **Model Accuracy**: Validated against experimental data with high correlation
- **Coverage**: Optimized for small molecule drugs and drug-like compounds
- **Limitations**: Best performance for molecules similar to training data
- **Updates**: Models continuously improved with new experimental data

---

*This tool provides computational predictions based on machine learning models. Results should be used as guidance alongside experimental validation. Predictions may vary for novel chemical scaffolds outside the training data domain.*
