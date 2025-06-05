# AlphaFold Annotations Viewer

A React component for viewing and visualizing AlphaFold protein annotations using UniProt accession keys.

## Features

- **Interactive Input**: Enter UniProt accession keys with optional residue range specifications
- **Annotations Visualization**: Interactive RCSB Saguaro-based viewer showing color-coded annotation tracks
- **Summary Display**: Quick overview of annotation types and counts
- **Complete Data Access**: Full API response data in a structured, collapsible viewer
- **Flexible Integration**: Can be used with pre-loaded data or as a standalone input component

## Installation & Dependencies

This component relies on the following dependencies:
- `@rcsb/rcsb-saguaro` - For sequence/annotations visualization
- `framer-motion` - For animations
- `react-icons` - For UI icons
- Custom components from the project (GlassyContainer, SimpleInputBox, DataViewer)

## Usage

### Basic Usage

```jsx
import AlphafoldAnnotationsViewer from './AlphafoldAnnotationsViewer';

// Simple usage with input box
<AlphafoldAnnotationsViewer />
```

### With Initial Data

```jsx
// Pre-populate with a UniProt accession
<AlphafoldAnnotationsViewer 
    initialAccessionKey="Q5VSL9"
/>

// With specific residue range
<AlphafoldAnnotationsViewer 
    initialAccessionKey="Q5VSL9"
    startPosition={1}
    endPosition={100}
/>
```

### With Pre-loaded Data

```jsx
// Display data without input box
<AlphafoldAnnotationsViewer 
    toolData={annotationsData}
    hideInputBox={true}
/>
```

## Props

| Prop | Type | Default | Description |
|------|------|---------|-------------|
| `toolData` | Object | `null` | Pre-loaded annotations data from API |
| `initialAccessionKey` | String | `""` | Initial UniProt accession to display |
| `hideInputBox` | Boolean | `false` | Whether to hide the input interface |
| `startPosition` | String/Number | `null` | Initial start position for residue range |
| `endPosition` | String/Number | `null` | Initial end position for residue range |

## API Integration

The component uses the `getAlphafoldAnnotations` function from `mcpToolsService`:

```javascript
import { getAlphafoldAnnotations } from './services/api/mcpToolsService';

const response = await getAlphafoldAnnotations({
    uniprot_accession: 'Q5VSL9',
    start: 1,        // Optional
    end: 100         // Optional
});
```

## Data Format

The component expects annotations data in the following format:

```javascript
{
    annotations: [
        {
            type: "Domain",
            start: 10,
            end: 50,
            description: "DNA-binding domain",
            label: "DBD"
        },
        {
            type: "Active Site", 
            start: 25,
            end: 25,
            description: "Catalytic residue",
            label: "Active"
        }
        // ... more annotations
    ]
}
```

## Styling

The component uses CSS classes from `UniProtViewer.css`. Key styles include:
- `.explore-container` - Main container
- `.explorer-row-*` - Row layouts
- `.loading-spinner` - Loading animation
- `.annotation-type-summary` - Summary card styling

## Example Integration

```jsx
import React, { useState } from 'react';
import AlphafoldAnnotationsViewer from './AlphafoldAnnotationsViewer';
import { getAlphafoldAnnotations } from './services/api/mcpToolsService';

const ProteinAnalysisPage = () => {
    const [annotationsData, setAnnotationsData] = useState(null);

    const handleFetchAnnotations = async (accession) => {
        try {
            const response = await getAlphafoldAnnotations({ 
                uniprot_accession: accession 
            });
            setAnnotationsData(response.result);
        } catch (error) {
            console.error('Failed to fetch annotations:', error);
        }
    };

    return (
        <div>
            <h1>Protein Analysis Dashboard</h1>
            
            {/* Standalone viewer */}
            <AlphafoldAnnotationsViewer />
            
            {/* Integrated with custom logic */}
            {annotationsData && (
                <AlphafoldAnnotationsViewer 
                    toolData={annotationsData}
                    hideInputBox={true}
                />
            )}
        </div>
    );
};
```

## Development Notes

- The component uses `useRef` and `useEffect` to manage the RCSB Saguaro viewer lifecycle
- Color palette for annotation tracks is automatically assigned based on annotation types
- The viewer automatically calculates position ranges from annotation data
- Error handling includes both API errors and data processing errors

## Related Components

- `UniProtSummaryViewer` - For protein summary information
- `UniProtViewer` - For basic protein structure viewing
- `DataViewer` - For displaying raw API response data
- `SimpleInputBox` - For user input interface
