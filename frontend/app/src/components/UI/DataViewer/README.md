# DataViewer Component

A React component for displaying complex JavaScript objects and arrays with collapsible/expandable functionality.

## Features

- ✅ **Nested Data Support**: Display deeply nested objects and arrays
- ✅ **Collapsible Keys**: Click to expand/collapse any object or array
- ✅ **Type-Specific Styling**: Different colors for strings, numbers, booleans, null, etc.
- ✅ **Expand/Collapse Controls**: Buttons to expand or collapse all keys at once
- ✅ **Data Statistics**: Shows count of keys, objects, and arrays
- ✅ **Configurable Depth**: Set maximum depth to prevent performance issues
- ✅ **Initial State Control**: Choose whether data starts expanded or collapsed
- ✅ **Responsive Design**: Works on mobile and desktop
- ✅ **Theme Support**: Supports both light and dark themes
- ✅ **Glass Morphism**: Beautiful modern styling with backdrop blur

## Installation

The component is already included in your project structure:

```
src/components/UI/DataViewer/
├── DataViewer.jsx          # Main component
├── DataViewer.css          # Styles
├── DataViewerExample.jsx   # Usage examples
└── index.js               # Export file
```

## Basic Usage

```jsx
import DataViewer from './components/UI/DataViewer/DataViewer';

const MyComponent = () => {
    const data = {
        user: {
            name: "John Doe",
            preferences: {
                theme: "dark",
                notifications: ["email", "push"]
            }
        }
    };

    return (
        <DataViewer 
            data={data}
            title="User Data"
        />
    );
};
```

## Props

| Prop | Type | Default | Description |
|------|------|---------|-------------|
| `data` | any | - | The data to display (object, array, or primitive) |
| `title` | string | "Data Viewer" | Title shown in the header |
| `maxDepth` | number | 10 | Maximum nesting depth to display |
| `initiallyExpanded` | boolean | true | Whether all keys start expanded |

## Examples

### Simple Object
```jsx
<DataViewer 
    data={{ name: "John", age: 30 }}
    title="User Profile"
/>
```

### Array of Objects
```jsx
<DataViewer 
    data={[
        { id: 1, name: "Item 1" },
        { id: 2, name: "Item 2" }
    ]}
    title="Items List"
/>
```

### Complex Nested Data
```jsx
<DataViewer 
    data={complexApiResponse}
    title="API Response"
    initiallyExpanded={false}
    maxDepth={5}
/>
```

### ADMET Integration Example
```jsx
// In your ADMET component
const handlePredict = async () => {
    const response = await getAdmetPrediction({ smiles: activeMol });
    setPredictionResults(response);
};

// In your render method
{predictionResults && (
    <DataViewer 
        data={predictionResults}
        title="ADMET Prediction Results"
        initiallyExpanded={false}
    />
)}
```

## Styling

The component uses CSS custom properties and supports both light and dark themes automatically based on user preference.

Key CSS classes:
- `.data-viewer` - Main container
- `.value-string` - String values (green)
- `.value-number` - Number values (red)
- `.value-boolean` - Boolean values (purple)
- `.value-null` - Null values (gray)
- `.object-container` - Object wrappers
- `.array-container` - Array wrappers

## Performance Considerations

- Use `maxDepth` prop to limit deep nesting
- For very large datasets, consider pagination or virtualization
- The component efficiently handles expand/collapse state without re-rendering the entire tree

## Browser Support

- Modern browsers with CSS backdrop-filter support
- Graceful fallback for older browsers
- Mobile-responsive design

## Integration with Existing Projects

The DataViewer component has been integrated into your ADMET prediction component and will show:

1. **Summary View**: Traditional bar chart for quick visualization
2. **Detailed View**: DataViewer component showing complete API response data

This allows users to see both a high-level overview and drill down into the detailed prediction data as needed.
