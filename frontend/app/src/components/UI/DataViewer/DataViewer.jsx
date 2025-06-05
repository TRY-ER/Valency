import React, { useState } from 'react';
import './DataViewer.css';

const DataViewer = ({ data, title = "Research Data Explorer", maxDepth = 10, initiallyExpanded = true }) => {
    const [expandedKeys, setExpandedKeys] = useState(new Set());
    const [expandedStrings, setExpandedStrings] = useState(new Set());
    const [copiedItems, setCopiedItems] = useState(new Set());

    // Initialize all keys as expanded if initiallyExpanded is true
    React.useEffect(() => {
        if (initiallyExpanded && data) {
            const allKeys = new Set();
            const allStringPaths = new Set();
            const collectKeys = (obj, path = '') => {
                if (Array.isArray(obj)) {
                    // Add the array path itself
                    if (path) allKeys.add(path);
                    // Add each array item
                    obj.forEach((item, index) => {
                        const itemPath = path ? `${path}[${index}]` : `[${index}]`;
                        if (typeof item === 'object' && item !== null) {
                            collectKeys(item, itemPath);
                        } else if (typeof item === 'string' && item.length > 50) {
                            allStringPaths.add(`${itemPath}_string`);
                        }
                    });
                } else if (typeof obj === 'object' && obj !== null) {
                    // Add the object path itself
                    if (path) allKeys.add(path);
                    // Add each object property
                    Object.keys(obj).forEach(key => {
                        const fullPath = path ? `${path}.${key}` : key;
                        allKeys.add(fullPath);
                        if (typeof obj[key] === 'object' && obj[key] !== null) {
                            collectKeys(obj[key], fullPath);
                        } else if (typeof obj[key] === 'string' && obj[key].length > 50) {
                            allStringPaths.add(`${fullPath}_string`);
                        }
                    });
                }
            };
            collectKeys(data, 'root');
            setExpandedKeys(allKeys);
            setExpandedStrings(allStringPaths);
        }
    }, [data, initiallyExpanded]);

    const toggleKey = (keyPath) => {
        const newExpandedKeys = new Set(expandedKeys);
        if (newExpandedKeys.has(keyPath)) {
            newExpandedKeys.delete(keyPath);
            // Also collapse all nested keys
            Array.from(newExpandedKeys).forEach(key => {
                if (key.startsWith(keyPath + '.')) {
                    newExpandedKeys.delete(key);
                }
            });
        } else {
            newExpandedKeys.add(keyPath);
        }
        setExpandedKeys(newExpandedKeys);
    };

    const expandAll = () => {
        const allKeys = new Set();
        const allStringPaths = new Set();
        const collectKeys = (obj, path = '') => {
            if (Array.isArray(obj)) {
                // Add the array path itself
                if (path) allKeys.add(path);
                // Add each array item
                obj.forEach((item, index) => {
                    const itemPath = path ? `${path}[${index}]` : `[${index}]`;
                    if (typeof item === 'object' && item !== null) {
                        collectKeys(item, itemPath);
                    } else if (typeof item === 'string' && item.length > 50) {
                        allStringPaths.add(`${itemPath}_string`);
                    }
                });
            } else if (typeof obj === 'object' && obj !== null) {
                // Add the object path itself
                if (path) allKeys.add(path);
                // Add each object property
                Object.keys(obj).forEach(key => {
                    const fullPath = path ? `${path}.${key}` : key;
                    allKeys.add(fullPath);
                    if (typeof obj[key] === 'object' && obj[key] !== null) {
                        collectKeys(obj[key], fullPath);
                    } else if (typeof obj[key] === 'string' && obj[key].length > 50) {
                        allStringPaths.add(`${fullPath}_string`);
                    }
                });
            }
        };
        collectKeys(data, 'root');
        setExpandedKeys(allKeys);
        setExpandedStrings(allStringPaths);
    };

    const collapseAll = () => {
        setExpandedKeys(new Set());
        setExpandedStrings(new Set());
    };

    const copyToClipboard = async (text, stringPath) => {
        try {
            await navigator.clipboard.writeText(text);
            // Visual feedback
            const newCopiedItems = new Set(copiedItems);
            newCopiedItems.add(stringPath);
            setCopiedItems(newCopiedItems);
            // Clear feedback after 2 seconds
            setTimeout(() => {
                setCopiedItems(prev => {
                    const updated = new Set(prev);
                    updated.delete(stringPath);
                    return updated;
                });
            }, 2000);
        } catch (err) {
            // Fallback for older browsers
            const textArea = document.createElement('textarea');
            textArea.value = text;
            document.body.appendChild(textArea);
            textArea.select();
            document.execCommand('copy');
            document.body.removeChild(textArea);
            // Same visual feedback for fallback
            const newCopiedItems = new Set(copiedItems);
            newCopiedItems.add(stringPath);
            setCopiedItems(newCopiedItems);
            setTimeout(() => {
                setCopiedItems(prev => {
                    const updated = new Set(prev);
                    updated.delete(stringPath);
                    return updated;
                });
            }, 2000);
        }
    };

    const toggleStringExpansion = (stringPath) => {
        const newExpandedStrings = new Set(expandedStrings);
        if (newExpandedStrings.has(stringPath)) {
            newExpandedStrings.delete(stringPath);
        } else {
            newExpandedStrings.add(stringPath);
        }
        setExpandedStrings(newExpandedStrings);
    };

    const getValueType = (value) => {
        if (value === null) return 'null';
        if (Array.isArray(value)) return 'array';
        return typeof value;
    };

    const getValueTypeClass = (value) => {
        const type = getValueType(value);
        return `value-${type}`;
    };

    const isValidUrl = (string) => {
        // Basic checks first
        if (!string || typeof string !== 'string' || string.length < 7) {
            return false;
        }
        
        // Check if it starts with http:// or https://
        if (!string.match(/^https?:\/\//i)) {
            return false;
        }
        
        try {
            const url = new URL(string);
            return url.protocol === 'http:' || url.protocol === 'https:';
        } catch {
            return false;
        }
    };

    const renderValue = (value, keyPath = '', depth = 0) => {
        if (depth > maxDepth) {
            return <span className="value-truncated">... (max depth reached)</span>;
        }

        if (value === null) {
            return <span className="value-null">null</span>;
        }

        if (value === undefined) {
            return <span className="value-undefined">undefined</span>;
        }

        if (typeof value === 'string') {
            const isUrl = isValidUrl(value);
            const isLongString = value.length > 50;
            const stringPath = `${keyPath}_string`;
            const isExpanded = expandedStrings.has(stringPath);
            const isCopied = copiedItems.has(stringPath);
            const displayValue = isLongString && !isExpanded ? `${value.substring(0, 47)}...` : value;
            
            if (isUrl) {
                return (
                    <div className="string-value-container">
                        <span className={`value-string value-url ${isExpanded ? 'expanded' : ''} ${isLongString && !isExpanded ? 'truncated' : ''}`}>
                            "<a 
                                href={value} 
                                target="_blank" 
                                rel="noopener noreferrer"
                                title={isLongString ? value : `Open ${value} in new tab`}
                                className="url-link"
                            >
                                {displayValue}
                            </a>"
                        </span>
                        {isLongString && (
                            <div className="string-controls">
                                <button 
                                    className="string-control-btn expand-btn"
                                    onClick={() => toggleStringExpansion(stringPath)}
                                    title={isExpanded ? "Collapse" : "Expand"}
                                >
                                    {isExpanded ? "▲" : "▼"}
                                </button>
                                <button 
                                    className={`string-control-btn copy-btn ${isCopied ? 'copied' : ''}`}
                                    onClick={() => copyToClipboard(value, stringPath)}
                                    title={isCopied ? "Copied!" : "Copy to clipboard"}
                                >
                                    {isCopied ? "✓" : "📋"}
                                </button>
                            </div>
                        )}
                    </div>
                );
            }
            
            return (
                <div className="string-value-container">
                    <span className={`value-string ${isExpanded ? 'expanded' : ''} ${isLongString && !isExpanded ? 'truncated' : ''}`} title={isLongString && !isExpanded ? value : undefined}>
                        "{displayValue}"
                    </span>
                    {isLongString && (
                        <div className="string-controls">
                            <button 
                                className="string-control-btn expand-btn"
                                onClick={() => toggleStringExpansion(stringPath)}
                                title={isExpanded ? "Collapse" : "Expand"}
                            >
                                {isExpanded ? "▲" : "▼"}
                            </button>
                            <button 
                                className={`string-control-btn copy-btn ${isCopied ? 'copied' : ''}`}
                                onClick={() => copyToClipboard(value, stringPath)}
                                title={isCopied ? "Copied!" : "Copy to clipboard"}
                            >
                                {isCopied ? "✓" : "📋"}
                            </button>
                        </div>
                    )}
                </div>
            );
        }

        if (typeof value === 'number') {
            const isLargeNumber = Math.abs(value) >= 1000000 || (Math.abs(value) < 0.001 && value !== 0);
            const formattedValue = isLargeNumber ? value.toExponential(2) : value.toString();
            return <span className={`value-number ${isLargeNumber ? 'scientific' : ''}`}>{formattedValue}</span>;
        }

        if (typeof value === 'boolean') {
            return (
                <span className="value-boolean" title={`Boolean value: ${value}`}>
                    {value.toString()}
                </span>
            );
        }

        if (Array.isArray(value)) {
            const isExpanded = expandedKeys.has(keyPath);
            return (
                <div className="array-container">
                    <div className="array-header" onClick={() => toggleKey(keyPath)}>
                        <span className="collapse-icon">{isExpanded ? '▼' : '▶'}</span>
                        <span className="array-label">Dataset ({value.length} entries)</span>
                    </div>
                    {isExpanded && (
                        <div className="array-content">
                            {value.map((item, index) => (
                                <div key={index} className="array-item">
                                    <span className="array-index">#{index + 1}:</span>
                                    {renderValue(item, `${keyPath}[${index}]`, depth + 1)}
                                </div>
                            ))}
                        </div>
                    )}
                </div>
            );
        }

        if (typeof value === 'object') {
            const isExpanded = expandedKeys.has(keyPath);
            const keys = Object.keys(value);
            
            return (
                <div className="object-container">
                    <div className="object-header" onClick={() => toggleKey(keyPath)}>
                        <span className="collapse-icon">{isExpanded ? '▼' : '▶'}</span>
                        <span className="object-label">Record ({keys.length} fields)</span>
                    </div>
                    {isExpanded && (
                        <div className="object-content">
                            {keys.map(key => {
                                const childKeyPath = keyPath ? `${keyPath}.${key}` : key;
                                return (
                                    <div key={key} className="object-property">
                                        <span className="property-key">{key}:</span>
                                        <div className="property-value">
                                            {renderValue(value[key], childKeyPath, depth + 1)}
                                        </div>
                                    </div>
                                );
                            })}
                        </div>
                    )}
                </div>
            );
        }

        return <span className="value-unknown">{String(value)}</span>;
    };

    const renderStats = () => {
        if (!data) return null;
        
        const getStats = (obj) => {
            let totalKeys = 0;
            let totalArrays = 0;
            let totalObjects = 0;
            
            const traverse = (value) => {
                if (Array.isArray(value)) {
                    totalArrays++;
                    value.forEach(traverse);
                } else if (typeof value === 'object' && value !== null) {
                    totalObjects++;
                    const keys = Object.keys(value);
                    totalKeys += keys.length;
                    keys.forEach(key => traverse(value[key]));
                }
            };
            
            traverse(obj);
            return { totalKeys, totalArrays, totalObjects };
        };

        const stats = getStats(data);
        return (
            <div className="data-stats">
                <span>Fields: {stats.totalKeys}</span>
                <span>Records: {stats.totalObjects}</span>
                <span>Datasets: {stats.totalArrays}</span>
            </div>
        );
    };

    if (!data) {
        return (
            <div className="data-viewer">
                <div className="data-viewer-header">
                    <h3>{title}</h3>
                </div>
                <div className="no-data">No research data available for analysis</div>
            </div>
        );
    }

    return (
        <div className="data-viewer">
            <div className="data-viewer-header">
                <h3>{title}</h3>
                <div className="data-viewer-controls">
                    {renderStats()}
                    <button className="control-btn" onClick={expandAll}>Expand All</button>
                    <button className="control-btn" onClick={collapseAll}>Collapse All</button>
                </div>
            </div>
            <div className="data-viewer-content">
                {renderValue(data, 'root')}
            </div>
        </div>
    );
};

export default DataViewer;
