import React from "react";
import ReactMarkdown from "react-markdown";

// Custom renderer for links. Adjust condition if needed.
const CustomLink = ({ href, children, ...props }) => {
  // If the link text is "Redirect", render a button.
  if (children && children === "Redirect") {
    return (
      <button onClick={() => window.location.href = href} {...props}>
        {children}
      </button>
    );
  }
  return (
    <a href={href} {...props}>
      {children}
    </a>
  );
};

const CustomMarkdownRenderer = ({ content }) => {
  return (
    <ReactMarkdown
      components={{
        a: CustomLink,
      }}
    >
      {content}
    </ReactMarkdown>
  );
};

export default CustomMarkdownRenderer;