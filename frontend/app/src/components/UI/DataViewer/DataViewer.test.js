import { render, screen } from '@testing-library/react';
import DataViewer from './DataViewer';

describe('DataViewer URL handling', () => {
  test('renders clickable link for valid URL strings', () => {
    const testData = {
      website: 'https://example.com',
      api: 'http://api.example.com/v1',
      regularString: 'This is just text'
    };

    render(<DataViewer data={testData} title="Test Data" />);

    // Check that URLs are rendered as links
    const httpsLink = screen.getByRole('link', { name: /https:\/\/example\.com/ });
    const httpLink = screen.getByRole('link', { name: /http:\/\/api\.example\.com\/v1/ });
    
    expect(httpsLink).toBeInTheDocument();
    expect(httpsLink).toHaveAttribute('href', 'https://example.com');
    expect(httpsLink).toHaveAttribute('target', '_blank');
    expect(httpsLink).toHaveAttribute('rel', 'noopener noreferrer');

    expect(httpLink).toBeInTheDocument();
    expect(httpLink).toHaveAttribute('href', 'http://api.example.com/v1');
    expect(httpLink).toHaveAttribute('target', '_blank');
    expect(httpLink).toHaveAttribute('rel', 'noopener noreferrer');

    // Check that regular strings are not rendered as links
    const regularText = screen.getByText(/"This is just text"/);
    expect(regularText).toBeInTheDocument();
    expect(regularText.tagName).toBe('SPAN');
  });

  test('handles invalid URLs as regular strings', () => {
    const testData = {
      invalidUrl1: 'not-a-url',
      invalidUrl2: 'ftp://example.com', // FTP is not supported
      invalidUrl3: 'javascript:alert("test")', // Potential XSS
    };

    render(<DataViewer data={testData} title="Test Invalid URLs" />);

    // These should all be rendered as regular strings, not links
    expect(screen.queryByRole('link')).not.toBeInTheDocument();
    expect(screen.getByText(/"not-a-url"/)).toBeInTheDocument();
    expect(screen.getByText(/"ftp:\/\/example\.com"/)).toBeInTheDocument();
    expect(screen.getByText(/"javascript:alert\(\"test\"\)"/)).toBeInTheDocument();
  });

  test('handles long URLs correctly', () => {
    const longUrl = 'https://example.com/very/long/path/that/exceeds/fifty/characters/in/length/and/should/be/truncated';
    const testData = {
      longUrl: longUrl
    };

    render(<DataViewer data={testData} title="Test Long URL" />);

    const link = screen.getByRole('link');
    expect(link).toBeInTheDocument();
    expect(link).toHaveAttribute('href', longUrl);
    expect(link).toHaveAttribute('title', `Open ${longUrl} in new tab`);
    // The displayed text should be truncated
    expect(link.textContent).toMatch(/^https:\/\/example\.com\/very\/long\/path\/that\/exceed\.\.\./);
  });
});
