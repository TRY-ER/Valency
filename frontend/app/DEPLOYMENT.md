# Valency Frontend Deployment Guide

This guide covers building and deploying the Valency React SPA to nginx web server.

## Quick Start

### 1. One-Command Setup
```bash
# Build and deploy in one go
./setup.sh --domain your-domain.com

# With SSL support
./setup.sh --domain your-domain.com --ssl

# With bundle analysis
./setup.sh --domain your-domain.com --analyze
```

### 2. Step-by-Step Deployment

#### Build the Application
```bash
# Production build
./build.sh --type production

# Development build
./build.sh --type development

# Production build with analysis
./build.sh --type production --analyze
```

#### Deploy to Nginx
```bash
# Basic deployment
./deploy.sh --domain your-domain.com

# With SSL
./deploy.sh --domain your-domain.com --ssl

# Custom web root
./deploy.sh --domain your-domain.com --web-root /opt/valency

# Dry run (see what would happen)
./deploy.sh --domain your-domain.com --dry-run
```

## Prerequisites

### System Requirements
- **Node.js**: 18 or later
- **nginx**: Latest stable version
- **Linux/Unix system** with systemd
- **sudo/root access** for nginx configuration

### Install Dependencies
```bash
# Ubuntu/Debian
sudo apt update
sudo apt install nodejs npm nginx

# CentOS/RHEL/Fedora
sudo dnf install nodejs npm nginx

# or using Node Version Manager (recommended)
curl -o- https://raw.githubusercontent.com/nvm-sh/nvm/v0.39.0/install.sh | bash
nvm install 20
nvm use 20
```

## Script Documentation

### build.sh

Builds the React application with optimizations for nginx deployment.

**Options:**
- `-t, --type TYPE`: Build type (development, production)
- `-a, --analyze`: Analyze bundle after build
- `-c, --clean`: Clean cache before build
- `-d, --deploy TARGET`: Deploy after build (netlify, vercel, docker)
- `-h, --help`: Show help

**Examples:**
```bash
./build.sh --type production --analyze
./build.sh --clean --deploy docker
./build.sh --type development
```

### deploy.sh

Deploys the built application to nginx web server.

**Options:**
- `-w, --web-root PATH`: Web root directory (default: /var/www/valency)
- `-s, --site NAME`: Nginx site name (default: valency)
- `-d, --domain DOMAIN`: Domain name (default: localhost)
- `-n, --dry-run`: Preview changes without executing
- `--no-backup`: Skip backup of existing files
- `--no-restart`: Don't restart nginx after deployment
- `--ssl`: Configure SSL/HTTPS
- `-f, --force`: Force deployment even if checks fail
- `-h, --help`: Show help

**Examples:**
```bash
./deploy.sh --domain myapp.com --ssl
./deploy.sh --web-root /opt/valency --dry-run
./deploy.sh --site valency-app --no-restart
```

### setup.sh

Combines build and deployment in a single command.

**Options:**
- `--domain DOMAIN`: Domain name (default: localhost)
- `--ssl`: Enable SSL/HTTPS
- `--analyze`: Analyze bundle after build
- `--help`: Show help

## Configuration Files

### Environment Files
- `.env.development`: Development configuration
- `.env.production`: Production configuration
- `.env.local`: Local overrides (gitignored)

### Build Configuration
- `config-overrides.js`: Webpack customizations
- `jsconfig.json`: TypeScript/JavaScript configuration
- `package.json`: Dependencies and scripts

### Deployment Configuration
- `nginx.conf`: Nginx configuration template
- `_redirects`: Netlify deployment configuration
- `Dockerfile`: Docker deployment

## SPA Routing Configuration

The application uses React Router with the following routes:

- `/` - Home page
- `/explorer` - Structure Analysis tools
- `/identification` - Identification tools (formerly /discriminators)
- `/optimization` - Optimization tools (formerly /generators)
- `/chatbot` - Master Agent interface
- `/activities` - Activities dashboard
- `/profile` - User profile
- `/settings` - Application settings
- `/login`, `/signup` - Authentication
- `/about` - About page

### Nginx SPA Configuration

The deployment script automatically configures nginx for SPA routing:

```nginx
# Handle client-side routing
location / {
    try_files $uri $uri/ @fallback;
}

# Fallback for SPA routing
location @fallback {
    rewrite ^.*$ /index.html last;
}

# Handle specific SPA routes explicitly
location ~ ^/(explorer|identification|optimization|chatbot|activities|profile|settings|login|signup|about)(/.*)?$ {
    try_files $uri $uri/ /index.html;
}
```

## SSL Configuration

When using the `--ssl` option, you need to provide SSL certificates:

```bash
# Place your certificates here:
/etc/ssl/certs/your-domain.com.crt
/etc/ssl/private/your-domain.com.key

# Or use Let's Encrypt:
sudo certbot --nginx -d your-domain.com
```

## Troubleshooting

### Common Issues

**1. Build Fails**
```bash
# Clean and rebuild
rm -rf node_modules package-lock.json
npm install
./build.sh --clean --type production
```

**2. Nginx Permission Errors**
```bash
# Fix ownership
sudo chown -R www-data:www-data /var/www/valency
sudo chmod -R 755 /var/www/valency
```

**3. 404 Errors on Routes**
```bash
# Check nginx configuration
sudo nginx -t
sudo systemctl reload nginx

# Verify SPA routing is configured
curl -I http://your-domain.com/explorer
```

**4. Static Assets Not Loading**
```bash
# Check file permissions
ls -la /var/www/valency/static/

# Verify cache headers
curl -I http://your-domain.com/static/js/main.js
```

### Logs and Monitoring

```bash
# Nginx error logs
sudo tail -f /var/log/nginx/error.log

# Nginx access logs
sudo tail -f /var/log/nginx/access.log

# Check nginx status
sudo systemctl status nginx

# Test nginx configuration
sudo nginx -t
```

## Performance Optimization

### Build Optimizations
- Code splitting by route and library
- Asset optimization with content hashing
- Gzip compression enabled
- Cache headers configured for static assets

### Runtime Optimizations
- React.lazy() for route-based code splitting
- Service worker for caching (production)
- Optimized bundle sizes with webpack

### Nginx Optimizations
- Gzip compression enabled
- Static asset caching with proper headers
- Security headers configured
- HTTP/2 support (with SSL)

## Security Considerations

### Headers Configured
- `X-Frame-Options: DENY`
- `X-XSS-Protection: 1; mode=block`
- `X-Content-Type-Options: nosniff`
- `Referrer-Policy: strict-origin-when-cross-origin`
- `Permissions-Policy: geolocation=(), microphone=(), camera=()`

### Best Practices
- Use HTTPS in production (`--ssl` option)
- Keep nginx and Node.js updated
- Regular security audits: `npm audit`
- Environment-specific configurations

## Backup and Recovery

The deployment script automatically creates backups:

```bash
# Backups are stored in:
/tmp/valency-backup-YYYYMMDD-HHMMSS/

# To restore from backup:
sudo cp -r /tmp/valency-backup-*/var/www/valency /var/www/
sudo systemctl reload nginx
```

## CI/CD Integration

### GitHub Actions Example
```yaml
name: Deploy to Production
on:
  push:
    branches: [main]
jobs:
  deploy:
    runs-on: ubuntu-latest
    steps:
      - uses: actions/checkout@v2
      - name: Setup Node.js
        uses: actions/setup-node@v2
        with:
          node-version: '20'
      - name: Build
        run: ./build.sh --type production
      - name: Deploy
        run: ./deploy.sh --domain ${{ secrets.DOMAIN }}
```

## Support

For issues and questions:
1. Check the troubleshooting section above
2. Review nginx error logs
3. Verify all prerequisites are met
4. Test with `--dry-run` first

## License

This deployment configuration is part of the Valency project.
