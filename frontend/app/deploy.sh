#!/bin/bash

# Valency Frontend Deployment Script
# Deploys the React SPA to nginx web server

set -e

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
CYAN='\033[0;36m'
NC='\033[0m' # No Color

# Configuration - Modify these for your server
DEFAULT_WEB_ROOT="/var/www/valency"
DEFAULT_NGINX_SITE="valency"
DEFAULT_DOMAIN="localhost"
NGINX_SITES_DIR="/etc/nginx/sites-available"
NGINX_ENABLED_DIR="/etc/nginx/sites-enabled"
BUILD_DIR="build"

# Default values
WEB_ROOT="$DEFAULT_WEB_ROOT"
NGINX_SITE="$DEFAULT_NGINX_SITE"
DOMAIN="$DEFAULT_DOMAIN"
DRY_RUN=false
BACKUP=true
RESTART_NGINX=true
SSL=false
FORCE=false

# Function to print colored output
print_status() {
    echo -e "${BLUE}[INFO]${NC} $1"
}

print_success() {
    echo -e "${GREEN}[SUCCESS]${NC} $1"
}

print_warning() {
    echo -e "${YELLOW}[WARNING]${NC} $1"
}

print_error() {
    echo -e "${RED}[ERROR]${NC} $1"
}

print_header() {
    echo -e "${CYAN}================================${NC}"
    echo -e "${CYAN}$1${NC}"
    echo -e "${CYAN}================================${NC}"
}

# Help function
show_help() {
    echo "Usage: $0 [OPTIONS]"
    echo ""
    echo "Options:"
    echo "  -w, --web-root PATH      Web root directory (default: $DEFAULT_WEB_ROOT)"
    echo "  -s, --site NAME          Nginx site name (default: $DEFAULT_NGINX_SITE)"
    echo "  -d, --domain DOMAIN      Domain name (default: $DEFAULT_DOMAIN)"
    echo "  -n, --dry-run           Show what would be done without executing"
    echo "  --no-backup             Skip backup of existing files"
    echo "  --no-restart            Don't restart nginx after deployment"
    echo "  --ssl                   Configure SSL/HTTPS"
    echo "  -f, --force             Force deployment even if checks fail"
    echo "  -h, --help              Show this help message"
    echo ""
    echo "Examples:"
    echo "  $0 --domain myapp.com --ssl"
    echo "  $0 --web-root /opt/valency --dry-run"
    echo "  $0 --site valency-app --no-restart"
    echo ""
    echo "Prerequisites:"
    echo "  - Nginx installed and running"
    echo "  - Build directory exists (run ./build.sh first)"
    echo "  - Appropriate permissions for web directory"
}

# Parse arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        -w|--web-root)
            WEB_ROOT="$2"
            shift 2
            ;;
        -s|--site)
            NGINX_SITE="$2"
            shift 2
            ;;
        -d|--domain)
            DOMAIN="$2"
            shift 2
            ;;
        -n|--dry-run)
            DRY_RUN=true
            shift
            ;;
        --no-backup)
            BACKUP=false
            shift
            ;;
        --no-restart)
            RESTART_NGINX=false
            shift
            ;;
        --ssl)
            SSL=true
            shift
            ;;
        -f|--force)
            FORCE=true
            shift
            ;;
        -h|--help)
            show_help
            exit 0
            ;;
        *)
            print_error "Unknown option: $1"
            show_help
            exit 1
            ;;
    esac
done

# Function to run commands with dry-run support
run_cmd() {
    if [[ "$DRY_RUN" == true ]]; then
        echo -e "${YELLOW}[DRY-RUN]${NC} $1"
    else
        eval "$1"
    fi
}

# Function to check if command exists
command_exists() {
    command -v "$1" >/dev/null 2>&1
}

print_header "Valency Frontend Deployment"

print_status "Configuration:"
echo "  Web root: $WEB_ROOT"
echo "  Nginx site: $NGINX_SITE"
echo "  Domain: $DOMAIN"
echo "  SSL: $SSL"
echo "  Dry run: $DRY_RUN"
echo ""

# Pre-deployment checks
print_header "Pre-deployment Checks"

# Check if running as root or with sudo
if [[ $EUID -ne 0 ]] && [[ -z "$SUDO_USER" ]]; then
    if [[ "$FORCE" != true ]]; then
        print_error "This script requires root privileges to modify nginx and web directories."
        print_error "Run with sudo or use --force to skip this check."
        exit 1
    else
        print_warning "Running without root privileges - some operations may fail"
    fi
fi

# Check if build directory exists
if [[ ! -d "$BUILD_DIR" ]]; then
    print_error "Build directory '$BUILD_DIR' not found."
    print_error "Run './build.sh' first to build the application."
    exit 1
fi

print_success "✓ Build directory found"

# Check if index.html exists
if [[ ! -f "$BUILD_DIR/index.html" ]]; then
    print_error "index.html not found in build directory"
    exit 1
fi

print_success "✓ index.html found"

# Check if nginx is installed
if ! command_exists nginx; then
    print_error "Nginx is not installed or not in PATH"
    exit 1
fi

print_success "✓ Nginx is installed"

# Check nginx configuration
if ! nginx -t >/dev/null 2>&1; then
    if [[ "$FORCE" != true ]]; then
        print_error "Nginx configuration test failed"
        print_error "Fix nginx configuration or use --force to skip this check"
        exit 1
    else
        print_warning "Nginx configuration test failed - continuing anyway"
    fi
else
    print_success "✓ Nginx configuration is valid"
fi

# Backup existing site if it exists
if [[ "$BACKUP" == true ]] && [[ -d "$WEB_ROOT" ]]; then
    print_header "Creating Backup"
    BACKUP_DIR="/tmp/valency-backup-$(date +%Y%m%d-%H%M%S)"
    print_status "Backing up existing site to $BACKUP_DIR"
    run_cmd "mkdir -p '$BACKUP_DIR'"
    run_cmd "cp -r '$WEB_ROOT' '$BACKUP_DIR/'"
    print_success "✓ Backup created"
fi

# Create web root directory
print_header "Preparing Web Directory"
print_status "Creating web root directory: $WEB_ROOT"
run_cmd "mkdir -p '$WEB_ROOT'"
run_cmd "chown -R www-data:www-data '$WEB_ROOT'" 2>/dev/null || true

# Deploy application files
print_header "Deploying Application"
print_status "Copying build files to $WEB_ROOT"
run_cmd "cp -r '$BUILD_DIR'/* '$WEB_ROOT/'"
run_cmd "chown -R www-data:www-data '$WEB_ROOT'" 2>/dev/null || true
run_cmd "chmod -R 755 '$WEB_ROOT'"
print_success "✓ Application files deployed"

# Generate nginx configuration
print_header "Configuring Nginx"

NGINX_CONFIG="$NGINX_SITES_DIR/$NGINX_SITE"

if [[ "$SSL" == true ]]; then
    SSL_CONFIG="
    listen 443 ssl http2;
    listen [::]:443 ssl http2;
    
    ssl_certificate /etc/ssl/certs/${DOMAIN}.crt;
    ssl_certificate_key /etc/ssl/private/${DOMAIN}.key;
    ssl_protocols TLSv1.2 TLSv1.3;
    ssl_ciphers ECDHE-RSA-AES128-GCM-SHA256:ECDHE-RSA-AES256-GCM-SHA384;
    ssl_prefer_server_ciphers off;
    
    # Redirect HTTP to HTTPS
    if (\$scheme != \"https\") {
        return 301 https://\$server_name\$request_uri;
    }"
    LISTEN_CONFIG="listen 80;"
else
    SSL_CONFIG=""
    LISTEN_CONFIG="listen 80;
    listen [::]:80;"
fi

print_status "Generating nginx configuration"

cat > /tmp/nginx-valency.conf << EOF
server {
    $LISTEN_CONFIG
    $SSL_CONFIG
    
    server_name $DOMAIN;
    root $WEB_ROOT;
    index index.html;

    # Enable gzip compression
    gzip on;
    gzip_vary on;
    gzip_min_length 1024;
    gzip_proxied any;
    gzip_comp_level 6;
    gzip_types
        text/plain
        text/css
        text/xml
        text/javascript
        application/javascript
        application/xml+rss
        application/json
        application/manifest+json
        image/svg+xml;

    # Handle client-side routing for React SPA
    location / {
        try_files \$uri \$uri/ @fallback;
    }

    # Fallback for SPA routing
    location @fallback {
        rewrite ^.*\$ /index.html last;
    }

    # Handle specific SPA routes explicitly
    location ~ ^/(explorer|identification|optimization|chatbot|activities|profile|settings|login|signup|about)(/.*)?$ {
        try_files \$uri \$uri/ /index.html;
    }

    # Cache static JavaScript and CSS files
    location /static/js/ {
        expires 1y;
        add_header Cache-Control "public, immutable";
        add_header Vary "Accept-Encoding";
    }

    location /static/css/ {
        expires 1y;
        add_header Cache-Control "public, immutable";
        add_header Vary "Accept-Encoding";
    }

    # Cache other static assets
    location /static/ {
        expires 1y;
        add_header Cache-Control "public, immutable";
    }

    # Cache images and models
    location /images/ {
        expires 1y;
        add_header Cache-Control "public";
        add_header Vary "Accept-Encoding";
    }

    location /models/ {
        expires 1y;
        add_header Cache-Control "public";
        add_header Vary "Accept-Encoding";
    }

    location /icons/ {
        expires 1y;
        add_header Cache-Control "public";
    }

    # Handle manifest and service worker
    location /manifest.json {
        expires 1d;
        add_header Cache-Control "public";
    }

    location /service-worker.js {
        expires 0;
        add_header Cache-Control "no-cache, no-store, must-revalidate";
    }

    # Security headers
    add_header X-Frame-Options "DENY" always;
    add_header X-XSS-Protection "1; mode=block" always;
    add_header X-Content-Type-Options "nosniff" always;
    add_header Referrer-Policy "strict-origin-when-cross-origin" always;
    add_header Permissions-Policy "geolocation=(), microphone=(), camera=()" always;

    # Optional: API proxy if backend is on different server
    # location /api/ {
    #     proxy_pass http://backend-server:8000/api/;
    #     proxy_set_header Host \$host;
    #     proxy_set_header X-Real-IP \$remote_addr;
    #     proxy_set_header X-Forwarded-For \$proxy_add_x_forwarded_for;
    #     proxy_set_header X-Forwarded-Proto \$scheme;
    # }
}
EOF

run_cmd "cp /tmp/nginx-valency.conf '$NGINX_CONFIG'"
print_success "✓ Nginx configuration created"

# Enable the site
if [[ ! -L "$NGINX_ENABLED_DIR/$NGINX_SITE" ]]; then
    print_status "Enabling nginx site"
    run_cmd "ln -sf '$NGINX_CONFIG' '$NGINX_ENABLED_DIR/$NGINX_SITE'"
    print_success "✓ Site enabled"
else
    print_success "✓ Site already enabled"
fi

# Test nginx configuration
print_status "Testing nginx configuration"
if [[ "$DRY_RUN" == false ]]; then
    if nginx -t; then
        print_success "✓ Nginx configuration test passed"
    else
        print_error "Nginx configuration test failed"
        exit 1
    fi
else
    print_warning "[DRY-RUN] Would test nginx configuration"
fi

# Restart nginx
if [[ "$RESTART_NGINX" == true ]]; then
    print_header "Restarting Nginx"
    run_cmd "systemctl reload nginx"
    print_success "✓ Nginx reloaded"
fi

# Final checks and status
print_header "Deployment Summary"

if [[ "$DRY_RUN" == false ]]; then
    # Check if the site is accessible
    print_status "Verifying deployment..."
    
    if command_exists curl; then
        if curl -s -o /dev/null -w "%{http_code}" "http://$DOMAIN" | grep -q "200"; then
            print_success "✓ Site is accessible at http://$DOMAIN"
        else
            print_warning "⚠ Site may not be accessible - check nginx logs"
        fi
    fi
    
    # Show nginx status
    if systemctl is-active --quiet nginx; then
        print_success "✓ Nginx is running"
    else
        print_error "✗ Nginx is not running"
    fi
    
    echo ""
    print_status "Deployment completed successfully!"
    echo ""
    echo -e "${GREEN}Your Valency app is now deployed at:${NC}"
    echo -e "  ${CYAN}http://$DOMAIN${NC}"
    if [[ "$SSL" == true ]]; then
        echo -e "  ${CYAN}https://$DOMAIN${NC}"
    fi
    echo ""
    echo -e "${BLUE}File locations:${NC}"
    echo -e "  Web root: ${CYAN}$WEB_ROOT${NC}"
    echo -e "  Nginx config: ${CYAN}$NGINX_CONFIG${NC}"
    if [[ "$BACKUP" == true ]] && [[ -n "$BACKUP_DIR" ]]; then
        echo -e "  Backup: ${CYAN}$BACKUP_DIR${NC}"
    fi
    echo ""
    echo -e "${BLUE}Available routes:${NC}"
    echo -e "  / (Home)"
    echo -e "  /explorer (Structure Analysis)" 
    echo -e "  /identification (Identification Tools)"
    echo -e "  /optimization (Optimization Tools)"
    echo -e "  /chatbot (Master Agent)"
    echo -e "  /activities (Activities)"
    echo -e "  /profile (User Profile)"
    echo -e "  /settings (Settings)"
    echo -e "  /about (About Page)"
    
else
    echo ""
    print_status "Dry run completed - no changes were made"
    print_status "Run without --dry-run to perform actual deployment"
fi

echo ""
print_status "Useful commands:"
echo "  Check nginx status: sudo systemctl status nginx"
echo "  View nginx logs: sudo tail -f /var/log/nginx/error.log"
echo "  Test nginx config: sudo nginx -t"
echo "  Reload nginx: sudo systemctl reload nginx"
