#!/bin/bash

# Build script for Valency React SPA
# Handles different build scenarios for optimal deployment

set -e

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

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

# Default values
BUILD_TYPE="production"
ANALYZE=false
CLEAN=false
DEPLOY_TARGET=""

# Help function
show_help() {
    echo "Usage: $0 [OPTIONS]"
    echo ""
    echo "Options:"
    echo "  -t, --type TYPE        Build type: development, production (default: production)"
    echo "  -a, --analyze         Analyze bundle after build"
    echo "  -c, --clean           Clean cache before build"
    echo "  -d, --deploy TARGET   Deploy after build (netlify, vercel, docker)"
    echo "  -h, --help            Show this help message"
    echo ""
    echo "Examples:"
    echo "  $0 --type production --analyze"
    echo "  $0 --clean --deploy netlify"
    echo "  $0 --type development"
}

# Parse arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        -t|--type)
            BUILD_TYPE="$2"
            shift 2
            ;;
        -a|--analyze)
            ANALYZE=true
            shift
            ;;
        -c|--clean)
            CLEAN=true
            shift
            ;;
        -d|--deploy)
            DEPLOY_TARGET="$2"
            shift 2
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

# Validate build type
if [[ "$BUILD_TYPE" != "development" && "$BUILD_TYPE" != "production" ]]; then
    print_error "Invalid build type: $BUILD_TYPE"
    print_error "Supported types: development, production"
    exit 1
fi

print_status "Starting build process for Valency React SPA"
print_status "Build type: $BUILD_TYPE"

# Clean cache if requested
if [[ "$CLEAN" == true ]]; then
    print_status "Cleaning cache and build artifacts..."
    npm run clean 2>/dev/null || {
        rm -rf build node_modules/.cache 2>/dev/null || true
    }
    print_success "Cache cleaned"
fi

# Check if node_modules exists
if [[ ! -d "node_modules" ]]; then
    print_status "Installing dependencies..."
    npm ci
    print_success "Dependencies installed"
fi

# Set environment based on build type
if [[ "$BUILD_TYPE" == "production" ]]; then
    export NODE_ENV=production
    BUILD_COMMAND="build:prod"
else
    export NODE_ENV=development
    BUILD_COMMAND="build"
fi

# Run the build
print_status "Building application..."
if npm run $BUILD_COMMAND; then
    print_success "Build completed successfully"
else
    print_error "Build failed"
    exit 1
fi

# Analyze bundle if requested
if [[ "$ANALYZE" == true ]]; then
    print_status "Analyzing bundle..."
    npm run build:analyze
fi

# Deploy if requested
if [[ -n "$DEPLOY_TARGET" ]]; then
    print_status "Deploying to $DEPLOY_TARGET..."
    case $DEPLOY_TARGET in
        netlify)
            npm run deploy:netlify
            ;;
        vercel)
            npm run deploy:vercel
            ;;
        docker)
            print_status "Building Docker image..."
            docker build -t valency-frontend .
            print_success "Docker image built: valency-frontend"
            ;;
        *)
            print_error "Unknown deploy target: $DEPLOY_TARGET"
            print_error "Supported targets: netlify, vercel, docker"
            exit 1
            ;;
    esac
    print_success "Deployment completed"
fi

print_success "Build process completed successfully!"

# Show build size information
if [[ -d "build" ]]; then
    print_status "Build size information:"
    du -sh build/
    echo ""
    print_status "Static assets:"
    find build/static -name "*.js" -o -name "*.css" | head -10 | while read file; do
        size=$(du -h "$file" | cut -f1)
        echo "  $size - $(basename "$file")"
    done
    echo ""
    
    # Validate SPA build for nginx deployment
    print_status "Validating SPA build for nginx deployment..."
    
    # Check critical files
    if [[ -f "build/index.html" ]]; then
        print_success "✓ index.html found"
    else
        print_error "✗ index.html missing"
        exit 1
    fi
    
    # Check for router configuration
    if grep -q "react-router" build/static/js/*.js 2>/dev/null; then
        print_success "✓ React Router detected"
    else
        print_warning "⚠ React Router not detected - check routing configuration"
    fi
    
    # Check for proper asset hashing
    if ls build/static/js/*.*.js 1> /dev/null 2>&1; then
        print_success "✓ JavaScript files have content hash"
    else
        print_warning "⚠ JavaScript files may not have content hash"
    fi
    
    # Create deployment info
    BUILD_TIME=$(date -u +"%Y-%m-%dT%H:%M:%SZ")
    GIT_HASH=$(git rev-parse --short HEAD 2>/dev/null || echo "unknown")
    GIT_BRANCH=$(git branch --show-current 2>/dev/null || echo "unknown")
    
    cat > build/deployment-info.json << EOF
{
  "buildTime": "$BUILD_TIME",
  "buildType": "$BUILD_TYPE",
  "gitHash": "$GIT_HASH",
  "gitBranch": "$GIT_BRANCH",
  "nodeVersion": "$(node --version)",
  "routes": [
    "/",
    "/explorer",
    "/identification", 
    "/optimization",
    "/chatbot",
    "/activities",
    "/profile",
    "/settings",
    "/login",
    "/signup",
    "/about"
  ]
}
EOF
    
    print_success "✓ Deployment info created"
    
    # Copy nginx config for reference
    if [[ -f "nginx.conf" ]]; then
        cp nginx.conf build/nginx.conf.example
        print_success "✓ Nginx configuration copied to build directory"
    fi
    
    echo ""
    print_status "Build ready for nginx deployment!"
    print_status "Next step: Run ./deploy.sh to deploy to your web server"
fi
