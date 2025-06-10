#!/bin/bash

# Valency Frontend - Quick Setup and Deployment Script
# This script builds and deploys the React SPA in one go

set -e

# Colors for output
GREEN='\033[0;32m'
BLUE='\033[0;34m'
YELLOW='\033[1;33m'
NC='\033[0m'

echo -e "${BLUE}🚀 Valency Frontend - Quick Setup & Deployment${NC}"
echo -e "${BLUE}===============================================${NC}"
echo ""

# Default values
DOMAIN="localhost"
SSL=false
ANALYZE=false

# Parse arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        --domain)
            DOMAIN="$2"
            shift 2
            ;;
        --ssl)
            SSL=true
            shift
            ;;
        --analyze)
            ANALYZE=true
            shift
            ;;
        --help)
            echo "Usage: $0 [OPTIONS]"
            echo ""
            echo "Options:"
            echo "  --domain DOMAIN    Domain name (default: localhost)"
            echo "  --ssl              Enable SSL/HTTPS"
            echo "  --analyze          Analyze bundle after build"
            echo "  --help             Show this help"
            echo ""
            echo "Example:"
            echo "  $0 --domain myapp.com --ssl --analyze"
            exit 0
            ;;
        *)
            echo "Unknown option: $1"
            echo "Use --help for usage information"
            exit 1
            ;;
    esac
done

echo -e "${YELLOW}Configuration:${NC}"
echo "  Domain: $DOMAIN"
echo "  SSL: $SSL"
echo "  Bundle analysis: $ANALYZE"
echo ""

# Step 1: Build the application
echo -e "${GREEN}Step 1: Building application...${NC}"
BUILD_ARGS="--type production"
if [[ "$ANALYZE" == true ]]; then
    BUILD_ARGS="$BUILD_ARGS --analyze"
fi

if ./build.sh $BUILD_ARGS; then
    echo -e "${GREEN}✓ Build completed successfully${NC}"
else
    echo "❌ Build failed"
    exit 1
fi

echo ""

# Step 2: Deploy to nginx
echo -e "${GREEN}Step 2: Deploying to nginx...${NC}"
DEPLOY_ARGS="--domain $DOMAIN"
if [[ "$SSL" == true ]]; then
    DEPLOY_ARGS="$DEPLOY_ARGS --ssl"
fi

if ./deploy.sh $DEPLOY_ARGS; then
    echo -e "${GREEN}✓ Deployment completed successfully${NC}"
else
    echo "❌ Deployment failed"
    exit 1
fi

echo ""
echo -e "${GREEN}🎉 Valency Frontend is now live!${NC}"
echo ""
echo -e "${BLUE}Access your application at:${NC}"
echo -e "  http://$DOMAIN"
if [[ "$SSL" == true ]]; then
    echo -e "  https://$DOMAIN"
fi
echo ""
echo -e "${BLUE}Key Features Available:${NC}"
echo "  🧬 Structure Analysis (/explorer)"
echo "  🔍 Identification Tools (/identification)"
echo "  ⚗️  Optimization Tools (/optimization)"
echo "  🤖 Master Agent (/chatbot)"
echo "  📊 Activities Dashboard (/activities)"
echo ""
echo -e "${YELLOW}Next Steps:${NC}"
echo "  1. Test all routes to ensure proper SPA routing"
echo "  2. Configure SSL certificates if using --ssl"
echo "  3. Set up your backend API endpoints"
echo "  4. Monitor nginx logs for any issues"
echo ""
echo -e "${BLUE}Useful Commands:${NC}"
echo "  sudo systemctl status nginx    # Check nginx status"
echo "  sudo tail -f /var/log/nginx/error.log  # View logs"
echo "  sudo nginx -t                 # Test configuration"
