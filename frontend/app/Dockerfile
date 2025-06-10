# Multi-stage build for production-ready React SPA
# Stage 1: Build the React app
FROM node:20-alpine AS builder

# Set working directory
WORKDIR /app

# Install dependencies for native modules
RUN apk add --no-cache python3 make g++

# Copy package files
COPY package*.json ./

# Install dependencies with optimizations
RUN npm ci --only=production --silent && \
    npm cache clean --force

# Copy source code
COPY . .

# Build the application with production optimizations
RUN npm run build:prod

# Stage 2: Production server with Nginx
FROM nginx:alpine AS production

# Install additional tools for better performance
RUN apk add --no-cache \
    curl \
    brotli

# Remove default nginx config
RUN rm /etc/nginx/conf.d/default.conf

# Copy custom nginx configuration
COPY nginx.conf /etc/nginx/conf.d/default.conf

# Copy built app from builder stage
COPY --from=builder /app/build /usr/share/nginx/html

# Add nginx user and set permissions
RUN chown -R nginx:nginx /usr/share/nginx/html && \
    chmod -R 755 /usr/share/nginx/html

# Create directory for nginx logs
RUN mkdir -p /var/log/nginx

# Expose port
EXPOSE 80

# Health check for container
HEALTHCHECK --interval=30s --timeout=3s --start-period=5s --retries=3 \
    CMD curl -f http://localhost/ || exit 1

# Start nginx
CMD ["nginx", "-g", "daemon off;"]
