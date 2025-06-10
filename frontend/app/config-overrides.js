const path = require('path');

module.exports = function override(config, env) {
  // Enable source maps in production for better debugging
  if (env === 'production') {
    config.devtool = 'source-map';
  }

  // Optimize chunks for better caching and route-based code splitting
  if (env === 'production') {
    config.optimization = {
      ...config.optimization,
      splitChunks: {
        chunks: 'all',
        maxInitialRequests: 6,
        maxAsyncRequests: 8,
        cacheGroups: {
          // Core React and ReactDOM
          react: {
            test: /[\\/]node_modules[\\/](react|react-dom)[\\/]/,
            name: 'react',
            priority: 20,
            chunks: 'all',
          },
          // React Router - critical for SPA routing
          router: {
            test: /[\\/]node_modules[\\/](react-router|react-router-dom)[\\/]/,
            name: 'router',
            priority: 18,
            chunks: 'all',
          },
          // Authentication-related dependencies
          auth: {
            test: /[\\/]node_modules[\\/](axios|crypto-browserify)[\\/]/,
            name: 'auth-libs',
            priority: 16,
            chunks: 'all',
          },
          // UI and Animation libraries
          ui: {
            test: /[\\/]node_modules[\\/](framer-motion|react-icons|@react-three)[\\/]/,
            name: 'ui-libs',
            priority: 14,
            chunks: 'all',
          },
          // Scientific libraries (3D modeling, molecular visualization)
          scientific: {
            test: /[\\/]node_modules[\\/](three|3dmol|@rcsb)[\\/]/,
            name: 'scientific-libs',
            priority: 12,
            chunks: 'all',
          },
          // Other vendor libraries
          vendor: {
            test: /[\\/]node_modules[\\/]/,
            name: 'vendors',
            priority: 10,
            chunks: 'all',
          },
          // Common application code (shared across routes)
          common: {
            name: 'common',
            minChunks: 2,
            priority: 8,
            chunks: 'all',
            enforce: true,
          },
          // Route-specific chunks for lazy loading
          explorer: {
            test: /[\\/]src[\\/]pages[\\/]Explorer[\\/]/,
            name: 'explorer-pages',
            priority: 6,
            chunks: 'async',
            minSize: 20000,
          },
          generator: {
            test: /[\\/]src[\\/]pages[\\/]Generator[\\/]/,
            name: 'generator-pages',
            priority: 6,
            chunks: 'async',
            minSize: 20000,
          },
          discriminator: {
            test: /[\\/]src[\\/]pages[\\/]Discriminator[\\/]/,
            name: 'discriminator-pages',
            priority: 6,
            chunks: 'async',
            minSize: 20000,
          },
        },
      },
      // Runtime chunk for better caching
      runtimeChunk: {
        name: 'runtime',
      },
    };
  }

  // Add alias for better imports
  config.resolve.alias = {
    ...config.resolve.alias,
    '@': path.resolve(__dirname, 'src'),
    '@components': path.resolve(__dirname, 'src/components'),
    '@pages': path.resolve(__dirname, 'src/pages'),
    '@contexts': path.resolve(__dirname, 'src/contexts'),
    '@utils': path.resolve(__dirname, 'src/utils'),
    '@services': path.resolve(__dirname, 'src/services'),
    '@routes': path.resolve(__dirname, 'src/routes'),
    '@contents': path.resolve(__dirname, 'src/contents'),
  };

  // Ensure proper handling of public assets and SPA routing
  config.resolve.fallback = {
    ...config.resolve.fallback,
    "crypto": require.resolve("crypto-browserify"),
    "stream": require.resolve("stream-browserify"),
    "path": require.resolve("path-browserify"),
    "buffer": require.resolve("buffer"),
    "process": require.resolve("process/browser"),
  };

  // Optimize asset handling for SPA
  if (env === 'production') {
    // Ensure proper public path for SPA routing
    config.output.publicPath = '/';
    
    // Optimize asset file names for better caching
    config.output.filename = 'static/js/[name].[contenthash:8].js';
    config.output.chunkFilename = 'static/js/[name].[contenthash:8].chunk.js';
    
    // Configure asset modules for better caching
    config.module.rules.forEach(rule => {
      if (rule.oneOf) {
        rule.oneOf.forEach(oneOf => {
          if (oneOf.test && oneOf.test.toString().includes('\\.(png|jpe?g|gif|webp)')) {
            oneOf.generator = {
              filename: 'static/media/[name].[contenthash:8][ext]'
            };
          }
          if (oneOf.test && oneOf.test.toString().includes('\\.(svg)')) {
            oneOf.generator = {
              filename: 'static/media/[name].[contenthash:8][ext]'
            };
          }
        });
      }
    });
  }

  // Performance optimizations for development
  if (env === 'development') {
    config.optimization = {
      ...config.optimization,
      removeAvailableModules: false,
      removeEmptyChunks: false,
      splitChunks: false,
    };
  }

  return config;
};
