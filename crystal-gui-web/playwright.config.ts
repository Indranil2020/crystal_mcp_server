import { defineConfig, devices } from '@playwright/test';

/**
 * Playwright Configuration for Crystal GUI Web E2E Tests
 * 
 * Uses production build (vite preview) for reliable testing.
 * Dev server has HMR timing issues in headless mode.
 */
export default defineConfig({
    testDir: './tests/e2e',
    fullyParallel: false, // Run sequentially to avoid resource conflicts
    forbidOnly: !!process.env.CI,
    retries: process.env.CI ? 2 : 1,
    workers: 1, // Single worker for LLM/MCP resource sharing
    timeout: 60000, // 60s per test (LLM can be slow)
    reporter: [
        ['html', { outputFolder: 'playwright-report' }],
        ['list'],
    ],

    use: {
        baseURL: 'http://localhost:4173', // Vite preview port
        trace: 'on-first-retry',
        screenshot: 'only-on-failure',
        video: 'on-first-retry',
    },

    projects: [
        {
            name: 'chromium',
            use: {
                ...devices['Desktop Chrome'],
                // Enable WebGL for MolStar 3D rendering in headless mode
                launchOptions: {
                    args: [
                        '--enable-webgl',
                        '--use-gl=swiftshader',  // Software WebGL rendering
                        '--disable-web-security',
                        '--allow-running-insecure-content',
                    ],
                },
            },
        },
    ],

    // Web server configuration - use production build for reliability
    webServer: [
        {
            // Build and serve production bundle
            command: 'npm run build && npm run preview',
            url: 'http://localhost:4173',
            reuseExistingServer: !process.env.CI,
            timeout: 60000, // Allow time for build
        },
        {
            command: 'cd bridge && python3 server.py',
            url: 'http://localhost:8080/health',
            reuseExistingServer: !process.env.CI,
            timeout: 15000,
        },
    ],
});
