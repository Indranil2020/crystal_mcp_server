/**
 * Diagnostic Test - Traces full pipeline step by step
 * 
 * Run this test to see exactly where the pipeline breaks:
 * npx playwright test tests/e2e/diagnostic.spec.ts --reporter=list
 */

import { test, expect, Page } from '@playwright/test';

test.describe('🔍 DIAGNOSTIC: Step-by-step Pipeline Test', () => {

    test.beforeAll(async () => {
        console.log('\n' + '='.repeat(80));
        console.log('DIAGNOSTIC TEST SUITE');
        console.log('Purpose: Find exact point of failure in the pipeline');
        console.log('='.repeat(80) + '\n');
    });

    test('Step 1: Can reach the web application?', async ({ page }) => {
        console.log('\n[STEP 1] Testing web app accessibility...');

        const response = await page.goto('/', { waitUntil: 'domcontentloaded' });
        console.log(`  - URL: ${response?.url()}`);
        console.log(`  - Status: ${response?.status()}`);

        expect(response?.status()).toBe(200);

        const title = await page.title();
        console.log(`  - Page title: ${title}`);

        console.log('  ✅ PASSED: Web app is accessible');
    });

    test('Step 2: Is the MCP Bridge reachable?', async ({ page }) => {
        console.log('\n[STEP 2] Testing MCP bridge health endpoint...');

        await page.goto('/');

        // Direct fetch to bridge
        const bridgeHealth = await page.evaluate(async () => {
            try {
                const resp = await fetch('http://localhost:8080/health');
                return { ok: resp.ok, status: resp.status, data: await resp.json() };
            } catch (e) {
                return { ok: false, error: String(e) };
            }
        });

        console.log(`  - Bridge response:`, JSON.stringify(bridgeHealth));
        expect(bridgeHealth.ok).toBe(true);

        console.log('  ✅ PASSED: MCP bridge is healthy');
    });

    test('Step 3: Can we initialize MCP?', async ({ page }) => {
        console.log('\n[STEP 3] Testing MCP initialization...');

        await page.goto('/');

        const initResult = await page.evaluate(async () => {
            try {
                const resp = await fetch('http://localhost:8080/mcp/initialize', {
                    method: 'POST',
                    headers: { 'Content-Type': 'application/json' },
                    body: JSON.stringify({
                        protocolVersion: '2024-11-05',
                        capabilities: {},
                        clientInfo: { name: 'diagnostic-test', version: '1.0.0' }
                    })
                });
                return { ok: resp.ok, status: resp.status, data: await resp.json() };
            } catch (e) {
                return { ok: false, error: String(e) };
            }
        });

        console.log(`  - Initialize response:`, JSON.stringify(initResult).slice(0, 500));
        expect(initResult.ok).toBe(true);

        console.log('  ✅ PASSED: MCP initialized');
    });

    test('Step 4: Can we list MCP tools?', async ({ page }) => {
        console.log('\n[STEP 4] Testing MCP tool listing...');

        await page.goto('/');

        // Initialize first
        await page.evaluate(async () => {
            await fetch('http://localhost:8080/mcp/initialize', {
                method: 'POST',
                headers: { 'Content-Type': 'application/json' },
                body: JSON.stringify({
                    protocolVersion: '2024-11-05',
                    capabilities: {},
                    clientInfo: { name: 'diagnostic-test', version: '1.0.0' }
                })
            });
        });

        const toolsResult = await page.evaluate(async () => {
            try {
                const resp = await fetch('http://localhost:8080/mcp/tools');
                const data = await resp.json();
                return {
                    ok: resp.ok,
                    count: data.tools?.length || 0,
                    names: data.tools?.map((t: { name: string }) => t.name) || []
                };
            } catch (e) {
                return { ok: false, error: String(e) };
            }
        });

        console.log(`  - Tools found: ${toolsResult.count}`);
        console.log(`  - Tool names: ${toolsResult.names?.join(', ')}`);
        expect(toolsResult.ok).toBe(true);
        expect(toolsResult.count).toBeGreaterThan(0);

        console.log('  ✅ PASSED: MCP tools listed');
    });

    test('Step 5: Is Ollama reachable?', async ({ page }) => {
        console.log('\n[STEP 5] Testing Ollama connectivity...');

        await page.goto('/');

        const ollamaResult = await page.evaluate(async () => {
            try {
                const resp = await fetch('http://localhost:11434/api/tags');
                const data = await resp.json();
                return {
                    ok: resp.ok,
                    models: data.models?.map((m: { name: string }) => m.name) || []
                };
            } catch (e) {
                return { ok: false, error: String(e) };
            }
        });

        console.log(`  - Ollama response ok: ${ollamaResult.ok}`);
        console.log(`  - Available models: ${ollamaResult.models?.join(', ')}`);
        expect(ollamaResult.ok).toBe(true);

        console.log('  ✅ PASSED: Ollama is reachable');
    });

    test('Step 6: Can Ollama respond to a simple prompt?', async ({ page }) => {
        console.log('\n[STEP 6] Testing basic Ollama chat...');

        await page.goto('/');

        const chatResult = await page.evaluate(async () => {
            try {
                const startTime = Date.now();
                const resp = await fetch('http://localhost:11434/api/chat', {
                    method: 'POST',
                    headers: { 'Content-Type': 'application/json' },
                    body: JSON.stringify({
                        model: 'qwen2.5:7b',
                        messages: [{ role: 'user', content: 'Say hello in one word' }],
                        stream: false
                    })
                });
                const elapsed = Date.now() - startTime;
                const data = await resp.json();
                return {
                    ok: resp.ok,
                    elapsed: elapsed,
                    response: data.message?.content?.slice(0, 100)
                };
            } catch (e) {
                return { ok: false, error: String(e) };
            }
        });

        console.log(`  - Chat response time: ${chatResult.elapsed}ms`);
        console.log(`  - Response: "${chatResult.response}"`);
        expect(chatResult.ok).toBe(true);

        console.log('  ✅ PASSED: Ollama responds to prompts');
    });

    test('Step 7: Can Ollama use tools?', async ({ page }) => {
        console.log('\n[STEP 7] Testing Ollama tool calling...');

        await page.goto('/');

        // First get the tools
        const tools = await page.evaluate(async () => {
            await fetch('http://localhost:8080/mcp/initialize', {
                method: 'POST',
                headers: { 'Content-Type': 'application/json' },
                body: JSON.stringify({
                    protocolVersion: '2024-11-05',
                    capabilities: {},
                    clientInfo: { name: 'diagnostic-test', version: '1.0.0' }
                })
            });
            const resp = await fetch('http://localhost:8080/mcp/tools');
            const data = await resp.json();
            return data.tools || [];
        });

        console.log(`  - Loaded ${tools.length} tools`);

        // Convert to Ollama format
        const ollamaTools = tools.map((t: any) => ({
            type: 'function',
            function: {
                name: t.name,
                description: t.description || '',
                parameters: t.inputSchema || { type: 'object', properties: {} }
            }
        }));

        const toolCallResult = await page.evaluate(async (ollamaTools) => {
            try {
                const startTime = Date.now();
                const resp = await fetch('http://localhost:11434/api/chat', {
                    method: 'POST',
                    headers: { 'Content-Type': 'application/json' },
                    body: JSON.stringify({
                        model: 'qwen2.5:7b',
                        messages: [{ role: 'user', content: 'Generate a benzene molecule' }],
                        stream: false,
                        tools: ollamaTools
                    })
                });
                const elapsed = Date.now() - startTime;
                const data = await resp.json();
                return {
                    ok: resp.ok,
                    elapsed: elapsed,
                    hasToolCalls: !!data.message?.tool_calls,
                    toolCalls: data.message?.tool_calls,
                    textResponse: data.message?.content?.slice(0, 200)
                };
            } catch (e) {
                return { ok: false, error: String(e) };
            }
        }, ollamaTools);

        console.log(`  - Response time: ${toolCallResult.elapsed}ms`);
        console.log(`  - Has tool_calls: ${toolCallResult.hasToolCalls}`);
        if (toolCallResult.hasToolCalls) {
            console.log(`  - Tool calls: ${JSON.stringify(toolCallResult.toolCalls)}`);
        } else {
            console.log(`  - Text response: "${toolCallResult.textResponse}"`);
        }

        // This is where we might see the issue
        if (!toolCallResult.hasToolCalls) {
            console.log('  ⚠️ WARNING: LLM did NOT return tool calls, returned text instead');
            console.log('  This means the LLM is not properly invoking tools.');
        }

        expect(toolCallResult.ok).toBe(true);
        console.log('  ✅ PASSED: Ollama responded (check tool call status above)');
    });

    test('Step 8: Can MCP execute build_molecule directly?', async ({ page }) => {
        console.log('\n[STEP 8] Testing direct MCP tool execution...');

        await page.goto('/');

        // Initialize first
        await page.evaluate(async () => {
            await fetch('http://localhost:8080/mcp/initialize', {
                method: 'POST',
                headers: { 'Content-Type': 'application/json' },
                body: JSON.stringify({
                    protocolVersion: '2024-11-05',
                    capabilities: {},
                    clientInfo: { name: 'diagnostic-test', version: '1.0.0' }
                })
            });
        });

        const callResult = await page.evaluate(async () => {
            try {
                const startTime = Date.now();
                const resp = await fetch('http://localhost:8080/mcp/call', {
                    method: 'POST',
                    headers: { 'Content-Type': 'application/json' },
                    body: JSON.stringify({
                        name: 'build_molecule',
                        arguments: { name: 'benzene' }  // Correct parameter is 'name', not 'molecule_source'
                    })
                });
                const elapsed = Date.now() - startTime;
                const data = await resp.json();

                // Check for structure in response
                let hasStructure = false;
                let atomCount = 0;
                if (data.content) {
                    for (const item of data.content) {
                        if (item.text && item.text.includes('<json-data>')) {
                            hasStructure = true;
                            const match = item.text.match(/<json-data>([\s\S]*?)<\/json-data>/);
                            if (match) {
                                try {
                                    const parsed = JSON.parse(match[1]);
                                    atomCount = parsed.structure?.atoms?.length || parsed.atoms?.length || 0;
                                } catch { }
                            }
                        }
                    }
                }

                return {
                    ok: resp.ok,
                    isError: data.isError,
                    elapsed: elapsed,
                    hasStructure: hasStructure,
                    atomCount: atomCount,
                    contentPreview: JSON.stringify(data).slice(0, 500)
                };
            } catch (e) {
                return { ok: false, error: String(e) };
            }
        });

        console.log(`  - Response time: ${callResult.elapsed}ms`);
        console.log(`  - isError: ${callResult.isError}`);
        console.log(`  - Has structure: ${callResult.hasStructure}`);
        console.log(`  - Atom count: ${callResult.atomCount}`);
        console.log(`  - Content preview: ${callResult.contentPreview}`);

        expect(callResult.ok).toBe(true);
        expect(callResult.hasStructure).toBe(true);
        expect(callResult.atomCount).toBe(12); // C6H6
        console.log('  ✅ PASSED: MCP tool execution works');
    });

    test('Step 9: UI elements present?', async ({ page }) => {
        console.log('\n[STEP 9] Testing UI component presence...');

        // Collect console errors
        const consoleErrors: string[] = [];
        page.on('console', msg => {
            if (msg.type() === 'error') {
                consoleErrors.push(msg.text());
            }
        });

        // Navigate and wait for React to mount
        await page.goto('/', { waitUntil: 'networkidle' });

        // Wait for actual React content to appear (not just the HTML shell)
        // The app renders "Crystal" text once React mounts
        await page.waitForSelector('text=Crystal', { timeout: 15000 }).catch(() => {
            console.log('  ⚠️ Timed out waiting for React content');
        });

        // Additional brief wait for components
        await page.waitForTimeout(1000);

        const uiElements = await page.evaluate(() => {
            const root = document.getElementById('root');
            return {
                title: document.querySelector('h1')?.textContent,
                hasTextarea: !!document.querySelector('textarea'),
                hasSendButton: !!Array.from(document.querySelectorAll('button')).some(b => b.textContent?.includes('Send')),
                hasCanvas: !!document.querySelector('canvas'),
                hasCrystalText: !!document.body.textContent?.includes('Crystal'),
                rootChildCount: root?.childElementCount || 0,
                bodyText: document.body.textContent?.slice(0, 200) || ''
            };
        });

        console.log(`  - Title: ${uiElements.title}`);
        console.log(`  - Root children: ${uiElements.rootChildCount}`);
        console.log(`  - Has Crystal text: ${uiElements.hasCrystalText}`);
        console.log(`  - Has textarea: ${uiElements.hasTextarea}`);
        console.log(`  - Has send button: ${uiElements.hasSendButton}`);
        console.log(`  - Has canvas: ${uiElements.hasCanvas}`);
        console.log(`  - Body text preview: ${uiElements.bodyText.slice(0, 100)}`);

        if (consoleErrors.length > 0) {
            console.log('  ⚠️ Console errors:', consoleErrors.slice(0, 3));
        }

        // Assert React mounted (has any children)
        expect(uiElements.rootChildCount).toBeGreaterThan(0);

        // If UI rendered, check for expected elements
        if (uiElements.rootChildCount > 0) {
            expect(uiElements.hasCrystalText).toBe(true);
            console.log('  ✅ PASSED: UI elements present');
        }
    });

    test.afterAll(async () => {
        console.log('\n' + '='.repeat(80));
        console.log('DIAGNOSTIC SUMMARY');
        console.log('='.repeat(80));
        console.log('\nIf all tests passed, the pipeline is working correctly.');
        console.log('If Step 7 shows "LLM did NOT return tool calls", the issue is in LLM config.');
        console.log('Check browser console for detailed debug logs.\n');
    });
});
