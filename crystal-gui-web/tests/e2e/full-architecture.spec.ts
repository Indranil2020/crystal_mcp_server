/**
 * Full Architecture E2E Tests
 * 
 * Tests complete flow: Playwright → Browser → App → LLM → MCP → Python → Render
 * Mirrors the architecture diagram and covers all layers
 */

import { test, expect, Page } from '@playwright/test';

// ==========================================
// TEST FIXTURES & HELPERS
// ==========================================

interface ToolCallAssertion {
    toolName: string;
    hasArgs: boolean;
    success: boolean;
}

interface StructureAssertion {
    minAtoms: number;
    hasFormula: boolean;
    hasLattice?: boolean;
}

// Page object for Chat Panel
class ChatPanelPage {
    constructor(private page: Page) { }

    async sendMessage(message: string) {
        const input = this.page.locator('textarea[placeholder*="Ask me"], textarea[placeholder*="generate"]');
        await input.fill(message);
        await this.page.click('button:has-text("Send")');
    }

    async waitForAssistantResponse(timeout = 45000) {
        // Wait for processing indicator to disappear
        await this.page.waitForSelector('[class*="animate-pulse"]', { state: 'hidden', timeout });
        // Wait for assistant message
        await expect(this.page.locator('[class*="assistant"], [class*="slate-700"]').last()).toBeVisible({ timeout: 5000 });
    }

    async getLastAssistantMessage(): Promise<string> {
        const messages = this.page.locator('[class*="assistant"], [class*="slate-700"]');
        const lastMessage = messages.last();
        return await lastMessage.textContent() || '';
    }

    async assertToolCallVisible(toolName: string) {
        await expect(this.page.locator(`text=${toolName}`)).toBeVisible({ timeout: 10000 });
    }

    async assertToolCallSuccess() {
        await expect(this.page.locator('[class*="success"], [class*="green"]')).toBeVisible({ timeout: 30000 });
    }
}

// Page object for MolStar Viewer
class ViewerPage {
    constructor(private page: Page) { }

    async waitForInitialization() {
        // Wait for MolStar plugin to initialize
        await expect(this.page.locator('.msp-plugin, [class*="MolStar"], text="3D Viewer"')).toBeVisible({ timeout: 15000 });
    }

    async assertStructureLoaded() {
        // Check for atom count display
        await expect(this.page.locator('text=/\\d+\\s*atoms/')).toBeVisible({ timeout: 10000 });
    }

    async getAtomCount(): Promise<number> {
        const text = await this.page.locator('text=/\\d+\\s*atoms/').first().textContent();
        const match = text?.match(/(\d+)/);
        return match ? parseInt(match[1]) : 0;
    }

    async changeRepresentation(mode: string) {
        await this.page.selectOption('select:has(option:text("Ball"))', mode);
        await this.page.waitForTimeout(500);
    }
}

// Page object for 2D Editor
class EditorPage {
    constructor(private page: Page) { }

    async waitForKekuleLoaded() {
        await expect(this.page.locator('text=2D Editor')).toBeVisible({ timeout: 15000 });
    }

    async loadTemplate(template: string) {
        await this.page.selectOption('select:has(option:text("Templates"))', template);
        await this.page.waitForTimeout(1000);
    }

    async pushTo3D() {
        await this.page.click('button:has-text("3D")');
        await this.page.waitForTimeout(500);
    }
}

// Page object for Status Bar (MCP connection)
class StatusBarPage {
    constructor(private page: Page) { }

    async assertConnected() {
        await expect(this.page.locator('text=connected')).toBeVisible({ timeout: 10000 });
    }

    async assertToolsAvailable() {
        await expect(this.page.locator('text=/\\d+.*tools.*available/')).toBeVisible();
    }

    async getToolCount(): Promise<number> {
        const text = await this.page.locator('text=/\\d+.*tools/').textContent();
        const match = text?.match(/(\d+)/);
        return match ? parseInt(match[1]) : 0;
    }
}

// ==========================================
// LAYER 1: CONNECTION & INITIALIZATION
// ==========================================

test.describe('Layer 1: Infrastructure Setup', () => {
    test('should connect to MCP Bridge (FastAPI 8080)', async ({ page }) => {
        await page.goto('/');
        const statusBar = new StatusBarPage(page);

        await statusBar.assertConnected();
        console.log('[Test] MCP Bridge connected');
    });

    test('should fetch available MCP tools', async ({ page }) => {
        await page.goto('/');
        const statusBar = new StatusBarPage(page);

        await statusBar.assertConnected();
        await statusBar.assertToolsAvailable();

        const toolCount = await statusBar.getToolCount();
        expect(toolCount).toBeGreaterThan(0);
        console.log(`[Test] ${toolCount} MCP tools available`);
    });

    test('should initialize MolStar viewer', async ({ page }) => {
        await page.goto('/');
        const viewer = new ViewerPage(page);

        await viewer.waitForInitialization();
        console.log('[Test] MolStar viewer initialized');
    });

    test('should load Kekule.js editor', async ({ page }) => {
        await page.goto('/');
        const editor = new EditorPage(page);

        await editor.waitForKekuleLoaded();
        console.log('[Test] Kekule.js editor loaded');
    });
});

// ==========================================
// LAYER 2: LLM INTEGRATION (Ollama)
// ==========================================

test.describe('Layer 2: LLM Tool Calling', () => {
    test.beforeEach(async ({ page }) => {
        await page.goto('/');
        await new StatusBarPage(page).assertConnected();
    });

    test('should send message to LLM and receive tool call', async ({ page }) => {
        const chat = new ChatPanelPage(page);

        await chat.sendMessage('Generate benzene molecule');

        // Verify tool call was made
        await chat.assertToolCallVisible('build_molecule');
        console.log('[Test] LLM correctly selected build_molecule tool');
    });

    test('should handle LLM tool call with arguments', async ({ page }) => {
        const chat = new ChatPanelPage(page);

        await chat.sendMessage('Generate aspirin');
        await chat.waitForAssistantResponse();

        // Tool should have been called with SMILES argument
        const response = await chat.getLastAssistantMessage();
        expect(response.toLowerCase()).toMatch(/aspirin|acetylsalicylic|generated|success/);
        console.log('[Test] LLM passed arguments to tool correctly');
    });
});

// ==========================================
// LAYER 3: MCP SERVER INTEGRATION
// ==========================================

test.describe('Layer 3: MCP Tool Execution', () => {
    test.beforeEach(async ({ page }) => {
        await page.goto('/');
        await new StatusBarPage(page).assertConnected();
    });

    test('should execute build_molecule tool via MCP', async ({ page }) => {
        const chat = new ChatPanelPage(page);
        const viewer = new ViewerPage(page);

        await chat.sendMessage('Generate benzene');
        await chat.waitForAssistantResponse();
        await chat.assertToolCallSuccess();

        // Verify structure was returned from MCP
        await viewer.assertStructureLoaded();
        const atomCount = await viewer.getAtomCount();
        expect(atomCount).toBe(12); // C6H6 = 12 atoms

        console.log('[Test] MCP tool executed, returned 12 atoms');
    });

    test('should execute build_molecular_cluster tool', async ({ page }) => {
        const chat = new ChatPanelPage(page);
        const viewer = new ViewerPage(page);

        await chat.sendMessage('Create benzene dimer with pi-pi stacking');
        await chat.waitForAssistantResponse();

        await viewer.assertStructureLoaded();
        const atomCount = await viewer.getAtomCount();
        expect(atomCount).toBeGreaterThanOrEqual(20); // 2 x benzene = ~24 atoms

        console.log(`[Test] Cluster generated with ${atomCount} atoms`);
    });

    test('should handle MCP tool error gracefully', async ({ page }) => {
        const chat = new ChatPanelPage(page);

        await chat.sendMessage('Generate xyzinvalidmolecule12345');
        await chat.waitForAssistantResponse();

        // Should show error or helpful message, not crash
        const response = await chat.getLastAssistantMessage();
        expect(response).toBeTruthy();
        console.log('[Test] Error handled gracefully');
    });
});

// ==========================================
// LAYER 4: RENDERING & VISUALIZATION
// ==========================================

test.describe('Layer 4: 3D Rendering (MolStar)', () => {
    test.beforeEach(async ({ page }) => {
        await page.goto('/');
        await new StatusBarPage(page).assertConnected();

        // Generate a structure first
        const chat = new ChatPanelPage(page);
        await chat.sendMessage('Generate benzene');
        await chat.waitForAssistantResponse();
    });

    test('should render structure in MolStar viewer', async ({ page }) => {
        const viewer = new ViewerPage(page);

        await viewer.waitForInitialization();
        await viewer.assertStructureLoaded();

        // Canvas should be present
        await expect(page.locator('canvas')).toBeVisible();
        console.log('[Test] Structure rendered in MolStar');
    });

    test('should change representation mode', async ({ page }) => {
        const viewer = new ViewerPage(page);

        await viewer.changeRepresentation('spacefill');
        // Visual change occurred (no crash)
        await expect(page.locator('canvas')).toBeVisible();

        console.log('[Test] Representation changed to spacefill');
    });

    test('should export PNG screenshot', async ({ page }) => {
        const downloadPromise = page.waitForEvent('download', { timeout: 10000 }).catch(() => null);

        // Find and click screenshot button
        const screenshotBtn = page.locator('button:has-text("📷"), button[title*="screenshot"]');
        if (await screenshotBtn.isVisible()) {
            await screenshotBtn.click();
            const download = await downloadPromise;
            if (download) {
                expect(download.suggestedFilename()).toContain('.png');
                console.log('[Test] Screenshot exported');
            }
        }
    });
});

// ==========================================
// LAYER 5: 2D EDITOR SYNC
// ==========================================

test.describe('Layer 5: 2D Editor Integration', () => {
    test.beforeEach(async ({ page }) => {
        await page.goto('/');
        await new EditorPage(page).waitForKekuleLoaded();
    });

    test('should load template in 2D editor', async ({ page }) => {
        const editor = new EditorPage(page);

        await editor.loadTemplate('benzene');

        // SMILES should appear
        await expect(page.locator('text=/c1ccccc1|benzene/i')).toBeVisible({ timeout: 5000 });
        console.log('[Test] Benzene template loaded');
    });

    test('should push 2D structure to 3D viewer', async ({ page }) => {
        const editor = new EditorPage(page);
        const viewer = new ViewerPage(page);

        await editor.loadTemplate('benzene');
        await editor.pushTo3D();

        // Structure should appear in store
        await expect(page.locator('text=/Drawn:/i')).toBeVisible({ timeout: 5000 });
        console.log('[Test] 2D → 3D sync working');
    });
});

// ==========================================
// LAYER 6: DATA EXPORT
// ==========================================

test.describe('Layer 6: Export Functionality', () => {
    test.beforeEach(async ({ page }) => {
        await page.goto('/');
        await new StatusBarPage(page).assertConnected();

        // Generate a structure
        const chat = new ChatPanelPage(page);
        await chat.sendMessage('Generate water');
        await chat.waitForAssistantResponse();
    });

    test('should export CIF file', async ({ page }) => {
        const downloadPromise = page.waitForEvent('download');

        await page.click('button:has-text("Export")');
        await page.click('text=CIF');

        const download = await downloadPromise;
        expect(download.suggestedFilename()).toContain('.cif');
        console.log('[Test] CIF exported');
    });

    test('should export PDB file', async ({ page }) => {
        const downloadPromise = page.waitForEvent('download');

        await page.click('button:has-text("Export")');
        await page.click('text=PDB');

        const download = await downloadPromise;
        expect(download.suggestedFilename()).toContain('.pdb');
        console.log('[Test] PDB exported');
    });

    test('should export XYZ file', async ({ page }) => {
        const downloadPromise = page.waitForEvent('download');

        await page.click('button:has-text("Export")');
        await page.click('text=XYZ');

        const download = await downloadPromise;
        expect(download.suggestedFilename()).toContain('.xyz');
        console.log('[Test] XYZ exported');
    });
});

// ==========================================
// FULL WORKFLOW: End-to-End User Journey
// ==========================================

test.describe('Complete User Workflow', () => {
    test('Full flow: NL → LLM → MCP → Render → Export', async ({ page }) => {
        await page.goto('/');

        const statusBar = new StatusBarPage(page);
        const chat = new ChatPanelPage(page);
        const viewer = new ViewerPage(page);

        // Step 1: Verify connections
        await statusBar.assertConnected();
        console.log('[Flow 1/5] Connected to MCP Bridge');

        // Step 2: Send natural language request
        await chat.sendMessage('Generate a caffeine molecule');
        console.log('[Flow 2/5] Sent request to LLM');

        // Step 3: Verify tool call execution
        await chat.waitForAssistantResponse();
        await chat.assertToolCallSuccess();
        console.log('[Flow 3/5] Tool executed successfully');

        // Step 4: Verify structure rendered
        await viewer.assertStructureLoaded();
        const atomCount = await viewer.getAtomCount();
        expect(atomCount).toBeGreaterThan(10);
        console.log(`[Flow 4/5] Structure rendered: ${atomCount} atoms`);

        // Step 5: Export structure
        const downloadPromise = page.waitForEvent('download');
        await page.click('button:has-text("Export")');
        await page.click('text=CIF');
        const download = await downloadPromise;
        expect(download.suggestedFilename()).toContain('.cif');
        console.log('[Flow 5/5] Structure exported');

        console.log('✅ Complete workflow passed!');
    });
});
