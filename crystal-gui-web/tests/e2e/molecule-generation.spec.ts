/**
 * Molecule Generation E2E Tests
 * 
 * Tests the complete workflow: NL input → LLM → MCP → render
 * Mirrors the Python end_to_end_test.py but through browser automation
 */

import { test, expect, Page } from '@playwright/test';

// Test fixtures
const MOLECULES = [
    { name: 'benzene', prompt: 'Generate benzene', expectedAtoms: 12, formula: 'C6H6' },
    { name: 'water', prompt: 'Generate water molecule', expectedAtoms: 3, formula: 'H2O' },
    { name: 'methane', prompt: 'Generate methane', expectedAtoms: 5, formula: 'CH4' },
];

// Page object helpers
async function sendChatMessage(page: Page, message: string) {
    const input = page.locator('textarea[placeholder*="Ask me"]');
    await input.fill(message);
    await page.click('button:has-text("Send")');
}

async function waitForToolCompletion(page: Page, timeout = 30000) {
    // Wait for processing indicator to disappear
    await expect(page.locator('[class*="animate-pulse"]')).toBeHidden({ timeout });
    // Wait for tool success indicator
    await expect(page.locator('.tool-call.success, [class*="text-green"]')).toBeVisible({ timeout });
}

async function getStructureInfo(page: Page) {
    const infoText = await page.locator('[class*="atoms"]').textContent();
    const atomMatch = infoText?.match(/(\d+)\s*atoms/);
    return {
        atomCount: atomMatch ? parseInt(atomMatch[1]) : 0,
    };
}

// Tests
test.describe('Molecule Generation', () => {
    test.beforeEach(async ({ page }) => {
        await page.goto('/');
        // Wait for connection
        await expect(page.locator('[class*="bg-green"]')).toBeVisible({ timeout: 10000 });
    });

    test('should display welcome message', async ({ page }) => {
        await expect(page.getByText('Crystal GUI')).toBeVisible();
        await expect(page.getByText('Welcome')).toBeVisible();
    });

    test('should show connected status', async ({ page }) => {
        await expect(page.locator('text=connected')).toBeVisible();
        await expect(page.locator('text=tools available')).toBeVisible();
    });

    for (const molecule of MOLECULES) {
        test(`should generate ${molecule.name}`, async ({ page }) => {
            // Send generation request
            await sendChatMessage(page, molecule.prompt);

            // Wait for tool execution
            await waitForToolCompletion(page);

            // Verify structure loaded
            await expect(page.locator(`text=${molecule.formula}`)).toBeVisible();

            // Check atom count in info
            const info = await getStructureInfo(page);
            expect(info.atomCount).toBe(molecule.expectedAtoms);

            // Verify 3D viewer shows structure
            await expect(page.locator('.msp-plugin, [class*="MolStar"]')).toBeVisible();
        });
    }

    test('should handle invalid molecule gracefully', async ({ page }) => {
        await sendChatMessage(page, 'Generate xyzabc123invalid');

        // Wait for response
        await page.waitForTimeout(5000);

        // Should show error or helpful message
        const hasError = await page.locator('text=/error|failed|unable|unknown/i').isVisible();
        const hasResponse = await page.locator('[class*="assistant"]').isVisible();

        expect(hasError || hasResponse).toBe(true);
    });
});

test.describe('Tool Visualization', () => {
    test.beforeEach(async ({ page }) => {
        await page.goto('/');
        await expect(page.locator('[class*="bg-green"]')).toBeVisible({ timeout: 10000 });
    });

    test('should show loading state during tool execution', async ({ page }) => {
        await sendChatMessage(page, 'Generate benzene');

        // Should show loading indicator
        await expect(page.locator('[class*="animate"], [class*="loading"]')).toBeVisible();
    });

    test('should show tool call card with details', async ({ page }) => {
        await sendChatMessage(page, 'Generate benzene');

        // Wait for tool card to appear
        await expect(page.locator('.tool-call, [class*="tool"]')).toBeVisible({ timeout: 10000 });

        // Should show tool name
        await expect(page.locator('text=build_molecule')).toBeVisible();
    });
});

test.describe('Export Functionality', () => {
    test.beforeEach(async ({ page }) => {
        await page.goto('/');
        await expect(page.locator('[class*="bg-green"]')).toBeVisible({ timeout: 10000 });

        // Generate a structure first
        await sendChatMessage(page, 'Generate benzene');
        await waitForToolCompletion(page);
    });

    test('should enable export button when structure loaded', async ({ page }) => {
        const exportBtn = page.locator('button:has-text("Export")');
        await expect(exportBtn).toBeEnabled();
    });

    test('should show export menu with format options', async ({ page }) => {
        await page.click('button:has-text("Export")');

        await expect(page.locator('text=CIF')).toBeVisible();
        await expect(page.locator('text=PDB')).toBeVisible();
        await expect(page.locator('text=XYZ')).toBeVisible();
    });

    test('should download CIF file', async ({ page }) => {
        const downloadPromise = page.waitForEvent('download');

        await page.click('button:has-text("Export")');
        await page.click('text=CIF');

        const download = await downloadPromise;
        expect(download.suggestedFilename()).toContain('.cif');
    });
});
