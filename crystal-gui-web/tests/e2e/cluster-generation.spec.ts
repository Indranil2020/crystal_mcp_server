/**
 * Cluster Generation E2E Tests
 * 
 * Tests molecular cluster arrangements: dimers, pi-stacking, H-bonding
 */

import { test, expect, Page } from '@playwright/test';

const CLUSTERS = [
    {
        name: 'benzene dimer',
        prompt: 'Create benzene dimer with pi-pi stacking',
        expectedMolecules: 2,
        minAtoms: 20,
    },
    {
        name: 'water tetramer',
        prompt: 'Build water tetramer in circular arrangement',
        expectedMolecules: 4,
        minAtoms: 10,
    },
];

async function sendChatMessage(page: Page, message: string) {
    const input = page.locator('textarea[placeholder*="Ask me"]');
    await input.fill(message);
    await page.click('button:has-text("Send")');
}

async function waitForToolCompletion(page: Page, timeout = 45000) {
    await expect(page.locator('[class*="animate-pulse"]')).toBeHidden({ timeout });
}

test.describe('Cluster Generation', () => {
    test.beforeEach(async ({ page }) => {
        await page.goto('/');
        await expect(page.locator('[class*="bg-green"]')).toBeVisible({ timeout: 10000 });
    });

    for (const cluster of CLUSTERS) {
        test(`should generate ${cluster.name}`, async ({ page }) => {
            await sendChatMessage(page, cluster.prompt);
            await waitForToolCompletion(page);

            // Check structure was created
            await expect(page.locator('[class*="success"], text=/generated|created/i')).toBeVisible();

            // Verify atom count
            const atomText = await page.locator('text=/\\d+\\s*atoms/').first().textContent();
            const atomMatch = atomText?.match(/(\d+)/);
            const atomCount = atomMatch ? parseInt(atomMatch[1]) : 0;

            expect(atomCount).toBeGreaterThanOrEqual(cluster.minAtoms);
        });
    }

    test('should show cluster metadata', async ({ page }) => {
        await sendChatMessage(page, 'Create benzene dimer with pi stacking');
        await waitForToolCompletion(page);

        // Should show n_molecules or cluster info
        const content = await page.content();
        expect(
            content.includes('molecules') ||
            content.includes('dimer') ||
            content.includes('cluster')
        ).toBe(true);
    });
});

test.describe('Arrangement Types', () => {
    test.beforeEach(async ({ page }) => {
        await page.goto('/');
        await expect(page.locator('[class*="bg-green"]')).toBeVisible({ timeout: 10000 });
    });

    test('should create parallel displaced arrangement', async ({ page }) => {
        await sendChatMessage(page, 'Create benzene dimer parallel displaced 3.5 angstrom');
        await waitForToolCompletion(page);

        await expect(page.locator('text=/success|generated|created/i')).toBeVisible();
    });

    test('should create t-shaped arrangement', async ({ page }) => {
        await sendChatMessage(page, 'Create benzene dimer t-shaped');
        await waitForToolCompletion(page);

        await expect(page.locator('text=/success|generated|created/i')).toBeVisible();
    });

    test('should create hydrogen-bonded cluster', async ({ page }) => {
        await sendChatMessage(page, 'Create water dimer with hydrogen bond');
        await waitForToolCompletion(page);

        await expect(page.locator('text=/success|generated|created/i')).toBeVisible();
    });
});
