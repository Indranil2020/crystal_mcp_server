/**
 * 2D Editor E2E Tests
 * 
 * Tests Kekule.js editor integration and push-to-3D workflow
 */

import { test, expect, Page } from '@playwright/test';

async function waitForKekuleLoaded(page: Page) {
    await expect(page.locator('[class*="Kekule"], text=2D Editor')).toBeVisible({ timeout: 15000 });
}

test.describe('2D Editor', () => {
    test.beforeEach(async ({ page }) => {
        await page.goto('/');
        await waitForKekuleLoaded(page);
    });

    test('should load Kekule editor', async ({ page }) => {
        await expect(page.locator('text=2D Editor')).toBeVisible();
    });

    test('should show template dropdown', async ({ page }) => {
        const templateSelect = page.locator('select:has(option:text("Templates"))');
        await expect(templateSelect).toBeVisible();
    });

    test('should load benzene template', async ({ page }) => {
        // Select benzene template
        await page.selectOption('select:has(option:text("Benzene"))', 'benzene');

        // Wait for template to load
        await page.waitForTimeout(1000);

        // SMILES should appear
        await expect(page.locator('text=/c1ccccc1|benzene/i')).toBeVisible();
    });

    test('should push structure to 3D viewer', async ({ page }) => {
        // Load a template
        await page.selectOption('select:has(option:text("Benzene"))', 'benzene');
        await page.waitForTimeout(500);

        // Click push to 3D
        await page.click('button:has-text("3D")');

        // Should create a structure
        await expect(page.locator('text=/Drawn:|structure/i')).toBeVisible({ timeout: 5000 });
    });

    test('should clear editor', async ({ page }) => {
        // Load template first
        await page.selectOption('select:has(option:text("Benzene"))', 'benzene');
        await page.waitForTimeout(500);

        // Clear
        await page.click('button:has-text("Clear")');

        // Wait for clear
        await page.waitForTimeout(500);
    });
});
