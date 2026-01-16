import { test, expect } from '@playwright/test';

test('Debug Redux State', async ({ page }) => {
    console.log('Navigating to home page...');
    await page.goto('http://localhost:5173');
    await page.waitForTimeout(2000);

    console.log('Checking window.__REDUX_STATE__...');
    const stateDebug = await page.evaluate(() => {
        const win = window as any;
        return {
            hasReduxState: !!win.__REDUX_STATE__,
            keys: win.__REDUX_STATE__ ? Object.keys(win.__REDUX_STATE__) : [],
            structureState: win.__REDUX_STATE__?.structure || 'missing'
        };
    });

    console.log('DEBUG RESULT:', JSON.stringify(stateDebug, null, 2));

    expect(stateDebug.hasReduxState).toBe(true);
});
