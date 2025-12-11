import { handleGenerateCrystal, generateCrystal } from "../src/tools/generation/generate-crystal";
import { handleCreateAlloy } from "../src/tools/transformation/alloy";
import { handleCreateHeterostructure } from "../src/tools/transformation/heterostructure";
import { handleAddAdsorbate } from "../src/tools/transformation/adsorbate";
import { handleApplyStrain } from "../src/tools/transformation/strain";
import { generateNanostructure } from "../src/tools/generation/nanostructure";

// Set timeout to 60s for real Python execution
jest.setTimeout(60000);

describe("MCP Server Integration Tests", () => {

    test("handleGenerateCrystal (GaAs) returns formatted text", async () => {
        const input = {
            composition: ["Ga", "As"],
            space_group: 216,
            lattice_params: { a: 5.653 }
        };

        // Test Handler
        const result = await handleGenerateCrystal(input);

        expect(result.isError).toBeFalsy();
        // Check for Markdown headers or key info
        expect(result.content[0].text).toMatch(/## 💎 Crystal Generated/);
        expect(result.content[0].text).toContain("GaAs");
        expect(result.content[0].text).toContain("Space Group: 216");
    });

    test("handleCreateAlloy (Si -> SiGe)", async () => {
        // 1. Generate Base Si (Logic)
        const baseInput = {
            composition: ["Si"],
            space_group: 227,
            num_atoms: 8
        };
        const baseResult = await generateCrystal(baseInput);
        if (!baseResult.success) throw new Error("Base Si generation failed");
        const structure = baseResult.data.structure;

        // 2. Transmute to Alloy (Handler)
        const alloyInput = {
            structure: structure,
            substitutions: { "Si": { element: "Ge", concentration: 0.5 } },
            seed: 42
        };

        const result = await handleCreateAlloy(alloyInput);

        expect(result.isError).toBeFalsy();
        expect(result.content[0].text).toMatch(/## 🧪 Alloy Created/);
        expect(result.content[0].text).toContain("Si4 Ge4"); // Expect 50/50 split of 8 atoms
    });

    test("handleCreateHeterostructure (Graphene Stack)", async () => {
        // 1. Generate Graphene (Logic)
        const gRes = await generateNanostructure({ type: "graphene", params: { size: [1, 1, 1] } });
        if (!gRes.success) throw new Error("Graphene generation failed");
        const graphene = gRes.data.structure;

        // 2. Stack (Handler)
        const stackInput = {
            substrate: graphene,
            overlayer: graphene,
            interface_distance: 3.4
        };

        const result = await handleCreateHeterostructure(stackInput);

        expect(result.isError).toBeFalsy();
        expect(result.content[0].text).toMatch(/## 🥞 Heterostructure Created/);
        // Check for updated formula or lattice c
        // The text usually contains the structure summary
        expect(result.content[0].text).toMatch(/Formula:.*C4/); // 2+2=4 atoms
    });

    test("handleAddAdsorbate (H on Graphene)", async () => {
        // 1. Generate Graphene (Logic)
        const gRes = await generateNanostructure({ type: "graphene", params: { size: [1, 1, 1] } });
        if (!gRes.success) throw new Error("Graphene generation failed");
        const graphene = gRes.data.structure;

        // 2. Add H (Handler)
        const adsInput = {
            structure: graphene,
            molecule: "H", // Can pass string name? Schema allows string or object. backend supports "H" -> molecule("H")
            site_index: 0,
            distance: 1.5
        };

        // Need to ensure backend `add_adsorbate` handles string molecule name? 
        // `structure_tools.py`: `molecule = atoms_from_dict(molecule_dict)`
        // If I pass "H", `atoms_from_dict("H")` might fail if it expects dict.
        // `structure_utils.py`: `atoms_from_dict` usually expects dict.
        // Wait, `add_adsorbate` implementation in `structure_tools.py` lines 860+:
        // I should check if it handles string.
        // If not, I should create a molecule dict or fix backend to support string.
        // `BaseSchema` says `z.union([z.string(), z.any()])`.
        // I'll check `structure_tools.py`. If it doesn't handle string, I'll pass a dict.

        // For safety in this test, I'll pass a dict.
        const h_mol = {
            lattice: { matrix: [[10, 0, 0], [0, 10, 0], [0, 0, 10]], a: 10, b: 10, c: 10, alpha: 90, beta: 90, gamma: 90, volume: 1000 },
            atoms: [{ element: "H", coords: [0.5, 0.5, 0.5], cartesian: [5, 5, 5] }]
        };

        const result = await handleAddAdsorbate({
            structure: graphene,
            molecule: h_mol,
            site_index: 0,
            distance: 1.5
        });

        expect(result.isError).toBeFalsy();
        expect(result.content[0].text).toMatch(/## 📎 Adsorbate Added/);
        expect(result.content[0].text).toContain("Formula: C2 H1");
    });

    test("handleApplyStrain", async () => {
        // 1. Generate Si
        const baseResult = await generateCrystal({ composition: ["Si"], space_group: 227 });
        if (!baseResult.success) throw new Error("Si gen failed");
        const si = baseResult.data.structure;

        // 2. Apply Strain (Handler)
        const strainInput = {
            structure: si,
            strain_tensor: [0.01, 0, 0, 0, 0.01, 0, 0, 0, 0.01] // 1% hydrostatic
        };

        const result = await handleApplyStrain(strainInput);

        expect(result.isError).toBeFalsy();
        expect(result.content[0].text).toMatch(/## 📏 Strain Applied/);
        // Volume should increase by ~3%
    });

});
