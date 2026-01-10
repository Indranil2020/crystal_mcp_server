/**
 * Molecule Suggestion Tool
 *
 * Provides intelligent molecule name suggestions when exact match fails.
 * Uses multiple strategies:
 * - Fuzzy string matching (trigram similarity)
 * - Phonetic matching (sounds-like)
 * - IUPAC fragment recognition
 * - Chemical name normalization (misspellings, abbreviations)
 */

import { SuggestMoleculesSchema } from "../../types/tools.js";
import { Result, createSuccess, createFailure, createError, CrystalErrorCode } from "../../types/errors.js";
import { executePythonWithJSON } from "../../utils/python-bridge.js";

interface SuggestMoleculesInput {
    query: string;
    max_results?: number;
    include_smiles?: boolean;
    min_similarity?: number;
}

interface MoleculeSuggestion {
    name: string;
    smiles?: string;
    score: number;
    match_type: string;
    pubchem_cid?: number;
}

interface SuggestMoleculesResult {
    query: string;
    suggestions: MoleculeSuggestion[];
    normalization: {
        original: string;
        normalized: string;
        was_corrected: boolean;
        variations?: string[];
        misspelling_correction?: string;
        abbreviation_expansion?: string;
    };
    iupac_analysis?: {
        recognized: boolean;
        parent?: string;
        parent_smiles?: string;
        prefixes?: string[];
        suffixes?: string[];
    };
    alternative_queries: string[];
}

export async function suggestMolecules(input: unknown): Promise<Result<SuggestMoleculesResult>> {
    const parsed = SuggestMoleculesSchema.safeParse(input);
    if (!parsed.success) {
        return createFailure(createError(
            CrystalErrorCode.INVALID_INPUT,
            "Invalid input parameters",
            { zodErrors: parsed.error.errors },
            ["Check query parameter"],
            true
        ));
    }

    console.error(`[MCP DEBUG] Suggesting molecules for query: ${parsed.data.query}`);

    // Call Python with the suggest_molecules function
    const result = await executePythonWithJSON<SuggestMoleculesInput, SuggestMoleculesResult>(
        "molecule_suggester.py",
        {
            query: parsed.data.query,
            max_results: parsed.data.max_results ?? 10,
            include_smiles: parsed.data.include_smiles ?? true,
            min_similarity: parsed.data.min_similarity ?? 0.3,
        },
        { timeout: 15000 }
    );

    if (!result.success) {
        console.error(`[MCP DEBUG] Python execution failed: ${result.error.message}`);
        return createFailure(result.error);
    }

    return createSuccess(result.data);
}

export async function handleSuggestMolecules(args: unknown): Promise<any> {
    const result = await suggestMolecules(args);

    if (!result.success) {
        return {
            content: [{
                type: "text",
                text: `**Molecule Suggestion Failed**\n\n${result.error.message}`
            }],
            isError: true
        };
    }

    const data = result.data;
    const suggestions = data.suggestions || [];

    // Build response text
    let responseText = `### Molecule Suggestions for: "${data.query}"\n\n`;

    // Add normalization info
    if (data.normalization?.was_corrected) {
        responseText += `**Normalized**: ${data.normalization.normalized}\n`;
        if (data.normalization.misspelling_correction) {
            responseText += `**Corrected**: ${data.normalization.misspelling_correction}\n`;
        }
        if (data.normalization.abbreviation_expansion) {
            responseText += `**Expanded**: ${data.normalization.abbreviation_expansion}\n`;
        }
        responseText += "\n";
    }

    // Add IUPAC analysis if present
    if (data.iupac_analysis?.recognized) {
        responseText += "**IUPAC Analysis**:\n";
        if (data.iupac_analysis.parent) {
            responseText += `- Parent structure: ${data.iupac_analysis.parent}`;
            if (data.iupac_analysis.parent_smiles) {
                responseText += ` (${data.iupac_analysis.parent_smiles})`;
            }
            responseText += "\n";
        }
        if (data.iupac_analysis.prefixes?.length) {
            responseText += `- Substituents: ${data.iupac_analysis.prefixes.join(", ")}\n`;
        }
        if (data.iupac_analysis.suffixes?.length) {
            responseText += `- Functional groups: ${data.iupac_analysis.suffixes.join(", ")}\n`;
        }
        responseText += "\n";
    }

    // Add suggestions
    if (suggestions.length > 0) {
        responseText += "**Top Matches**:\n\n";
        const topSuggestions = suggestions.slice(0, 5);
        for (let i = 0; i < topSuggestions.length; i++) {
            const s = topSuggestions[i];
            if (!s) continue;
            responseText += `${i + 1}. **${s.name}** (score: ${s.score.toFixed(2)}, ${s.match_type})\n`;
            if (s.smiles) {
                responseText += `   - SMILES: \`${s.smiles}\`\n`;
            }
            if (s.pubchem_cid) {
                responseText += `   - PubChem CID: ${s.pubchem_cid}\n`;
            }
        }
        responseText += "\n";

        // Add hint for using the suggestions
        const firstSuggestion = suggestions[0];
        if (firstSuggestion) {
            responseText += "**To use a suggestion**, call `build_molecule` with:\n";
            responseText += `- The name: \`build_molecule({name: "${firstSuggestion.name}"})\`\n`;
            if (firstSuggestion.smiles) {
                responseText += `- Or the SMILES: \`build_molecule({name: "${firstSuggestion.smiles}", input_type: "smiles"})\`\n`;
            }
        }
    } else {
        responseText += "No similar molecules found in the database.\n\n";
        responseText += "**Suggestions**:\n";
        responseText += "- Try providing the SMILES string directly\n";
        responseText += "- Search PubChem for the correct identifier\n";
        responseText += "- Use the IUPAC systematic name\n";
    }

    // Add alternative queries if any
    if (data.alternative_queries?.length) {
        const validQueries = data.alternative_queries.filter(q => !q.startsWith("Try providing"));
        if (validQueries.length > 0) {
            responseText += "\n**Alternative search terms**: " + validQueries.slice(0, 5).join(", ");
        }
    }

    // Include JSON data for programmatic access
    const jsonData = JSON.stringify({
        success: true,
        ...data
    });

    return {
        content: [{
            type: "text",
            text: responseText
        }, {
            type: "text",
            text: `\n\n<json-data>\n${jsonData}\n</json-data>`
        }]
    };
}
