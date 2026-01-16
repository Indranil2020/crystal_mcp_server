"""
molecule/molecule_database.py - Universal Molecule Database

SQLite-based molecule database supporting:
- Local cache of user lookups
- ChEMBL (~2.4M bioactive molecules)
- PubChemLite (~500K curated)
- PubChem Complete (~130M, optional 15GB download)
- ZINC (purchasable compounds)
- DrugBank (approved drugs)

Enhanced with:
- Multi-strategy fuzzy matching (trigram, phonetic, prefix)
- IUPAC fragment recognition and synthesis
- Chemical name normalization
- Common misspelling corrections

The database automatically incorporates any downloaded databases
placed in the data directory.
"""

import sqlite3
import os
import gzip
import csv
import logging
import hashlib
import re
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Any, Union, Set
from dataclasses import dataclass, field
from datetime import datetime
import importlib.util

# Check if RDKit is available for InChIKey generation
RDKIT_AVAILABLE = importlib.util.find_spec("rdkit") is not None

if RDKIT_AVAILABLE:
    from rdkit import Chem
    from rdkit.Chem.inchi import MolFromInchi, InchiToInchiKey

logger = logging.getLogger(__name__)


# =============================================================================
# Chemical Name Utilities - Normalization, IUPAC Fragments, Phonetics
# =============================================================================

class ChemicalNameNormalizer:
    """
    Normalizes chemical names for better matching.
    Handles common variations, abbreviations, and formatting differences.
    """

    # Common abbreviations and their expansions
    ABBREVIATIONS = {
        'ptcda': 'perylene-3,4,9,10-tetracarboxylic dianhydride',
        'ntcda': 'naphthalene-1,4,5,8-tetracarboxylic dianhydride',
        'pmdi': 'pyromellitic diimide',
        'ptcdi': 'perylene-3,4,9,10-tetracarboxylic diimide',
        'tcnq': 'tetracyanoquinodimethane',
        'ttf': 'tetrathiafulvalene',
        'dna': 'deoxyribonucleic acid',
        'rna': 'ribonucleic acid',
        'atp': 'adenosine triphosphate',
        'adp': 'adenosine diphosphate',
        'amp': 'adenosine monophosphate',
        'nad': 'nicotinamide adenine dinucleotide',
        'nadp': 'nicotinamide adenine dinucleotide phosphate',
        'fad': 'flavin adenine dinucleotide',
        'coq10': 'coenzyme q10',
        'thf': 'tetrahydrofuran',
        'dmf': 'dimethylformamide',
        'dmso': 'dimethyl sulfoxide',
        'nmp': 'n-methyl-2-pyrrolidone',
        'dma': 'dimethylacetamide',
        'dcm': 'dichloromethane',
        'thc': 'tetrahydrocannabinol',
        'cbd': 'cannabidiol',
        'lsd': 'lysergic acid diethylamide',
        'mdma': '3,4-methylenedioxymethamphetamine',
        'peg': 'polyethylene glycol',
        'pva': 'polyvinyl alcohol',
        'pvp': 'polyvinylpyrrolidone',
        'edta': 'ethylenediaminetetraacetic acid',
        'tris': 'tris(hydroxymethyl)aminomethane',
        'hepes': '4-(2-hydroxyethyl)-1-piperazineethanesulfonic acid',
    }

    # Common misspellings
    MISSPELLINGS = {
        'asprin': 'aspirin',
        'cafeine': 'caffeine',
        'caffiene': 'caffeine',
        'acetomenophen': 'acetaminophen',
        'acetomeniphen': 'acetaminophen',
        'paracetamol': 'acetaminophen',  # Alternative name
        'ibuprofen': 'ibuprofen',  # Just for safety
        'ibuprofin': 'ibuprofen',
        'naproxan': 'naproxen',
        'penicillan': 'penicillin',
        'penicillin': 'penicillin',
        'methadrone': 'methadone',
        'morphene': 'morphine',
        'cocain': 'cocaine',
        'herion': 'heroin',
        'fentanyl': 'fentanyl',
        'fentanil': 'fentanyl',
        'methamphetamine': 'methamphetamine',
        'amphetamine': 'amphetamine',
        'amphetamin': 'amphetamine',
        'benzine': 'benzene',
        'toluol': 'toluene',
        'xylol': 'xylene',
        'napthalene': 'naphthalene',
        'napthalen': 'naphthalene',
        'antrasen': 'anthracene',
        'antracene': 'anthracene',
        'flourene': 'fluorene',
        'flourine': 'fluorine',
        'clorine': 'chlorine',
        'bromine': 'bromine',
        'iodene': 'iodine',
        'sulfer': 'sulfur',
        'sulfur': 'sulfur',
        'sulphur': 'sulfur',
        'phosphorus': 'phosphorus',
        'phosphorous': 'phosphorus',
        'nitrogen': 'nitrogen',
        'oxigen': 'oxygen',
        'hidrogen': 'hydrogen',
        'cabon': 'carbon',
        'silicon': 'silicon',
        'silicone': 'silicon',  # Common confusion
        'glycerin': 'glycerol',
        'glycerine': 'glycerol',
        'formaldehide': 'formaldehyde',
        'acetaldehide': 'acetaldehyde',
        'acetaldeyde': 'acetaldehyde',
        'metanol': 'methanol',
        'ethanol': 'ethanol',
        'etanol': 'ethanol',
        'propanol': 'propanol',
        'butanol': 'butanol',
        'acetone': 'acetone',
        'aceton': 'acetone',
        'toluen': 'toluene',
        'chloroform': 'chloroform',
        'cloroform': 'chloroform',
        'phenol': 'phenol',
        'fenol': 'phenol',
        'aniline': 'aniline',
        'anilin': 'aniline',
        'pyridine': 'pyridine',
        'pyridin': 'pyridine',
        'pyrole': 'pyrrole',
        'pyrrole': 'pyrrole',
        'furane': 'furan',
        'thiophene': 'thiophene',
        'thiophen': 'thiophene',
        # Typo aliases
        'peryne': 'perylene',  # Common typo
        # Formula aliases for common PAHs (Polycyclic Aromatic Hydrocarbons)
        'c6h6': 'benzene',
        'c10h8': 'naphthalene',
        'c14h10': 'anthracene',
        'c16h10': 'pyrene',
        'c18h12': 'triphenylene',  # or chrysene
        'c20h12': 'perylene',
        'c24h12': 'coronene',
        'c8h8': 'cubane',  # or styrene
        'c20h20': 'dodecahedrane',
        'c60': 'fullerene',
    }

    # Greek letter mappings (for IUPAC names)
    GREEK_LETTERS = {
        'alpha': 'α', 'α': 'alpha',
        'beta': 'β', 'β': 'beta',
        'gamma': 'γ', 'γ': 'gamma',
        'delta': 'δ', 'δ': 'delta',
        'epsilon': 'ε', 'ε': 'epsilon',
        'zeta': 'ζ', 'ζ': 'zeta',
        'eta': 'η', 'η': 'eta',
        'theta': 'θ', 'θ': 'theta',
        'iota': 'ι', 'ι': 'iota',
        'kappa': 'κ', 'κ': 'kappa',
        'lambda': 'λ', 'λ': 'lambda',
        'mu': 'μ', 'μ': 'mu',
        'nu': 'ν', 'ν': 'nu',
        'xi': 'ξ', 'ξ': 'xi',
        'omicron': 'ο', 'ο': 'omicron',
        'pi': 'π', 'π': 'pi',
        'rho': 'ρ', 'ρ': 'rho',
        'sigma': 'σ', 'σ': 'sigma',
        'tau': 'τ', 'τ': 'tau',
        'upsilon': 'υ', 'υ': 'upsilon',
        'phi': 'φ', 'φ': 'phi',
        'chi': 'χ', 'χ': 'chi',
        'psi': 'ψ', 'ψ': 'psi',
        'omega': 'ω', 'ω': 'omega',
    }

    @classmethod
    def normalize(cls, name: str) -> str:
        """
        Normalize a chemical name for matching.

        Steps:
        1. Convert to lowercase
        2. Expand abbreviations
        3. Fix common misspellings
        4. Normalize Greek letters
        5. Normalize whitespace and punctuation
        """
        if not name:
            return ""

        normalized = name.lower().strip()

        # Check abbreviations first (exact match)
        if normalized in cls.ABBREVIATIONS:
            normalized = cls.ABBREVIATIONS[normalized]

        # Fix common misspellings
        if normalized in cls.MISSPELLINGS:
            normalized = cls.MISSPELLINGS[normalized]

        # Note: Greek letter conversion is disabled to avoid false positives
        # (e.g., 'aspirin' contains 'pi' but shouldn't become 'asπrin')
        # Greek letters in chemical names are rare and usually explicit (α, β, γ)

        # Normalize punctuation variations
        normalized = normalized.replace('_', '-')
        normalized = normalized.replace(' - ', '-')
        normalized = re.sub(r'\s+', ' ', normalized)  # Collapse whitespace

        return normalized

    @classmethod
    def get_variations(cls, name: str) -> List[str]:
        """
        Generate possible variations of a chemical name.

        Returns list of names to try in order of likelihood.
        """
        variations = [name.lower().strip()]
        normalized = cls.normalize(name)

        if normalized not in variations:
            variations.append(normalized)

        # Try without hyphens
        no_hyphen = normalized.replace('-', '')
        if no_hyphen not in variations:
            variations.append(no_hyphen)

        # Try with spaces instead of hyphens
        with_spaces = normalized.replace('-', ' ')
        if with_spaces not in variations:
            variations.append(with_spaces)

        # Check for abbreviation expansion
        if name.lower() in cls.ABBREVIATIONS:
            expanded = cls.ABBREVIATIONS[name.lower()]
            if expanded not in variations:
                variations.append(expanded)

        # Check for misspelling correction
        if name.lower() in cls.MISSPELLINGS:
            corrected = cls.MISSPELLINGS[name.lower()]
            if corrected not in variations:
                variations.append(corrected)

        return variations


class IUPACFragmentParser:
    """
    Parses and recognizes IUPAC name fragments for smarter matching.

    Can identify:
    - Substituent prefixes (methyl-, ethyl-, chloro-, etc.)
    - Parent chain names (benzene, cyclohexane, etc.)
    - Functional group suffixes (-ol, -one, -oic acid, etc.)
    - Position indicators (1,2-, ortho-, meta-, para-)
    """

    # Common IUPAC prefixes with their meanings
    PREFIXES = {
        # Alkyl groups
        'methyl': {'type': 'alkyl', 'carbon': 1, 'smiles_fragment': 'C'},
        'ethyl': {'type': 'alkyl', 'carbon': 2, 'smiles_fragment': 'CC'},
        'propyl': {'type': 'alkyl', 'carbon': 3, 'smiles_fragment': 'CCC'},
        'isopropyl': {'type': 'alkyl', 'carbon': 3, 'smiles_fragment': 'C(C)C'},
        'butyl': {'type': 'alkyl', 'carbon': 4, 'smiles_fragment': 'CCCC'},
        'tert-butyl': {'type': 'alkyl', 'carbon': 4, 'smiles_fragment': 'C(C)(C)C'},
        'pentyl': {'type': 'alkyl', 'carbon': 5, 'smiles_fragment': 'CCCCC'},
        'hexyl': {'type': 'alkyl', 'carbon': 6, 'smiles_fragment': 'CCCCCC'},
        'phenyl': {'type': 'aryl', 'carbon': 6, 'smiles_fragment': 'c1ccccc1'},
        'benzyl': {'type': 'aryl', 'carbon': 7, 'smiles_fragment': 'Cc1ccccc1'},
        'vinyl': {'type': 'alkenyl', 'carbon': 2, 'smiles_fragment': 'C=C'},
        'allyl': {'type': 'alkenyl', 'carbon': 3, 'smiles_fragment': 'CC=C'},

        # Halogen substituents
        'fluoro': {'type': 'halogen', 'element': 'F', 'smiles_fragment': 'F'},
        'chloro': {'type': 'halogen', 'element': 'Cl', 'smiles_fragment': 'Cl'},
        'bromo': {'type': 'halogen', 'element': 'Br', 'smiles_fragment': 'Br'},
        'iodo': {'type': 'halogen', 'element': 'I', 'smiles_fragment': 'I'},

        # Functional groups as prefixes
        'amino': {'type': 'functional', 'group': 'amine', 'smiles_fragment': 'N'},
        'hydroxy': {'type': 'functional', 'group': 'alcohol', 'smiles_fragment': 'O'},
        'hydroxyl': {'type': 'functional', 'group': 'alcohol', 'smiles_fragment': 'O'},
        'nitro': {'type': 'functional', 'group': 'nitro', 'smiles_fragment': '[N+](=O)[O-]'},
        'cyano': {'type': 'functional', 'group': 'nitrile', 'smiles_fragment': 'C#N'},
        'carboxy': {'type': 'functional', 'group': 'carboxylic acid', 'smiles_fragment': 'C(=O)O'},
        'acetyl': {'type': 'functional', 'group': 'acyl', 'smiles_fragment': 'C(=O)C'},
        'formyl': {'type': 'functional', 'group': 'aldehyde', 'smiles_fragment': 'C=O'},
        'oxo': {'type': 'functional', 'group': 'ketone', 'smiles_fragment': '=O'},
        'thio': {'type': 'functional', 'group': 'thiol', 'smiles_fragment': 'S'},
        'mercapto': {'type': 'functional', 'group': 'thiol', 'smiles_fragment': 'S'},
        'sulfo': {'type': 'functional', 'group': 'sulfonic acid', 'smiles_fragment': 'S(=O)(=O)O'},
        'phospho': {'type': 'functional', 'group': 'phosphate', 'smiles_fragment': 'P(=O)(O)O'},
    }

    # Common parent structures
    PARENT_CHAINS = {
        # Alkanes
        'methane': {'type': 'alkane', 'carbon': 1, 'smiles': 'C'},
        'ethane': {'type': 'alkane', 'carbon': 2, 'smiles': 'CC'},
        'propane': {'type': 'alkane', 'carbon': 3, 'smiles': 'CCC'},
        'butane': {'type': 'alkane', 'carbon': 4, 'smiles': 'CCCC'},
        'pentane': {'type': 'alkane', 'carbon': 5, 'smiles': 'CCCCC'},
        'hexane': {'type': 'alkane', 'carbon': 6, 'smiles': 'CCCCCC'},
        'heptane': {'type': 'alkane', 'carbon': 7, 'smiles': 'CCCCCCC'},
        'octane': {'type': 'alkane', 'carbon': 8, 'smiles': 'CCCCCCCC'},
        'nonane': {'type': 'alkane', 'carbon': 9, 'smiles': 'CCCCCCCCC'},
        'decane': {'type': 'alkane', 'carbon': 10, 'smiles': 'CCCCCCCCCC'},

        # Cycloalkanes
        'cyclopropane': {'type': 'cycloalkane', 'carbon': 3, 'smiles': 'C1CC1'},
        'cyclobutane': {'type': 'cycloalkane', 'carbon': 4, 'smiles': 'C1CCC1'},
        'cyclopentane': {'type': 'cycloalkane', 'carbon': 5, 'smiles': 'C1CCCC1'},
        'cyclohexane': {'type': 'cycloalkane', 'carbon': 6, 'smiles': 'C1CCCCC1'},

        # Aromatics
        'benzene': {'type': 'aromatic', 'carbon': 6, 'smiles': 'c1ccccc1'},
        'toluene': {'type': 'aromatic', 'carbon': 7, 'smiles': 'Cc1ccccc1'},
        'naphthalene': {'type': 'aromatic', 'carbon': 10, 'smiles': 'c1ccc2ccccc2c1'},
        'anthracene': {'type': 'aromatic', 'carbon': 14, 'smiles': 'c1ccc2cc3ccccc3cc2c1'},
        'phenanthrene': {'type': 'aromatic', 'carbon': 14, 'smiles': 'c1ccc2c(c1)ccc1ccccc12'},
        'pyrene': {'type': 'aromatic', 'carbon': 16, 'smiles': 'c1cc2ccc3cccc4ccc(c1)c2c34'},
        'perylene': {'type': 'aromatic', 'carbon': 20, 'smiles': 'c1cc2ccc3ccc4cccc5ccc(c1)c2c3c45'},
        'coronene': {'type': 'aromatic', 'carbon': 24, 'smiles': 'c1cc2ccc3ccc4ccc5ccc6ccc1c7c2c3c4c5c67'},

        # Heterocycles
        'furan': {'type': 'heterocycle', 'carbon': 4, 'smiles': 'c1ccoc1'},
        'pyrrole': {'type': 'heterocycle', 'carbon': 4, 'smiles': 'c1cc[nH]c1'},
        'thiophene': {'type': 'heterocycle', 'carbon': 4, 'smiles': 'c1ccsc1'},
        'pyridine': {'type': 'heterocycle', 'carbon': 5, 'smiles': 'c1ccncc1'},
        'pyrimidine': {'type': 'heterocycle', 'carbon': 4, 'smiles': 'c1cncnc1'},
        'imidazole': {'type': 'heterocycle', 'carbon': 3, 'smiles': 'c1cnc[nH]1'},
        'oxazole': {'type': 'heterocycle', 'carbon': 3, 'smiles': 'c1cocn1'},
        'thiazole': {'type': 'heterocycle', 'carbon': 3, 'smiles': 'c1cscn1'},
        'indole': {'type': 'heterocycle', 'carbon': 8, 'smiles': 'c1ccc2[nH]ccc2c1'},
        'quinoline': {'type': 'heterocycle', 'carbon': 9, 'smiles': 'c1ccc2ncccc2c1'},
        'isoquinoline': {'type': 'heterocycle', 'carbon': 9, 'smiles': 'c1ccc2ccncc2c1'},
        'purine': {'type': 'heterocycle', 'carbon': 5, 'smiles': 'c1ncc2[nH]cnc2n1'},
    }

    # Functional group suffixes
    SUFFIXES = {
        '-ol': {'type': 'alcohol', 'smiles_mod': 'O'},
        '-diol': {'type': 'diol', 'smiles_mod': 'O...O'},
        '-triol': {'type': 'triol', 'smiles_mod': 'O...O...O'},
        '-one': {'type': 'ketone', 'smiles_mod': '=O'},
        '-al': {'type': 'aldehyde', 'smiles_mod': 'C=O'},
        '-oic acid': {'type': 'carboxylic acid', 'smiles_mod': 'C(=O)O'},
        '-carboxylic acid': {'type': 'carboxylic acid', 'smiles_mod': 'C(=O)O'},
        '-amine': {'type': 'amine', 'smiles_mod': 'N'},
        '-imine': {'type': 'imine', 'smiles_mod': '=N'},
        '-amide': {'type': 'amide', 'smiles_mod': 'C(=O)N'},
        '-nitrile': {'type': 'nitrile', 'smiles_mod': 'C#N'},
        '-thiol': {'type': 'thiol', 'smiles_mod': 'S'},
        '-ene': {'type': 'alkene', 'smiles_mod': '='},
        '-yne': {'type': 'alkyne', 'smiles_mod': '#'},
        '-anhydride': {'type': 'anhydride', 'smiles_mod': 'C(=O)OC(=O)'},
        '-imide': {'type': 'imide', 'smiles_mod': 'C(=O)NC(=O)'},
    }

    # Position indicators
    POSITION_PREFIXES = {
        'ortho': {'positions': [1, 2], 'abbrev': 'o'},
        'o-': {'positions': [1, 2], 'abbrev': 'o'},
        'meta': {'positions': [1, 3], 'abbrev': 'm'},
        'm-': {'positions': [1, 3], 'abbrev': 'm'},
        'para': {'positions': [1, 4], 'abbrev': 'p'},
        'p-': {'positions': [1, 4], 'abbrev': 'p'},
        'cis': {'geometry': 'cis'},
        'trans': {'geometry': 'trans'},
        'r-': {'stereochem': 'R'},
        's-': {'stereochem': 'S'},
        'd-': {'stereochem': 'D'},
        'l-': {'stereochem': 'L'},
    }

    @classmethod
    def parse(cls, name: str) -> Dict[str, Any]:
        """
        Parse an IUPAC name into its component parts.

        Returns:
            Dict with 'prefixes', 'parent', 'suffixes', 'positions'
        """
        name_lower = name.lower().strip()
        result = {
            'original': name,
            'prefixes': [],
            'parent': None,
            'suffixes': [],
            'positions': [],
            'recognized': False,
        }

        # Find parent chain
        for parent_name, info in cls.PARENT_CHAINS.items():
            if parent_name in name_lower:
                result['parent'] = {'name': parent_name, **info}
                result['recognized'] = True
                break

        # Find prefixes
        for prefix_name, info in cls.PREFIXES.items():
            if prefix_name in name_lower:
                result['prefixes'].append({'name': prefix_name, **info})
                result['recognized'] = True

        # Find suffixes
        for suffix, info in cls.SUFFIXES.items():
            if name_lower.endswith(suffix) or suffix.replace('-', '') in name_lower:
                result['suffixes'].append({'name': suffix, **info})
                result['recognized'] = True

        # Find position indicators
        for pos, info in cls.POSITION_PREFIXES.items():
            if name_lower.startswith(pos) or f'-{pos}' in name_lower:
                result['positions'].append({'indicator': pos, **info})

        # Extract numeric positions (e.g., "1,2-" or "3,4,9,10-")
        position_pattern = r'(\d+(?:,\d+)*)-'
        matches = re.findall(position_pattern, name_lower)
        for match in matches:
            positions = [int(p) for p in match.split(',')]
            result['positions'].append({'numeric': positions})

        return result

    @classmethod
    def suggest_related_names(cls, parsed: Dict[str, Any]) -> List[str]:
        """
        Given a parsed IUPAC name, suggest related molecule names that might match.

        E.g., for "1,4-phenylene-2,3-dicarboximide":
        - Try "phthalimide" (common name for phenyl imide)
        - Try variations with "benzene" instead of "phenylene"
        """
        suggestions = []

        # If we found prefixes and suffixes, try combining them differently
        if parsed.get('prefixes'):
            for prefix_info in parsed['prefixes']:
                prefix = prefix_info['name']
                # Suggest parent + prefix combinations
                for parent_name in cls.PARENT_CHAINS.keys():
                    suggestions.append(f"{prefix}{parent_name}")
                    suggestions.append(f"{prefix} {parent_name}")

        if parsed.get('parent'):
            parent_name = parsed['parent']['name']

            # If parent is found, suggest common derivatives
            if parent_name == 'benzene':
                suggestions.extend(['phenol', 'aniline', 'toluene', 'styrene', 'benzoic acid'])
            elif parent_name == 'naphthalene':
                suggestions.extend(['naphthol', 'naphthylamine', 'naphthoic acid'])
            elif parent_name == 'perylene':
                suggestions.extend(['ptcda', 'ptcdi', 'perylene-3,4,9,10-tetracarboxylic dianhydride'])

        # If we see "imide" or "anhydride" in suffixes
        for suffix_info in parsed.get('suffixes', []):
            if suffix_info['type'] == 'imide':
                suggestions.extend(['phthalimide', 'succinimide', 'maleimide', 'naphthalimide'])
            elif suffix_info['type'] == 'anhydride':
                suggestions.extend(['phthalic anhydride', 'maleic anhydride', 'succinic anhydride'])

        return list(set(suggestions))


class PhoneticMatcher:
    """
    Implements phonetic matching for chemical names using Double Metaphone.

    This helps match names that sound similar but are spelled differently,
    e.g., "asprin" → "aspirin", "caffiene" → "caffeine"
    """

    # Vowels
    VOWELS = 'AEIOU'

    @classmethod
    def double_metaphone(cls, word: str) -> Tuple[str, str]:
        """
        Generate Double Metaphone codes for a word.

        Returns two codes: (primary, secondary)
        The secondary may be empty if there's only one pronunciation.

        This is a simplified implementation focused on chemical names.
        """
        if not word:
            return ('', '')

        word = word.upper()
        word = re.sub(r'[^A-Z]', '', word)  # Remove non-letters

        if not word:
            return ('', '')

        # Common chemical name transformations
        primary = ''
        secondary = ''

        # Skip leading vowels for code, but keep them for matching
        i = 0
        if word[0] in cls.VOWELS:
            primary += word[0]
            i = 1

        while i < len(word):
            char = word[i]
            prev_char = word[i-1] if i > 0 else ''
            next_char = word[i+1] if i < len(word) - 1 else ''
            next_next = word[i+2] if i < len(word) - 2 else ''

            # Consonant transformations
            if char == 'B':
                if not (i == len(word) - 1 and prev_char == 'M'):
                    primary += 'P'
            elif char == 'C':
                if next_char in 'EIY':
                    primary += 'S'
                elif next_char == 'H':
                    primary += 'X'  # CH -> X
                    i += 1
                else:
                    primary += 'K'
            elif char == 'D':
                if next_char == 'G' and next_next in 'EIY':
                    primary += 'J'
                    i += 1
                else:
                    primary += 'T'
            elif char == 'F':
                primary += 'F'
            elif char == 'G':
                if next_char in 'EIY':
                    primary += 'J'
                elif next_char == 'H':
                    # GH is usually silent
                    i += 1
                else:
                    primary += 'K'
            elif char == 'H':
                if prev_char in cls.VOWELS and next_char in cls.VOWELS:
                    pass  # Silent H between vowels
                elif i == 0 or prev_char in cls.VOWELS:
                    primary += 'H'
            elif char == 'J':
                primary += 'J'
            elif char == 'K':
                if prev_char != 'C':
                    primary += 'K'
            elif char == 'L':
                primary += 'L'
            elif char == 'M':
                primary += 'M'
            elif char == 'N':
                primary += 'N'
            elif char == 'P':
                if next_char == 'H':
                    primary += 'F'
                    i += 1
                else:
                    primary += 'P'
            elif char == 'Q':
                primary += 'K'
            elif char == 'R':
                primary += 'R'
            elif char == 'S':
                if next_char == 'H':
                    primary += 'X'
                    i += 1
                elif next_char == 'I' and next_next in 'AO':
                    primary += 'X'
                else:
                    primary += 'S'
            elif char == 'T':
                if next_char == 'H':
                    primary += '0'  # TH
                    i += 1
                elif next_char == 'I' and next_next in 'AO':
                    primary += 'X'
                else:
                    primary += 'T'
            elif char == 'V':
                primary += 'F'
            elif char == 'W':
                if next_char in cls.VOWELS:
                    primary += 'W'
            elif char == 'X':
                primary += 'KS'
            elif char == 'Y':
                if next_char in cls.VOWELS:
                    primary += 'Y'
            elif char == 'Z':
                primary += 'S'
            elif char in cls.VOWELS:
                if i == 0:
                    primary += char
                # Internal vowels are generally skipped in metaphone

            i += 1

        return (primary, secondary if secondary else primary)

    @classmethod
    def match_score(cls, word1: str, word2: str) -> float:
        """
        Calculate phonetic similarity score between two words.

        Returns:
            Score from 0.0 (no match) to 1.0 (perfect match)
        """
        code1_pri, code1_sec = cls.double_metaphone(word1)
        code2_pri, code2_sec = cls.double_metaphone(word2)

        if not code1_pri or not code2_pri:
            return 0.0

        # Exact primary match
        if code1_pri == code2_pri:
            return 1.0

        # One primary matches other's secondary
        if code1_pri == code2_sec or code1_sec == code2_pri:
            return 0.9

        # Calculate similarity between codes
        def code_similarity(c1: str, c2: str) -> float:
            if not c1 or not c2:
                return 0.0
            common = sum(1 for a, b in zip(c1, c2) if a == b)
            return common / max(len(c1), len(c2))

        best_sim = max(
            code_similarity(code1_pri, code2_pri),
            code_similarity(code1_pri, code2_sec),
            code_similarity(code1_sec, code2_pri),
            code_similarity(code1_sec, code2_sec),
        )

        return best_sim * 0.8  # Scale down non-exact matches


@dataclass
class FuzzySuggestion:
    """A molecule suggestion with match details."""
    name: str
    smiles: str
    score: float
    match_type: str
    source: str = ""
    pubchem_cid: Optional[int] = None
    details: Dict[str, Any] = field(default_factory=dict)

    def to_dict(self) -> Dict[str, Any]:
        return {
            "name": self.name,
            "smiles": self.smiles,
            "score": self.score,
            "match_type": self.match_type,
            "source": self.source,
            "pubchem_cid": self.pubchem_cid,
            **self.details
        }


@dataclass
class MoleculeRecord:
    """Represents a molecule in the database."""
    id: int
    name: str
    smiles: str
    inchi: Optional[str]
    inchikey: Optional[str]
    formula: Optional[str]
    molecular_weight: Optional[float]
    source: str
    category: Optional[str]
    pubchem_cid: Optional[int]
    chembl_id: Optional[str]
    created_at: Optional[str]
    
    def to_dict(self) -> Dict[str, Any]:
        return {
            "id": self.id,
            "name": self.name,
            "smiles": self.smiles,
            "inchi": self.inchi,
            "inchikey": self.inchikey,
            "formula": self.formula,
            "molecular_weight": self.molecular_weight,
            "source": self.source,
            "category": self.category,
            "pubchem_cid": self.pubchem_cid,
            "chembl_id": self.chembl_id,
        }


class MoleculeDatabase:
    """
    Unified interface to all molecule databases.
    
    Provides fast SQLite-based lookup for millions of molecules.
    Automatically discovers and integrates downloaded databases.
    """
    
    # Database schema version for migrations
    SCHEMA_VERSION = 1
    
    # Default paths
    DEFAULT_DATA_DIR = Path(__file__).parent.parent.parent.parent / "data" / "molecule"
    DEFAULT_DB_NAME = "molecules.db"
    
    def __init__(self, db_path: Optional[Union[str, Path]] = None, data_dir: Optional[Union[str, Path]] = None):
        """
        Initialize molecule database.
        
        Args:
            db_path: Path to SQLite database file. If None, uses default location.
            data_dir: Directory for database files and downloads.
        """
        self.data_dir = Path(data_dir) if data_dir else self.DEFAULT_DATA_DIR
        self.db_path = Path(db_path) if db_path else self.data_dir / self.DEFAULT_DB_NAME
        
        # Ensure data directory exists
        self.data_dir.mkdir(parents=True, exist_ok=True)
        
        # Initialize database
        self._conn: Optional[sqlite3.Connection] = None
        self._ensure_database()
    
    @property
    def conn(self) -> sqlite3.Connection:
        """Get database connection, creating if needed."""
        if self._conn is None:
            self._conn = sqlite3.connect(str(self.db_path))
            self._conn.row_factory = sqlite3.Row
        return self._conn
    
    def _ensure_database(self) -> None:
        """Ensure database exists with correct schema."""
        cursor = self.conn.cursor()
        
        # Check if molecules table exists
        cursor.execute("""
            SELECT name FROM sqlite_master 
            WHERE type='table' AND name='molecules'
        """)
        
        if cursor.fetchone() is None:
            self._create_schema()
        
        # Check for and scan downloadable databases
        self._scan_downloaded_databases()
    
    def _create_schema(self) -> None:
        """Create database schema."""
        cursor = self.conn.cursor()
        
        # Main molecules table
        cursor.execute("""
            CREATE TABLE IF NOT EXISTS molecules (
                id INTEGER PRIMARY KEY AUTOINCREMENT,
                name TEXT,
                smiles TEXT NOT NULL,
                inchi TEXT,
                inchikey TEXT,
                formula TEXT,
                molecular_weight REAL,
                source TEXT NOT NULL DEFAULT 'user',
                category TEXT,
                pubchem_cid INTEGER,
                chembl_id TEXT,
                created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP
            )
        """)
        
        # Aliases table for name lookup
        cursor.execute("""
            CREATE TABLE IF NOT EXISTS aliases (
                id INTEGER PRIMARY KEY AUTOINCREMENT,
                alias TEXT NOT NULL COLLATE NOCASE,
                molecule_id INTEGER NOT NULL,
                FOREIGN KEY (molecule_id) REFERENCES molecules(id) ON DELETE CASCADE,
                UNIQUE(alias COLLATE NOCASE)
            )
        """)
        
        # Sources table - tracks imported databases
        cursor.execute("""
            CREATE TABLE IF NOT EXISTS sources (
                id INTEGER PRIMARY KEY AUTOINCREMENT,
                name TEXT UNIQUE NOT NULL,
                version TEXT,
                molecule_count INTEGER DEFAULT 0,
                file_path TEXT,
                imported_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP
            )
        """)
        
        # Schema metadata
        cursor.execute("""
            CREATE TABLE IF NOT EXISTS schema_meta (
                key TEXT PRIMARY KEY,
                value TEXT
            )
        """)
        
        # Create indexes for fast lookup
        cursor.execute("CREATE INDEX IF NOT EXISTS idx_molecules_name ON molecules(name COLLATE NOCASE)")
        cursor.execute("CREATE INDEX IF NOT EXISTS idx_molecules_smiles ON molecules(smiles)")
        cursor.execute("CREATE INDEX IF NOT EXISTS idx_molecules_inchikey ON molecules(inchikey)")
        cursor.execute("CREATE INDEX IF NOT EXISTS idx_molecules_formula ON molecules(formula)")
        cursor.execute("CREATE INDEX IF NOT EXISTS idx_molecules_source ON molecules(source)")
        cursor.execute("CREATE INDEX IF NOT EXISTS idx_molecules_pubchem_cid ON molecules(pubchem_cid)")
        cursor.execute("CREATE INDEX IF NOT EXISTS idx_molecules_chembl_id ON molecules(chembl_id)")
        cursor.execute("CREATE INDEX IF NOT EXISTS idx_aliases_alias ON aliases(alias COLLATE NOCASE)")
        
        # Store schema version
        cursor.execute(
            "INSERT OR REPLACE INTO schema_meta (key, value) VALUES (?, ?)",
            ("schema_version", str(self.SCHEMA_VERSION))
        )
        
        self.conn.commit()
        logger.info(f"Created molecule database schema at {self.db_path}")
    
    def _scan_downloaded_databases(self) -> None:
        """Scan data directory for downloadable databases and integrate them."""
        downloads_dir = self.data_dir / "downloads"
        if not downloads_dir.exists():
            downloads_dir.mkdir(parents=True, exist_ok=True)
            return
        
        # Look for known database files
        database_patterns = {
            "chembl_*.db": self._attach_chembl_database,
            "pubchem_cid_smiles*.gz": self._import_pubchem_smiles,
            "pubchemlite*.csv": self._import_pubchemlite,
        }
        
        for pattern, importer in database_patterns.items():
            for db_file in downloads_dir.glob(pattern):
                if not self._is_source_imported(db_file.name):
                    logger.info(f"Found new database: {db_file.name}")
                    # Import will be triggered on first lookup or explicit call
    
    def _is_source_imported(self, source_name: str) -> bool:
        """Check if a source has been imported."""
        cursor = self.conn.cursor()
        cursor.execute("SELECT 1 FROM sources WHERE name = ?", (source_name,))
        return cursor.fetchone() is not None
    
    def _register_source(self, name: str, version: str, count: int, file_path: str) -> None:
        """Register an imported source."""
        cursor = self.conn.cursor()
        cursor.execute("""
            INSERT OR REPLACE INTO sources (name, version, molecule_count, file_path, imported_at)
            VALUES (?, ?, ?, ?, ?)
        """, (name, version, count, file_path, datetime.now().isoformat()))
        self.conn.commit()
    
    # =========================================================================
    # Lookup Methods
    # =========================================================================
    
    def lookup(
        self, 
        identifier: str, 
        id_type: str = "auto"
    ) -> Optional[MoleculeRecord]:
        """
        Look up molecule by any identifier.
        
        Args:
            identifier: Name, SMILES, InChI, InChIKey, CID, or ChEMBL ID
            id_type: Type hint ("auto", "name", "smiles", "inchi", "inchikey", 
                     "cid", "chembl")
        
        Returns:
            MoleculeRecord if found, None otherwise
        """
        identifier = identifier.strip()
        
        if id_type == "auto":
            id_type = self._detect_identifier_type(identifier)
        
        lookup_methods = {
            "name": self._lookup_by_name,
            "smiles": self._lookup_by_smiles,
            "inchi": self._lookup_by_inchi,
            "inchikey": self._lookup_by_inchikey,
            "cid": self._lookup_by_cid,
            "chembl": self._lookup_by_chembl,
        }
        
        method = lookup_methods.get(id_type, self._lookup_by_name)
        return method(identifier)
    
    def _detect_identifier_type(self, identifier: str) -> str:
        """Detect the type of identifier."""
        # InChI
        if identifier.startswith("InChI="):
            return "inchi"
        
        # InChIKey (27 characters, AAA-BBB-C format)
        if len(identifier) == 27 and identifier.count("-") == 2:
            return "inchikey"
        
        # PubChem CID (pure numeric)
        if identifier.isdigit():
            return "cid"
        
        # ChEMBL ID
        if identifier.upper().startswith("CHEMBL"):
            return "chembl"
        
        # SMILES (contains special chars like =, #, @, [, ], (, ), /)
        smiles_chars = set("=#@[]()/%\\")
        if any(c in identifier for c in smiles_chars):
            return "smiles"
        
        # Default to name
        return "name"
    
    def _lookup_by_name(self, name: str) -> Optional[MoleculeRecord]:
        """Look up by name or alias."""
        cursor = self.conn.cursor()
        
        # Try direct name match
        cursor.execute("""
            SELECT * FROM molecules WHERE name = ? COLLATE NOCASE LIMIT 1
        """, (name,))
        row = cursor.fetchone()
        
        if row is None:
            # Try alias lookup
            cursor.execute("""
                SELECT m.* FROM molecules m
                JOIN aliases a ON m.id = a.molecule_id
                WHERE a.alias = ? COLLATE NOCASE
                LIMIT 1
            """, (name,))
            row = cursor.fetchone()
        
        if row:
            return self._row_to_record(row)
        return None
    
    def _lookup_by_smiles(self, smiles: str) -> Optional[MoleculeRecord]:
        """Look up by SMILES string."""
        cursor = self.conn.cursor()
        
        # Try exact match first
        cursor.execute("""
            SELECT * FROM molecules WHERE smiles = ? LIMIT 1
        """, (smiles,))
        row = cursor.fetchone()
        
        if row is None and RDKIT_AVAILABLE:
            # Try canonicalized SMILES
            mol = Chem.MolFromSmiles(smiles)
            if mol:
                canonical = Chem.MolToSmiles(mol, canonical=True)
                cursor.execute("""
                    SELECT * FROM molecules WHERE smiles = ? LIMIT 1
                """, (canonical,))
                row = cursor.fetchone()
        
        if row:
            return self._row_to_record(row)
        return None
    
    def _lookup_by_inchi(self, inchi: str) -> Optional[MoleculeRecord]:
        """Look up by InChI string."""
        cursor = self.conn.cursor()
        cursor.execute("""
            SELECT * FROM molecules WHERE inchi = ? LIMIT 1
        """, (inchi,))
        row = cursor.fetchone()
        
        if row:
            return self._row_to_record(row)
        return None
    
    def _lookup_by_inchikey(self, inchikey: str) -> Optional[MoleculeRecord]:
        """Look up by InChIKey."""
        cursor = self.conn.cursor()
        cursor.execute("""
            SELECT * FROM molecules WHERE inchikey = ? LIMIT 1
        """, (inchikey.upper(),))
        row = cursor.fetchone()
        
        if row:
            return self._row_to_record(row)
        return None
    
    def _lookup_by_cid(self, cid: str) -> Optional[MoleculeRecord]:
        """Look up by PubChem CID."""
        cursor = self.conn.cursor()
        cursor.execute("""
            SELECT * FROM molecules WHERE pubchem_cid = ? LIMIT 1
        """, (int(cid),))
        row = cursor.fetchone()
        
        if row:
            return self._row_to_record(row)
        return None
    
    def _lookup_by_chembl(self, chembl_id: str) -> Optional[MoleculeRecord]:
        """Look up by ChEMBL ID."""
        cursor = self.conn.cursor()
        cursor.execute("""
            SELECT * FROM molecules WHERE chembl_id = ? COLLATE NOCASE LIMIT 1
        """, (chembl_id,))
        row = cursor.fetchone()
        
        if row:
            return self._row_to_record(row)
        return None
    
    def _row_to_record(self, row: sqlite3.Row) -> MoleculeRecord:
        """Convert database row to MoleculeRecord."""
        return MoleculeRecord(
            id=row["id"],
            name=row["name"],
            smiles=row["smiles"],
            inchi=row["inchi"],
            inchikey=row["inchikey"],
            formula=row["formula"],
            molecular_weight=row["molecular_weight"],
            source=row["source"],
            category=row["category"],
            pubchem_cid=row["pubchem_cid"],
            chembl_id=row["chembl_id"],
            created_at=row["created_at"],
        )

    # =========================================================================
    # Enhanced Fuzzy Matching / Suggestions
    # =========================================================================

    def find_suggestions(
        self,
        query: str,
        max_results: int = 10,
        min_similarity: float = 0.3,
        include_iupac_analysis: bool = True,
    ) -> Dict[str, Any]:
        """
        Comprehensive molecule name resolution with fuzzy matching.

        Uses multiple strategies in priority order:
        1. Exact match after normalization (misspelling/abbreviation correction)
        2. Database substring matching
        3. Phonetic matching (sounds-like)
        4. IUPAC fragment analysis (structural hints)
        5. Word-based matching for compound names

        Args:
            query: The molecule name to search for
            max_results: Maximum number of suggestions
            min_similarity: Minimum similarity score (0-1)
            include_iupac_analysis: Whether to analyze IUPAC structure

        Returns:
            Dict with:
            - 'suggestions': List of FuzzySuggestion dicts
            - 'normalization': Normalization info (corrections made)
            - 'iupac_analysis': IUPAC fragment analysis (if enabled)
            - 'alternative_queries': Suggested alternative search terms
        """
        query = query.strip()
        if not query:
            return {
                "suggestions": [],
                "normalization": None,
                "iupac_analysis": None,
                "alternative_queries": []
            }

        result = {
            "original_query": query,
            "suggestions": [],
            "normalization": {},
            "iupac_analysis": None,
            "alternative_queries": [],
        }

        # Step 1: Normalize the query
        normalized = ChemicalNameNormalizer.normalize(query)
        variations = ChemicalNameNormalizer.get_variations(query)

        result["normalization"] = {
            "original": query,
            "normalized": normalized,
            "variations_tried": variations,
            "was_corrected": normalized != query.lower().strip(),
        }

        # Check if normalization found a correction
        if query.lower() in ChemicalNameNormalizer.MISSPELLINGS:
            result["normalization"]["correction_type"] = "misspelling"
            result["normalization"]["corrected_to"] = ChemicalNameNormalizer.MISSPELLINGS[query.lower()]

        if query.lower() in ChemicalNameNormalizer.ABBREVIATIONS:
            result["normalization"]["correction_type"] = "abbreviation"
            result["normalization"]["expanded_to"] = ChemicalNameNormalizer.ABBREVIATIONS[query.lower()]

        # Step 2: Try all variations against database
        seen_names: Set[str] = set()
        suggestions: List[FuzzySuggestion] = []

        for variation in variations:
            # Try exact match for each variation
            record = self._lookup_by_name(variation)
            if record and record.name and record.name.lower() not in seen_names:
                suggestions.append(FuzzySuggestion(
                    name=record.name,
                    smiles=record.smiles,
                    score=1.0,
                    match_type="exact_normalized",
                    source=record.source or "",
                    pubchem_cid=record.pubchem_cid,
                ))
                seen_names.add(record.name.lower())

        # Step 3: IUPAC Fragment Analysis
        if include_iupac_analysis:
            parsed = IUPACFragmentParser.parse(query)
            result["iupac_analysis"] = parsed

            if parsed.get("recognized"):
                # Get related name suggestions based on IUPAC structure
                iupac_suggestions = IUPACFragmentParser.suggest_related_names(parsed)
                result["alternative_queries"].extend(iupac_suggestions)

                # Try these suggestions against the database
                for iupac_name in iupac_suggestions[:5]:  # Limit to top 5
                    record = self._lookup_by_name(iupac_name)
                    if record and record.name and record.name.lower() not in seen_names:
                        suggestions.append(FuzzySuggestion(
                            name=record.name,
                            smiles=record.smiles,
                            score=0.85,
                            match_type="iupac_related",
                            source=record.source or "",
                            pubchem_cid=record.pubchem_cid,
                            details={"iupac_hint": iupac_name}
                        ))
                        seen_names.add(record.name.lower())

        # Step 4: Database fuzzy search (substring, prefix, word-based)
        db_suggestions = self._fuzzy_db_search(normalized, max_results * 2, seen_names)
        for sugg in db_suggestions:
            if sugg.score >= min_similarity:
                suggestions.append(sugg)
                seen_names.add(sugg.name.lower())

        # Step 5: Phonetic matching
        phonetic_suggestions = self._phonetic_search(query, max_results, seen_names)
        for sugg in phonetic_suggestions:
            if sugg.score >= min_similarity:
                suggestions.append(sugg)
                seen_names.add(sugg.name.lower())

        # Sort by score and limit
        suggestions.sort(key=lambda x: x.score, reverse=True)
        result["suggestions"] = [s.to_dict() for s in suggestions[:max_results]]

        # Add alternative query suggestions
        if not suggestions:
            # If no matches, suggest based on IUPAC fragments
            if result.get("iupac_analysis", {}).get("parent"):
                parent_name = result["iupac_analysis"]["parent"]["name"]
                if parent_name not in result["alternative_queries"]:
                    result["alternative_queries"].append(parent_name)

            # Add SMILES suggestion hint
            result["alternative_queries"].append(
                "Try providing the SMILES string directly for exact structure"
            )

        return result

    def _fuzzy_db_search(
        self,
        query: str,
        max_results: int,
        seen_names: Set[str]
    ) -> List[FuzzySuggestion]:
        """
        Perform fuzzy search in the database using multiple strategies.
        """
        suggestions = []
        cursor = self.conn.cursor()
        query_lower = query.lower()

        # Strategy 1: Substring match
        cursor.execute("""
            SELECT DISTINCT name, smiles, source, pubchem_cid
            FROM molecules
            WHERE name LIKE ? COLLATE NOCASE
            AND name IS NOT NULL
            LIMIT ?
        """, (f"%{query_lower}%", max_results))

        for row in cursor.fetchall():
            if row[0] and row[0].lower() not in seen_names:
                score = self._calculate_similarity(query_lower, row[0].lower())
                suggestions.append(FuzzySuggestion(
                    name=row[0],
                    smiles=row[1],
                    score=score,
                    match_type="substring",
                    source=row[2] or "",
                    pubchem_cid=row[3],
                ))

        # Strategy 2: Alias search
        cursor.execute("""
            SELECT DISTINCT m.name, m.smiles, m.source, m.pubchem_cid, a.alias
            FROM molecules m
            JOIN aliases a ON m.id = a.molecule_id
            WHERE a.alias LIKE ? COLLATE NOCASE
            AND m.name IS NOT NULL
            LIMIT ?
        """, (f"%{query_lower}%", max_results))

        for row in cursor.fetchall():
            if row[0] and row[0].lower() not in seen_names:
                score = self._calculate_similarity(query_lower, row[4].lower()) if row[4] else 0.5
                suggestions.append(FuzzySuggestion(
                    name=row[0],
                    smiles=row[1],
                    score=score * 0.95,  # Slight penalty for alias match
                    match_type="alias",
                    source=row[2] or "",
                    pubchem_cid=row[3],
                    details={"matched_alias": row[4]}
                ))

        # Strategy 3: Word-based matching
        words = query_lower.replace('-', ' ').replace(',', ' ').split()
        if len(words) > 1:
            for word in words:
                if len(word) >= 3:
                    cursor.execute("""
                        SELECT DISTINCT name, smiles, source, pubchem_cid
                        FROM molecules
                        WHERE name LIKE ? COLLATE NOCASE
                        AND name IS NOT NULL
                        LIMIT ?
                    """, (f"%{word}%", max_results // 2))

                    for row in cursor.fetchall():
                        if row[0] and row[0].lower() not in seen_names:
                            base_score = self._calculate_similarity(query_lower, row[0].lower())
                            suggestions.append(FuzzySuggestion(
                                name=row[0],
                                smiles=row[1],
                                score=base_score * 0.8,
                                match_type="word_match",
                                source=row[2] or "",
                                pubchem_cid=row[3],
                                details={"matched_word": word}
                            ))

        # Strategy 4: IUPAC prefix matching
        for prefix in IUPACFragmentParser.PREFIXES.keys():
            if prefix in query_lower:
                cursor.execute("""
                    SELECT DISTINCT name, smiles, source, pubchem_cid
                    FROM molecules
                    WHERE name LIKE ? COLLATE NOCASE
                    AND name IS NOT NULL
                    LIMIT ?
                """, (f"%{prefix}%", 5))

                for row in cursor.fetchall():
                    if row[0] and row[0].lower() not in seen_names:
                        base_score = self._calculate_similarity(query_lower, row[0].lower())
                        suggestions.append(FuzzySuggestion(
                            name=row[0],
                            smiles=row[1],
                            score=base_score * 0.7,
                            match_type="iupac_prefix",
                            source=row[2] or "",
                            pubchem_cid=row[3],
                            details={"matched_prefix": prefix}
                        ))

        return suggestions

    def _phonetic_search(
        self,
        query: str,
        max_results: int,
        seen_names: Set[str]
    ) -> List[FuzzySuggestion]:
        """
        Search for phonetically similar molecule names.
        """
        suggestions = []
        cursor = self.conn.cursor()

        # Get a sample of molecule names for phonetic comparison
        # This is expensive, so limit the scope
        cursor.execute("""
            SELECT DISTINCT name, smiles, source, pubchem_cid
            FROM molecules
            WHERE name IS NOT NULL
            AND LENGTH(name) BETWEEN ? AND ?
            LIMIT 1000
        """, (max(1, len(query) - 3), len(query) + 5))

        query_code = PhoneticMatcher.double_metaphone(query)

        for row in cursor.fetchall():
            if row[0] and row[0].lower() not in seen_names:
                name = row[0]
                phonetic_score = PhoneticMatcher.match_score(query, name)

                if phonetic_score >= 0.6:  # High phonetic threshold
                    suggestions.append(FuzzySuggestion(
                        name=name,
                        smiles=row[1],
                        score=phonetic_score * 0.85,  # Scale down phonetic matches
                        match_type="phonetic",
                        source=row[2] or "",
                        pubchem_cid=row[3],
                        details={"phonetic_code": query_code[0]}
                    ))

        return suggestions

    def suggest_similar_names(
        self,
        query: str,
        max_results: int = 10,
        min_similarity: float = 0.3
    ) -> List[Dict[str, Any]]:
        """
        Find similar molecule names using fuzzy matching.

        This is a simplified wrapper around find_suggestions for backward compatibility.

        Args:
            query: The search term
            max_results: Maximum number of suggestions
            min_similarity: Minimum similarity score (0-1)

        Returns:
            List of dicts with 'name', 'smiles', 'score', 'match_type'
        """
        result = self.find_suggestions(query, max_results, min_similarity)
        return result.get("suggestions", [])

    def _calculate_similarity(self, query: str, target: str) -> float:
        """
        Calculate similarity score between two strings.
        Uses a combination of Jaccard similarity, phonetic matching, and prefix bonus.
        """
        if not query or not target:
            return 0.0

        query = query.lower()
        target = target.lower()

        # Exact match
        if query == target:
            return 1.0

        # Check for misspelling correction
        if query in ChemicalNameNormalizer.MISSPELLINGS:
            if ChemicalNameNormalizer.MISSPELLINGS[query] == target:
                return 0.98  # High score for known correction

        # Substring match
        if query in target:
            # Score based on how much of target is the query
            return 0.85 + 0.1 * (len(query) / len(target))

        if target in query:
            return 0.8 + 0.1 * (len(target) / len(query))

        # Calculate Jaccard similarity on character trigrams
        def get_trigrams(s: str) -> Set[str]:
            s = f"  {s}  "  # Pad for edge trigrams
            return {s[i:i+3] for i in range(len(s) - 2)}

        q_trigrams = get_trigrams(query)
        t_trigrams = get_trigrams(target)

        if not q_trigrams or not t_trigrams:
            return 0.0

        intersection = len(q_trigrams & t_trigrams)
        union = len(q_trigrams | t_trigrams)

        jaccard = intersection / union if union > 0 else 0

        # Common prefix bonus
        common_prefix = 0
        for qc, tc in zip(query, target):
            if qc == tc:
                common_prefix += 1
            else:
                break
        prefix_bonus = (common_prefix / max(len(query), len(target))) * 0.2

        # Phonetic similarity bonus
        phonetic_score = PhoneticMatcher.match_score(query, target)
        phonetic_bonus = phonetic_score * 0.15

        return min(1.0, jaccard + prefix_bonus + phonetic_bonus)

    # =========================================================================
    # Add/Cache Methods
    # =========================================================================
    
    def add_molecule(
        self,
        smiles: str,
        name: Optional[str] = None,
        source: str = "user_cache",
        inchi: Optional[str] = None,
        inchikey: Optional[str] = None,
        formula: Optional[str] = None,
        molecular_weight: Optional[float] = None,
        category: Optional[str] = None,
        pubchem_cid: Optional[int] = None,
        chembl_id: Optional[str] = None,
        aliases: Optional[List[str]] = None,
    ) -> int:
        """
        Add a molecule to the database.
        
        Args:
            smiles: SMILES string (required)
            name: Common name
            source: Source identifier
            inchi: InChI string
            inchikey: InChIKey
            formula: Molecular formula
            molecular_weight: Molecular weight
            category: Category (drug, solvent, etc.)
            pubchem_cid: PubChem CID
            chembl_id: ChEMBL ID
            aliases: List of alternative names
        
        Returns:
            Database ID of inserted molecule
        """
        # Generate InChIKey if not provided and RDKit available
        if inchikey is None and RDKIT_AVAILABLE:
            mol = Chem.MolFromSmiles(smiles)
            if mol:
                inchi_str = Chem.MolToInchi(mol) if hasattr(Chem, 'MolToInchi') else None
                if inchi_str:
                    inchikey = Chem.InchiToInchiKey(inchi_str) if hasattr(Chem, 'InchiToInchiKey') else None
        
        cursor = self.conn.cursor()
        
        # Check for duplicate by InChIKey
        if inchikey:
            cursor.execute("SELECT id FROM molecules WHERE inchikey = ?", (inchikey,))
            existing = cursor.fetchone()
            if existing:
                # Update name if not set
                if name:
                    cursor.execute(
                        "UPDATE molecules SET name = COALESCE(name, ?) WHERE id = ?",
                        (name, existing["id"])
                    )
                    self.conn.commit()
                return existing["id"]
        
        # Insert new molecule
        cursor.execute("""
            INSERT INTO molecules (
                name, smiles, inchi, inchikey, formula, molecular_weight,
                source, category, pubchem_cid, chembl_id
            ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
        """, (
            name, smiles, inchi, inchikey, formula, molecular_weight,
            source, category, pubchem_cid, chembl_id
        ))
        
        molecule_id = cursor.lastrowid
        
        # Add aliases
        if aliases:
            for alias in aliases:
                cursor.execute("""
                    INSERT OR IGNORE INTO aliases (alias, molecule_id)
                    VALUES (?, ?)
                """, (alias.lower(), molecule_id))
        
        # Also add name as alias
        if name:
            cursor.execute("""
                INSERT OR IGNORE INTO aliases (alias, molecule_id)
                VALUES (?, ?)
            """, (name.lower(), molecule_id))
        
        self.conn.commit()
        return molecule_id
    
    def add_molecules_bulk(
        self,
        molecules: List[Dict[str, Any]],
        source: str,
        batch_size: int = 10000
    ) -> int:
        """
        Bulk insert molecules efficiently.
        
        Args:
            molecules: List of molecule dicts with at least 'smiles' key
            source: Source identifier
            batch_size: Commit after this many inserts
        
        Returns:
            Number of molecules inserted
        """
        cursor = self.conn.cursor()
        count = 0
        
        for i, mol in enumerate(molecules):
            smiles = mol.get("smiles")
            if not smiles:
                continue
            
            cursor.execute("""
                INSERT OR IGNORE INTO molecules (
                    name, smiles, inchi, inchikey, formula, molecular_weight,
                    source, category, pubchem_cid, chembl_id
                ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
            """, (
                mol.get("name"),
                smiles,
                mol.get("inchi"),
                mol.get("inchikey"),
                mol.get("formula"),
                mol.get("molecular_weight"),
                source,
                mol.get("category"),
                mol.get("pubchem_cid"),
                mol.get("chembl_id"),
            ))
            
            if cursor.rowcount > 0:
                count += 1
                
                # Add name as alias
                if mol.get("name"):
                    cursor.execute("""
                        INSERT OR IGNORE INTO aliases (alias, molecule_id)
                        VALUES (?, ?)
                    """, (mol["name"].lower(), cursor.lastrowid))
            
            # Commit periodically
            if (i + 1) % batch_size == 0:
                self.conn.commit()
                logger.info(f"Imported {i + 1} molecules...")
        
        self.conn.commit()
        return count
    
    # =========================================================================
    # Database Import Methods
    # =========================================================================
    
    def _attach_chembl_database(self, db_path: Path) -> int:
        """
        Attach ChEMBL SQLite database and import molecules.
        
        ChEMBL database can be downloaded from:
        https://ftp.ebi.ac.uk/pub/databases/chembl/ChEMBLdb/releases/
        
        Args:
            db_path: Path to ChEMBL SQLite database
        
        Returns:
            Number of molecules imported
        """
        if self._is_source_imported(db_path.name):
            logger.info(f"ChEMBL database already imported: {db_path.name}")
            return 0
        
        logger.info(f"Importing ChEMBL database: {db_path}")
        
        cursor = self.conn.cursor()
        
        # Attach ChEMBL database
        cursor.execute(f"ATTACH DATABASE ? AS chembl", (str(db_path),))
        
        # Import molecules from ChEMBL
        logger.info("Importing molecules from ChEMBL (this may take several minutes)...")
        cursor.execute("""
            INSERT OR IGNORE INTO molecules (
                name, smiles, chembl_id, source
            )
            SELECT 
                md.pref_name,
                cs.canonical_smiles,
                md.chembl_id,
                'chembl'
            FROM chembl.molecule_dictionary md
            JOIN chembl.compound_structures cs ON md.molregno = cs.molregno
            WHERE cs.canonical_smiles IS NOT NULL
        """)
        
        count = cursor.rowcount
        
        # Commit BEFORE detaching to release locks
        self.conn.commit()
        
        # Detach
        cursor.execute("DETACH DATABASE chembl")
        
        # Register source
        self._register_source(db_path.name, "chembl", count, str(db_path))
        
        logger.info(f"Imported {count} molecules from ChEMBL")
        
        return count
    
    def _import_pubchem_smiles(self, file_path: Path) -> int:
        """
        Import PubChem CID-SMILES file.
        
        Download from: ftp://ftp.ncbi.nlm.nih.gov/pubchem/Compound/Extras/CID-SMILES.gz
        
        Args:
            file_path: Path to CID-SMILES.gz file
        
        Returns:
            Number of molecules imported
        """
        if self._is_source_imported(file_path.name):
            logger.info(f"PubChem SMILES already imported: {file_path.name}")
            return 0
        
        logger.info(f"Importing PubChem SMILES: {file_path}")
        
        count = 0
        batch = []
        batch_size = 100000
        
        # Open gzipped file
        open_func = gzip.open if str(file_path).endswith('.gz') else open
        
        with open_func(file_path, 'rt') as f:
            for line in f:
                parts = line.strip().split('\t')
                if len(parts) >= 2:
                    cid, smiles = parts[0], parts[1]
                    if cid.isdigit():
                        batch.append({
                            "smiles": smiles,
                            "pubchem_cid": int(cid),
                        })
                
                if len(batch) >= batch_size:
                    count += self.add_molecules_bulk(batch, source="pubchem")
                    batch = []
                    logger.info(f"Imported {count} PubChem molecules...")
        
        if batch:
            count += self.add_molecules_bulk(batch, source="pubchem")
        
        self._register_source(file_path.name, "pubchem", count, str(file_path))
        logger.info(f"Imported {count} molecules from PubChem")
        
        return count
    
    def _import_pubchemlite(self, file_path: Path) -> int:
        """
        Import PubChemLite CSV file.
        
        Download from: https://zenodo.org/records/14251246
        
        Args:
            file_path: Path to PubChemLite CSV
        
        Returns:
            Number of molecules imported
        """
        if self._is_source_imported(file_path.name):
            logger.info(f"PubChemLite already imported: {file_path.name}")
            return 0
        
        logger.info(f"Importing PubChemLite: {file_path}")
        
        molecules = []
        
        with open(file_path, 'r', encoding='utf-8') as f:
            reader = csv.DictReader(f)
            for row in reader:
                mol = {
                    "smiles": row.get("SMILES") or row.get("smiles"),
                    "name": row.get("Name") or row.get("name"),
                    "inchikey": row.get("InChIKey") or row.get("inchikey"),
                    "formula": row.get("MolecularFormula") or row.get("formula"),
                }
                if mol["smiles"]:
                    molecules.append(mol)
        
        count = self.add_molecules_bulk(molecules, source="pubchemlite")
        self._register_source(file_path.name, "pubchemlite", count, str(file_path))
        
        logger.info(f"Imported {count} molecules from PubChemLite")
        return count
    
    def import_database(self, file_path: Union[str, Path]) -> int:
        """
        Auto-detect and import a database file.
        
        Supports:
        - ChEMBL SQLite (.db)
        - PubChem CID-SMILES (.gz)
        - CSV with SMILES column
        - SDF files
        
        Args:
            file_path: Path to database file
        
        Returns:
            Number of molecules imported
        """
        file_path = Path(file_path)
        
        if not file_path.exists():
            raise FileNotFoundError(f"Database file not found: {file_path}")
        
        name = file_path.name.lower()
        
        if "chembl" in name and name.endswith(".db"):
            return self._attach_chembl_database(file_path)
        
        if name.endswith(".gz") and "smiles" in name.lower():
            return self._import_pubchem_smiles(file_path)
        
        if "pubchemlite" in name and name.endswith(".csv"):
            return self._import_pubchemlite(file_path)
        
        if name.endswith(".csv"):
            return self._import_generic_csv(file_path)
        
        raise ValueError(f"Unknown database format: {file_path}")
    
    def _import_generic_csv(self, file_path: Path) -> int:
        """Import generic CSV with SMILES column."""
        molecules = []
        
        with open(file_path, 'r', encoding='utf-8') as f:
            reader = csv.DictReader(f)
            headers_lower = {h.lower(): h for h in reader.fieldnames or []}
            
            # Find SMILES column
            smiles_col = headers_lower.get("smiles") or headers_lower.get("canonical_smiles")
            name_col = headers_lower.get("name") or headers_lower.get("compound_name")
            
            if not smiles_col:
                raise ValueError("CSV must have a 'smiles' or 'canonical_smiles' column")
            
            for row in reader:
                mol = {"smiles": row.get(smiles_col)}
                if name_col:
                    mol["name"] = row.get(name_col)
                if mol["smiles"]:
                    molecules.append(mol)
        
        source_name = file_path.stem
        count = self.add_molecules_bulk(molecules, source=source_name)
        self._register_source(file_path.name, source_name, count, str(file_path))
        
        return count
    
    # =========================================================================
    # Statistics and Info
    # =========================================================================
    
    def stats(self) -> Dict[str, Any]:
        """Get database statistics."""
        cursor = self.conn.cursor()
        
        # Total molecules
        cursor.execute("SELECT COUNT(*) as count FROM molecules")
        total = cursor.fetchone()["count"]
        
        # By source
        cursor.execute("""
            SELECT source, COUNT(*) as count 
            FROM molecules 
            GROUP BY source
            ORDER BY count DESC
        """)
        by_source = {row["source"]: row["count"] for row in cursor.fetchall()}
        
        # Imported sources
        cursor.execute("SELECT name, version, molecule_count, imported_at FROM sources")
        sources = [dict(row) for row in cursor.fetchall()]
        
        return {
            "total_molecules": total,
            "by_source": by_source,
            "imported_sources": sources,
            "database_path": str(self.db_path),
            "database_size_mb": self.db_path.stat().st_size / (1024 * 1024) if self.db_path.exists() else 0,
        }
    
    def search(
        self, 
        query: str, 
        limit: int = 20,
        search_type: str = "prefix"
    ) -> List[MoleculeRecord]:
        """
        Search for molecules by name.
        
        Args:
            query: Search query
            limit: Maximum results
            search_type: "prefix" (starts with) or "contains"
        
        Returns:
            List of matching molecules
        """
        cursor = self.conn.cursor()
        
        if search_type == "prefix":
            pattern = f"{query}%"
        else:
            pattern = f"%{query}%"
        
        cursor.execute("""
            SELECT DISTINCT m.* FROM molecules m
            LEFT JOIN aliases a ON m.id = a.molecule_id
            WHERE m.name LIKE ? COLLATE NOCASE
               OR a.alias LIKE ? COLLATE NOCASE
            LIMIT ?
        """, (pattern, pattern, limit))
        
        return [self._row_to_record(row) for row in cursor.fetchall()]
    
    def close(self) -> None:
        """Close database connection."""
        if self._conn:
            self._conn.close()
            self._conn = None


# Global instance for convenience
_default_db: Optional[MoleculeDatabase] = None


def get_molecule_database() -> MoleculeDatabase:
    """Get the global molecule database instance."""
    global _default_db
    if _default_db is None:
        _default_db = MoleculeDatabase()
    return _default_db


def lookup_molecule(identifier: str, id_type: str = "auto") -> Optional[Dict]:
    """
    Convenience function to look up a molecule.
    
    Args:
        identifier: Any molecule identifier
        id_type: Type hint ("auto", "name", "smiles", etc.)
    
    Returns:
        Molecule dict or None
    """
    db = get_molecule_database()
    record = db.lookup(identifier, id_type)
    if record:
        return record.to_dict()
    return None
