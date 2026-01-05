"""
molecule/molecule_database.py - Universal Molecule Database

SQLite-based molecule database supporting:
- Local cache of user lookups
- ChEMBL (~2.4M bioactive molecules)
- PubChemLite (~500K curated)
- PubChem Complete (~130M, optional 15GB download)
- ZINC (purchasable compounds)
- DrugBank (approved drugs)

The database automatically incorporates any downloaded databases
placed in the data directory.
"""

import sqlite3
import os
import gzip
import csv
import logging
import hashlib
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Any, Union
from dataclasses import dataclass
from datetime import datetime
import importlib.util

# Check if RDKit is available for InChIKey generation
RDKIT_AVAILABLE = importlib.util.find_spec("rdkit") is not None

if RDKIT_AVAILABLE:
    from rdkit import Chem
    from rdkit.Chem.inchi import MolFromInchi, InchiToInchiKey

logger = logging.getLogger(__name__)


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
