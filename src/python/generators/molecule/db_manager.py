"""
molecule/db_manager.py - Database Manager CLI

Download, import, and manage molecule databases.

Usage:
    python -m generators.molecule.db_manager download chembl
    python -m generators.molecule.db_manager download pubchemlite
    python -m generators.molecule.db_manager download pubchem-complete
    python -m generators.molecule.db_manager import /path/to/database.csv
    python -m generators.molecule.db_manager status
    python -m generators.molecule.db_manager search "aspirin"
"""

import argparse
import sys
import os
import urllib.request
import tarfile
import gzip
import shutil
import logging
from pathlib import Path
from typing import Optional
import importlib.util

from .molecule_database import MoleculeDatabase, get_molecule_database

logging.basicConfig(level=logging.INFO, format='%(levelname)s: %(message)s')
logger = logging.getLogger(__name__)


# Database download sources
DOWNLOAD_SOURCES = {
    "chembl": {
        "name": "ChEMBL 34",
        "description": "Bioactive molecules with drug-like properties (~2.4M molecules)",
        "url": "https://ftp.ebi.ac.uk/pub/databases/chembl/ChEMBLdb/releases/chembl_34/chembl_34_sqlite.tar.gz",
        "size": "~1.5GB compressed, ~4GB extracted",
        "file_type": "sqlite",
    },
    "pubchemlite": {
        "name": "PubChemLite",
        "description": "Curated subset for environmental/metabolomics (~500K molecules)",
        "url": "https://zenodo.org/records/14251246/files/PubChemLite_CCSbase.csv",
        "size": "~50MB",
        "file_type": "csv",
    },
    "pubchem-complete": {
        "name": "PubChem Complete CID-SMILES",
        "description": "ALL PubChem compounds (~130M molecules)",
        "url": "ftp://ftp.ncbi.nlm.nih.gov/pubchem/Compound/Extras/CID-SMILES.gz",
        "size": "~15GB compressed",
        "file_type": "smiles_gz",
        "note": "Very large! Only download if you have sufficient disk space.",
    },
    "zinc-druglike": {
        "name": "ZINC Drug-like",
        "description": "Purchasable drug-like compounds",
        "url": "https://zinc.docking.org/tranches/download",
        "size": "Variable",
        "file_type": "smiles",
        "note": "Requires manual download from ZINC website.",
    },
}


class DatabaseManager:
    """Manage molecule database downloads and imports."""
    
    def __init__(self, data_dir: Optional[Path] = None):
        self.db = get_molecule_database()
        self.data_dir = data_dir or self.db.data_dir
        self.downloads_dir = self.data_dir / "downloads"
        self.downloads_dir.mkdir(parents=True, exist_ok=True)
    
    def download(self, source_key: str) -> bool:
        """
        Download a database.
        
        Args:
            source_key: Key from DOWNLOAD_SOURCES (chembl, pubchemlite, etc.)
        
        Returns:
            True if successful
        """
        if source_key not in DOWNLOAD_SOURCES:
            logger.error(f"Unknown source: {source_key}")
            logger.info(f"Available sources: {', '.join(DOWNLOAD_SOURCES.keys())}")
            return False
        
        source = DOWNLOAD_SOURCES[source_key]
        url = source["url"]
        
        logger.info(f"Downloading {source['name']}...")
        logger.info(f"  URL: {url}")
        logger.info(f"  Size: {source['size']}")
        
        if source.get("note"):
            logger.warning(f"  Note: {source['note']}")
        
        # Determine output filename
        filename = url.split("/")[-1]
        output_path = self.downloads_dir / filename
        
        # Download
        success = self._download_file(url, output_path)
        
        if not success:
            return False
        
        # Extract if needed
        extracted_path = self._extract_if_needed(output_path, source["file_type"])
        
        # Import into database
        if extracted_path:
            logger.info(f"Importing into database...")
            count = self.db.import_database(extracted_path)
            logger.info(f"Successfully imported {count} molecules")
        
        return True
    
    def _download_file(self, url: str, output_path: Path) -> bool:
        """Download a file with progress indicator."""
        if output_path.exists():
            logger.info(f"File already exists: {output_path}")
            return True
        
        logger.info(f"Downloading to {output_path}...")
        
        # Handle FTP URLs
        if url.startswith("ftp://"):
            # Convert FTP to HTTP for ncbi
            if "ftp.ncbi.nlm.nih.gov" in url:
                url = url.replace("ftp://ftp.ncbi.nlm.nih.gov", "https://ftp.ncbi.nlm.nih.gov")
        
        request = urllib.request.Request(url, headers={"User-Agent": "Crystal-MCP-Server/1.0"})
        
        response = urllib.request.urlopen(request)
        total_size = int(response.headers.get('Content-Length', 0))
        
        # Download with progress
        block_size = 8192
        downloaded = 0
        
        with open(output_path, 'wb') as f:
            while True:
                buffer = response.read(block_size)
                if not buffer:
                    break
                f.write(buffer)
                downloaded += len(buffer)
                
                if total_size > 0:
                    percent = (downloaded / total_size) * 100
                    mb_downloaded = downloaded / (1024 * 1024)
                    mb_total = total_size / (1024 * 1024)
                    print(f"\r  Progress: {mb_downloaded:.1f} / {mb_total:.1f} MB ({percent:.1f}%)", end="")
        
        print()  # Newline after progress
        logger.info(f"Download complete: {output_path}")
        return True
    
    def _extract_if_needed(self, file_path: Path, file_type: str) -> Optional[Path]:
        """Extract compressed files if needed."""
        if file_type == "sqlite" and str(file_path).endswith(".tar.gz"):
            # Extract tar.gz
            logger.info(f"Extracting {file_path}...")
            
            with tarfile.open(file_path, "r:gz") as tar:
                # Find the .db file
                for member in tar.getmembers():
                    if member.name.endswith(".db"):
                        tar.extract(member, self.downloads_dir)
                        extracted_path = self.downloads_dir / member.name
                        logger.info(f"Extracted: {extracted_path}")
                        return extracted_path
        
        return file_path
    
    def import_file(self, file_path: str) -> bool:
        """
        Import a database file.
        
        Args:
            file_path: Path to database file
        
        Returns:
            True if successful
        """
        path = Path(file_path)
        
        if not path.exists():
            logger.error(f"File not found: {path}")
            return False
        
        logger.info(f"Importing {path}...")
        count = self.db.import_database(path)
        logger.info(f"Imported {count} molecules")
        
        return True
    
    def status(self) -> None:
        """Print database status."""
        stats = self.db.stats()
        
        print("\n" + "=" * 60)
        print("MOLECULE DATABASE STATUS")
        print("=" * 60)
        print(f"Database path: {stats['database_path']}")
        print(f"Database size: {stats['database_size_mb']:.2f} MB")
        print(f"Total molecules: {stats['total_molecules']:,}")
        
        print("\nMolecules by source:")
        for source, count in stats["by_source"].items():
            print(f"  {source}: {count:,}")
        
        print("\nImported databases:")
        if stats["imported_sources"]:
            for src in stats["imported_sources"]:
                print(f"  - {src['name']} ({src['molecule_count']:,} molecules)")
        else:
            print("  (none)")
        
        print("\nAvailable downloads:")
        for key, source in DOWNLOAD_SOURCES.items():
            status = "✓ imported" if self._is_imported(key) else "  available"
            print(f"  {status} {key}: {source['name']} ({source['size']})")
        
        print()
    
    def _is_imported(self, source_key: str) -> bool:
        """Check if a source has been imported."""
        # Check downloads directory
        source = DOWNLOAD_SOURCES.get(source_key, {})
        url = source.get("url", "")
        filename = url.split("/")[-1] if url else ""
        
        # Check in sources table
        cursor = self.db.conn.cursor()
        cursor.execute(
            "SELECT 1 FROM sources WHERE name LIKE ?",
            (f"%{source_key}%",)
        )
        return cursor.fetchone() is not None
    
    def search(self, query: str, limit: int = 20) -> None:
        """Search for molecules."""
        results = self.db.search(query, limit=limit)
        
        print(f"\nSearch results for '{query}':")
        print("-" * 60)
        
        if not results:
            print("No results found.")
            return
        
        for mol in results:
            print(f"  {mol.name or '(unnamed)'}")
            print(f"    SMILES: {mol.smiles[:50]}..." if len(mol.smiles) > 50 else f"    SMILES: {mol.smiles}")
            print(f"    Source: {mol.source}")
            if mol.formula:
                print(f"    Formula: {mol.formula}")
            print()


def list_sources() -> None:
    """List available database sources."""
    print("\nAvailable molecule databases for download:")
    print("=" * 70)
    
    for key, source in DOWNLOAD_SOURCES.items():
        print(f"\n{key}:")
        print(f"  Name: {source['name']}")
        print(f"  Description: {source['description']}")
        print(f"  Size: {source['size']}")
        print(f"  URL: {source['url']}")
        if source.get("note"):
            print(f"  Note: {source['note']}")


def main():
    parser = argparse.ArgumentParser(
        description="Molecule Database Manager",
        epilog="Example: python -m generators.molecule.db_manager download chembl"
    )
    
    subparsers = parser.add_subparsers(dest="command", help="Commands")
    
    # Download command
    download_parser = subparsers.add_parser("download", help="Download a database")
    download_parser.add_argument(
        "source",
        choices=list(DOWNLOAD_SOURCES.keys()),
        help="Database source to download"
    )
    
    # Import command
    import_parser = subparsers.add_parser("import", help="Import a database file")
    import_parser.add_argument("file", help="Path to database file")
    
    # Status command
    subparsers.add_parser("status", help="Show database status")
    
    # List command
    subparsers.add_parser("list", help="List available database sources")
    
    # Search command
    search_parser = subparsers.add_parser("search", help="Search for molecules")
    search_parser.add_argument("query", help="Search query")
    search_parser.add_argument("-n", "--limit", type=int, default=20, help="Max results")
    
    args = parser.parse_args()
    
    if not args.command:
        parser.print_help()
        return
    
    manager = DatabaseManager()
    
    if args.command == "download":
        manager.download(args.source)
    elif args.command == "import":
        manager.import_file(args.file)
    elif args.command == "status":
        manager.status()
    elif args.command == "list":
        list_sources()
    elif args.command == "search":
        manager.search(args.query, args.limit)


if __name__ == "__main__":
    main()
