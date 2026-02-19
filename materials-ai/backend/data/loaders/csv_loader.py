"""
CSV data loader — loads material data from CSV/Excel files.
"""

import logging
from typing import Optional

logger = logging.getLogger(__name__)


class CSVLoader:
    """Load material data from CSV or Excel files."""

    # Expected column name mappings
    COLUMN_ALIASES = {
        "SMILES": "smiles",
        "smiles": "smiles",
        "Smiles": "smiles",
        "Name": "name",
        "name": "name",
        "Category": "category",
        "category": "category",
        "Thermal Stability": "thermal_stability",
        "Td": "thermal_stability",
        "Dielectric": "dielectric_constant",
        "Dk": "dielectric_constant",
        "Bandgap": "bandgap",
        "Band Gap": "bandgap",
        "Eg": "bandgap",
        "Solubility": "solubility",
        "Density": "density",
    }

    async def load(self, file_path: str) -> list[dict]:
        """
        Load data from CSV or Excel file.

        Supports:
          - .csv files (comma or tab delimited)
          - .xlsx, .xls files
          - Auto-detects column names and maps to canonical names
        """
        import pandas as pd

        try:
            if file_path.endswith((".xlsx", ".xls")):
                df = pd.read_excel(file_path)
            else:
                # Try comma first, then tab
                try:
                    df = pd.read_csv(file_path)
                except Exception:
                    df = pd.read_csv(file_path, sep="\t")

            # Rename columns to canonical names
            rename_map = {}
            for col in df.columns:
                canonical = self.COLUMN_ALIASES.get(col.strip())
                if canonical:
                    rename_map[col] = canonical
            df = df.rename(columns=rename_map)

            # Ensure SMILES column exists
            if "smiles" not in df.columns:
                # Try to find a column that looks like SMILES
                for col in df.columns:
                    sample = df[col].dropna().head(5)
                    if sample.astype(str).str.contains(r"[A-Za-z]\(", regex=True).any():
                        df = df.rename(columns={col: "smiles"})
                        break

            # Convert to list of dicts, drop NaN
            records = df.where(df.notna(), None).to_dict(orient="records")

            logger.info(f"Loaded {len(records)} records from {file_path}")
            return records

        except Exception as e:
            logger.error(f"CSV load failed: {e}")
            raise
