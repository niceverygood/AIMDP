"""
SDF/MOL file loader — loads molecular structures from SDF and MOL files.
Uses RDKit for parsing.
"""

import logging
from typing import Optional

logger = logging.getLogger(__name__)


class SDFLoader:
    """Load molecular data from SDF (Structure Data File) and MOL files."""

    async def load(self, file_path: str) -> list[dict]:
        """
        Load molecules from SDF or MOL file.

        SDF files can contain multiple molecules with associated properties.
        MOL files contain a single molecule.

        Returns list of dicts with 'smiles' and any available properties.
        """
        try:
            from rdkit import Chem
            from rdkit.Chem import Descriptors

            records = []

            if file_path.endswith(".mol"):
                # Single molecule MOL file
                mol = Chem.MolFromMolFile(file_path)
                if mol:
                    records.append(self._mol_to_record(mol))
            else:
                # SDF file (can contain multiple molecules)
                supplier = Chem.SDMolSupplier(file_path, removeHs=True)

                for idx, mol in enumerate(supplier):
                    if mol is None:
                        logger.warning(f"Could not parse molecule at index {idx}")
                        continue

                    record = self._mol_to_record(mol)

                    # Extract SDF properties
                    for prop_name in mol.GetPropsAsDict():
                        value = mol.GetPropsAsDict()[prop_name]
                        record[prop_name.lower().replace(" ", "_")] = value

                    records.append(record)

            logger.info(f"Loaded {len(records)} molecules from {file_path}")
            return records

        except ImportError:
            logger.error("RDKit not installed. pip install rdkit-pypi")
            raise
        except Exception as e:
            logger.error(f"SDF load failed: {e}")
            raise

    def _mol_to_record(self, mol) -> dict:
        """Convert an RDKit Mol object to a data record."""
        from rdkit import Chem
        from rdkit.Chem import Descriptors

        smiles = Chem.MolToSmiles(mol)

        record = {
            "smiles": smiles,
            "name": mol.GetProp("_Name") if mol.HasProp("_Name") else None,
            "molecular_weight": round(Descriptors.MolWt(mol), 2),
            "logp": round(Descriptors.MolLogP(mol), 3),
            "hbd": Descriptors.NumHDonors(mol),
            "hba": Descriptors.NumHAcceptors(mol),
            "tpsa": round(Descriptors.TPSA(mol), 2),
        }

        return record
