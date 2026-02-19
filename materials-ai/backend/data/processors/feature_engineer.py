"""
Feature engineering — compute molecular descriptors from SMILES.
Uses RDKit for descriptor calculation.
"""

import logging
from typing import Optional

logger = logging.getLogger(__name__)


class FeatureEngineer:
    """
    Compute molecular descriptors (features) from SMILES strings.
    Generates a rich feature set for ML model training.
    """

    def compute_features(self, smiles: str) -> Optional[dict]:
        """
        Compute molecular descriptors for a SMILES string.

        Returns dict with:
          - Basic descriptors (MW, LogP, TPSA, HBD, HBA, ...)
          - Morgan fingerprint (2048-bit)
          - Additional RDKit descriptors
        """
        try:
            from rdkit import Chem
            from rdkit.Chem import Descriptors, AllChem, rdMolDescriptors

            mol = Chem.MolFromSmiles(smiles)
            if mol is None:
                logger.warning(f"Invalid SMILES: {smiles}")
                return None

            features = {
                # Basic descriptors
                "molecular_weight": round(Descriptors.MolWt(mol), 2),
                "logp": round(Descriptors.MolLogP(mol), 3),
                "hbd": Descriptors.NumHDonors(mol),
                "hba": Descriptors.NumHAcceptors(mol),
                "tpsa": round(Descriptors.TPSA(mol), 2),
                "rotatable_bonds": Descriptors.NumRotatableBonds(mol),
                "aromatic_rings": Descriptors.NumAromaticRings(mol),

                # Additional descriptors
                "num_heavy_atoms": mol.GetNumHeavyAtoms(),
                "num_rings": Descriptors.RingCount(mol),
                "fraction_csp3": round(Descriptors.FractionCSP3(mol), 3),
                "num_heteroatoms": rdMolDescriptors.CalcNumHeteroatoms(mol),
                "num_amide_bonds": rdMolDescriptors.CalcNumAmideBonds(mol),

                # Complexity
                "bertz_ct": round(Descriptors.BertzCT(mol), 2),

                # Morgan fingerprint as list (for storage)
                "morgan_fp": list(
                    AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=2048)
                ),
            }

            return features
        except ImportError:
            logger.error("RDKit not installed. pip install rdkit-pypi")
            return None
        except Exception as e:
            logger.error(f"Feature computation failed for {smiles}: {e}")
            return None

    def compute_batch(self, smiles_list: list[str]) -> list[Optional[dict]]:
        """Compute features for a batch of SMILES strings."""
        return [self.compute_features(s) for s in smiles_list]
