"""
Property prediction service.
Predicts material properties (thermal stability, dielectric, bandgap, etc.)
from molecular structure (SMILES).

Uses a fine-tuned neural network head on top of molecular embeddings.
Falls back to RDKit descriptor-based heuristics if the model is not trained yet.
"""

import logging
import numpy as np
from typing import Optional

from app.models.schemas import PropertyPrediction

logger = logging.getLogger(__name__)


class PropertyPredictor:
    """
    Molecular property prediction service.

    Architecture:
        SMILES → Uni-Mol embedding (768-dim) → FC head → [thermal, dielectric, bandgap, solubility, density]

    If the trained model is not available, uses RDKit descriptor-based heuristic estimates.
    """

    def __init__(self, model_path: Optional[str] = None):
        self._model = None
        self._model_path = model_path
        self._initialized = False

    def _ensure_initialized(self):
        if self._initialized:
            return

        if self._model_path:
            try:
                import torch

                self._model = torch.load(self._model_path, map_location="cpu")
                self._model.eval()
                logger.info(f"Property prediction model loaded from {self._model_path}")
            except Exception as e:
                logger.warning(f"Could not load model: {e}. Using heuristic predictor.")

        self._initialized = True

    async def predict(self, smiles: str) -> PropertyPrediction:
        """
        Predict properties for a single molecule.
        """
        self._ensure_initialized()

        if self._model is not None:
            return await self._predict_with_model(smiles)
        else:
            return self._predict_heuristic(smiles)

    async def predict_batch(self, smiles_list: list[str]) -> list[PropertyPrediction]:
        """Predict properties for a batch of molecules."""
        self._ensure_initialized()
        return [await self.predict(s) for s in smiles_list]

    async def _predict_with_model(self, smiles: str) -> PropertyPrediction:
        """Use trained neural network for prediction."""
        import torch
        from app.services.embedding import EmbeddingService

        embedding_service = EmbeddingService()
        embedding = await embedding_service.get_embedding(smiles)

        with torch.no_grad():
            input_tensor = torch.tensor([embedding], dtype=torch.float32)
            output = self._model(input_tensor)
            predictions = output[0].numpy()

        return PropertyPrediction(
            smiles=smiles,
            thermal_stability=float(predictions[0]),
            dielectric_constant=float(predictions[1]),
            bandgap=float(predictions[2]),
            solubility=float(predictions[3]),
            density=float(predictions[4]),
            confidence=0.85,  # TODO: implement uncertainty estimation
        )

    def _predict_heuristic(self, smiles: str) -> PropertyPrediction:
        """
        Heuristic property estimation using RDKit descriptors.
        This is a rough approximation — to be replaced with the trained model.
        """
        try:
            from rdkit import Chem
            from rdkit.Chem import Descriptors

            mol = Chem.MolFromSmiles(smiles)
            if mol is None:
                return PropertyPrediction(smiles=smiles)

            mw = Descriptors.MolWt(mol)
            logp = Descriptors.MolLogP(mol)
            tpsa = Descriptors.TPSA(mol)
            n_aromatic = Descriptors.NumAromaticRings(mol)
            n_heavy = mol.GetNumHeavyAtoms()

            # Heuristic estimations (rough correlations)
            thermal_stability = min(500, max(100, 150 + n_aromatic * 50 + mw * 0.2))
            dielectric_constant = max(1, min(25, 2.0 + tpsa * 0.05 + logp * -0.3))
            bandgap = max(0, min(8, 4.0 - n_aromatic * 0.5 - logp * 0.1))
            solubility = max(0, min(1, 0.5 - logp * 0.1 + tpsa * 0.003))
            density = max(0.5, min(3.0, 0.8 + mw * 0.002 + n_heavy * 0.01))

            return PropertyPrediction(
                smiles=smiles,
                thermal_stability=round(thermal_stability, 1),
                dielectric_constant=round(dielectric_constant, 2),
                bandgap=round(bandgap, 2),
                solubility=round(solubility, 3),
                density=round(density, 3),
                confidence=0.3,  # Low confidence for heuristic
            )
        except Exception as e:
            logger.error(f"Heuristic prediction failed for {smiles}: {e}")
            return PropertyPrediction(smiles=smiles)
