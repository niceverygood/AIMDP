"""
Inference pipeline — load trained models and run predictions.
"""

import logging
from pathlib import Path
from typing import Optional

import numpy as np
import torch

logger = logging.getLogger(__name__)

PROPERTY_NAMES = [
    "thermal_stability",
    "dielectric_constant",
    "bandgap",
    "solubility",
    "density",
]


class InferencePipeline:
    """
    End-to-end inference pipeline:
    SMILES → Embedding → Property Prediction → Results
    """

    def __init__(
        self,
        model_path: str = "./ml/property_models/predictor.pt",
        device: Optional[str] = None,
    ):
        self._model_path = model_path
        self._device = device or ("cuda" if torch.cuda.is_available() else "cpu")
        self._model = None
        self._embedding_service = None

    def _ensure_loaded(self):
        if self._model is not None:
            return

        model_file = Path(self._model_path)
        if model_file.exists():
            self._model = torch.load(
                str(model_file), map_location=self._device
            )
            self._model.eval()
            logger.info(f"Model loaded from {self._model_path}")
        else:
            logger.warning(
                f"Model file not found at {self._model_path}. "
                "Using heuristic predictor."
            )

        from app.services.embedding import EmbeddingService
        self._embedding_service = EmbeddingService()

    async def predict(self, smiles: str) -> dict:
        """
        Predict properties for a single molecule.

        Returns:
            dict with property names as keys and predicted values.
        """
        self._ensure_loaded()

        # Get embedding
        embedding = await self._embedding_service.get_embedding(smiles)

        if self._model is not None:
            with torch.no_grad():
                input_tensor = torch.tensor(
                    [embedding], dtype=torch.float32, device=self._device
                )
                output = self._model(input_tensor)
                values = output[0].cpu().numpy()

            return {
                name: float(round(values[i], 4))
                for i, name in enumerate(PROPERTY_NAMES)
            }
        else:
            # Fallback to heuristic
            from app.services.predictor import PropertyPredictor
            predictor = PropertyPredictor()
            result = predictor._predict_heuristic(smiles)
            return result.model_dump(exclude={"smiles", "confidence"})

    async def predict_batch(self, smiles_list: list[str]) -> list[dict]:
        """Predict properties for multiple molecules."""
        return [await self.predict(s) for s in smiles_list]

    async def rank_candidates(
        self,
        smiles_list: list[str],
        target_properties: dict[str, float],
        top_k: int = 10,
    ) -> list[dict]:
        """
        Rank candidate molecules by how well they match target properties.

        Args:
            smiles_list: List of candidate SMILES
            target_properties: Dict of target property values
            top_k: Number of top candidates to return

        Returns:
            Sorted list of dicts with SMILES, predicted properties, and score
        """
        predictions = await self.predict_batch(smiles_list)

        scored = []
        for smiles, pred in zip(smiles_list, predictions):
            score = self._compute_match_score(pred, target_properties)
            scored.append({
                "smiles": smiles,
                "predictions": pred,
                "match_score": score,
            })

        scored.sort(key=lambda x: x["match_score"], reverse=True)
        return scored[:top_k]

    def _compute_match_score(
        self,
        predicted: dict,
        target: dict[str, float],
    ) -> float:
        """Compute how well predicted properties match target values."""
        ranges = {
            "thermal_stability": (0, 500),
            "dielectric_constant": (0, 25),
            "bandgap": (0, 8),
            "solubility": (0, 1),
            "density": (0.5, 3.0),
        }

        total = 0.0
        count = 0

        for prop, target_val in target.items():
            if prop in predicted and predicted[prop] is not None:
                low, high = ranges.get(prop, (0, 1))
                if high - low > 0:
                    pred_norm = (predicted[prop] - low) / (high - low)
                    target_norm = (target_val - low) / (high - low)
                    proximity = 1.0 - abs(pred_norm - target_norm)
                    total += max(0, proximity)
                    count += 1

        return round(total / max(count, 1), 4)
