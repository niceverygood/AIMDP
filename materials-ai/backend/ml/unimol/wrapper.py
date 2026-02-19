"""
Uni-Mol model wrapper — unified interface for molecular embeddings and property prediction.
Uni-Mol is a 3D-aware molecular representation model from DP Technology.
"""

import logging
from typing import Optional

import numpy as np

logger = logging.getLogger(__name__)


class UniMolWrapper:
    """
    Wrapper for the Uni-Mol model.

    Provides:
      1. Molecular embeddings (768-dim vectors)
      2. Property prediction via fine-tuned heads

    Usage:
        wrapper = UniMolWrapper()
        embedding = wrapper.get_embedding("CCO")  # ethanol
        properties = wrapper.predict_properties("CCO")
    """

    def __init__(self, model_name: str = "unimol_base"):
        self._model_name = model_name
        self._repr_model = None
        self._initialized = False

    def _ensure_initialized(self):
        """Lazy-load the Uni-Mol model."""
        if self._initialized:
            return

        try:
            from unimol_tools import UniMolRepr
            self._repr_model = UniMolRepr(model_name=self._model_name)
            logger.info(f"Uni-Mol ({self._model_name}) loaded successfully")
        except ImportError:
            logger.warning(
                "unimol-tools not installed. "
                "Install with: pip install unimol-tools"
            )
        except Exception as e:
            logger.error(f"Uni-Mol initialization failed: {e}")

        self._initialized = True

    def get_embedding(self, smiles: str) -> Optional[np.ndarray]:
        """Get 768-dim embedding for a single molecule."""
        self._ensure_initialized()

        if self._repr_model is None:
            return None

        try:
            reprs = self._repr_model.get_repr([smiles])
            return reprs[0]  # shape: (768,)
        except Exception as e:
            logger.error(f"Embedding failed for {smiles}: {e}")
            return None

    def get_embeddings_batch(self, smiles_list: list[str]) -> Optional[np.ndarray]:
        """Get embeddings for a batch of molecules."""
        self._ensure_initialized()

        if self._repr_model is None:
            return None

        try:
            reprs = self._repr_model.get_repr(smiles_list)
            return reprs  # shape: (N, 768)
        except Exception as e:
            logger.error(f"Batch embedding failed: {e}")
            return None

    @property
    def is_available(self) -> bool:
        """Check if Uni-Mol model is loaded and ready."""
        self._ensure_initialized()
        return self._repr_model is not None
