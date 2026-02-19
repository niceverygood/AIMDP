"""
Embedding service — generates molecular embeddings using Uni-Mol.
Falls back to Morgan fingerprints + random projection if Uni-Mol is not available.
"""

import logging
import numpy as np
from typing import Optional

logger = logging.getLogger(__name__)

# Dimensionality for all embeddings
EMBEDDING_DIM = 768


class EmbeddingService:
    """
    Molecular embedding service.
    Primary: Uni-Mol (768-dim 3D-aware molecular representations).
    Fallback: Morgan fingerprints projected to 768-dim via random matrix.
    """

    def __init__(self):
        self._unimol_model = None
        self._fallback_projection = None
        self._initialized = False

    def _ensure_initialized(self):
        if self._initialized:
            return

        try:
            from unimol_tools import UniMolRepr
            self._unimol_model = UniMolRepr(model_name="unimol_base")
            logger.info("Uni-Mol model loaded successfully")
        except Exception as e:
            logger.warning(f"Uni-Mol not available, using fallback: {e}")
            # Create a fixed random projection matrix: 2048 → 768
            rng = np.random.RandomState(42)
            self._fallback_projection = rng.randn(2048, EMBEDDING_DIM).astype(np.float32)
            self._fallback_projection /= np.linalg.norm(
                self._fallback_projection, axis=0, keepdims=True
            )

        self._initialized = True

    async def get_embedding(self, smiles: str) -> list[float]:
        """
        Generate a 768-dim embedding for a SMILES string.
        """
        self._ensure_initialized()

        if self._unimol_model is not None:
            return self._get_unimol_embedding(smiles)
        else:
            return self._get_fallback_embedding(smiles)

    async def get_embeddings_batch(self, smiles_list: list[str]) -> list[list[float]]:
        """
        Generate embeddings for a batch of SMILES strings.
        """
        self._ensure_initialized()

        if self._unimol_model is not None:
            return self._get_unimol_embeddings_batch(smiles_list)
        else:
            return [self._get_fallback_embedding(s) for s in smiles_list]

    def _get_unimol_embedding(self, smiles: str) -> list[float]:
        """Get embedding from Uni-Mol model."""
        reprs = self._unimol_model.get_repr([smiles])
        # reprs shape: (1, 768)
        return reprs[0].tolist()

    def _get_unimol_embeddings_batch(self, smiles_list: list[str]) -> list[list[float]]:
        """Get batch embeddings from Uni-Mol."""
        reprs = self._unimol_model.get_repr(smiles_list)
        return [r.tolist() for r in reprs]

    def _get_fallback_embedding(self, smiles: str) -> list[float]:
        """
        Fallback: Morgan fingerprint → random projection → 768-dim.
        This provides a reasonable baseline for similarity search.
        """
        try:
            from rdkit import Chem
            from rdkit.Chem import AllChem

            mol = Chem.MolFromSmiles(smiles)
            if mol is None:
                # Return zero vector for invalid SMILES
                return [0.0] * EMBEDDING_DIM

            # 2048-bit Morgan fingerprint
            fp = AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=2048)
            fp_array = np.array(fp, dtype=np.float32)

            # Project to 768-dim
            embedding = fp_array @ self._fallback_projection
            # L2 normalize
            norm = np.linalg.norm(embedding)
            if norm > 0:
                embedding = embedding / norm

            return embedding.tolist()
        except Exception as e:
            logger.error(f"Fallback embedding failed for {smiles}: {e}")
            return [0.0] * EMBEDDING_DIM
