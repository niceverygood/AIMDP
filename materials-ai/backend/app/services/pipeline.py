"""
Data pipeline orchestrator.
Transforms raw data (CSV, SDF, MOL, JSON) into AI-ready training datasets.

Pipeline steps:
  1. Collect → 2. Normalize → 3. Clean → 4. Feature Engineer → 5. Embed → 6. Dataset
"""

import asyncio
import logging
import time
from typing import Callable, Optional

logger = logging.getLogger(__name__)

ProgressCallback = Callable[[str, str, dict], None]


class PipelineOrchestrator:
    """
    Orchestrates the full data processing pipeline.
    Reports progress via callback (WebSocket-friendly).
    """

    async def run_pipeline(
        self,
        pipeline_id: str,
        source_type: str,
        file_path: str,
        options: dict,
        progress_callback: Optional[ProgressCallback] = None,
    ):
        """Execute the full pipeline with progress reporting."""
        steps = [
            ("collect", self._step_collect),
            ("normalize", self._step_normalize),
            ("clean", self._step_clean),
            ("engineer", self._step_feature_engineering),
            ("embed", self._step_generate_embeddings),
            ("dataset", self._step_create_dataset),
        ]

        context = {
            "source_type": source_type,
            "file_path": file_path,
            "options": options,
            "data": None,
            "rows": 0,
        }

        for step_name, step_func in steps:
            start = time.time()
            try:
                if progress_callback:
                    await progress_callback(pipeline_id, step_name, {
                        "status": "running",
                        "progress": 0,
                    })

                context = await step_func(context, pipeline_id, progress_callback)

                elapsed = time.time() - start
                if progress_callback:
                    await progress_callback(pipeline_id, step_name, {
                        "status": "completed",
                        "progress": 100,
                        "rows_processed": context.get("rows", 0),
                        "total_rows": context.get("total_rows", 0),
                        "quality_score": context.get("quality_score"),
                        "time_elapsed_sec": round(elapsed, 2),
                    })

            except Exception as e:
                elapsed = time.time() - start
                logger.error(f"Pipeline step {step_name} failed: {e}")
                if progress_callback:
                    await progress_callback(pipeline_id, step_name, {
                        "status": "failed",
                        "error": str(e),
                        "time_elapsed_sec": round(elapsed, 2),
                    })
                return  # Stop pipeline on failure

        logger.info(f"Pipeline {pipeline_id} completed successfully")

    async def _step_collect(
        self, context: dict, pipeline_id: str, callback: Optional[ProgressCallback]
    ) -> dict:
        """Step 1: Collect and load raw data from source files."""
        source_type = context["source_type"]
        file_path = context["file_path"]

        if source_type == "csv":
            from data.loaders.csv_loader import CSVLoader
            loader = CSVLoader()
        elif source_type in ("sdf", "mol"):
            from data.loaders.sdf_loader import SDFLoader
            loader = SDFLoader()
        else:
            raise ValueError(f"Unsupported source type: {source_type}")

        data = await loader.load(file_path)
        context["data"] = data
        context["rows"] = len(data)
        context["total_rows"] = len(data)

        if callback:
            await callback(pipeline_id, "collect", {
                "status": "running",
                "progress": 100,
                "rows_processed": len(data),
                "total_rows": len(data),
            })

        return context

    async def _step_normalize(
        self, context: dict, pipeline_id: str, callback: Optional[ProgressCallback]
    ) -> dict:
        """Step 2: Normalize data (unit conversion, format standardization)."""
        from data.processors.normalizer import DataNormalizer

        normalizer = DataNormalizer()
        data = context["data"]
        total = len(data)

        normalized = []
        for i, record in enumerate(data):
            normalized.append(normalizer.normalize(record))
            if callback and i % max(1, total // 10) == 0:
                await callback(pipeline_id, "normalize", {
                    "status": "running",
                    "progress": round((i + 1) / total * 100, 1),
                    "rows_processed": i + 1,
                    "total_rows": total,
                })

        context["data"] = normalized
        return context

    async def _step_clean(
        self, context: dict, pipeline_id: str, callback: Optional[ProgressCallback]
    ) -> dict:
        """Step 3: Clean data (remove duplicates, handle missing values)."""
        data = context["data"]
        original_count = len(data)

        # Remove invalid entries
        cleaned = [d for d in data if d.get("smiles") and d["smiles"].strip()]

        # Remove duplicates by SMILES
        seen = set()
        unique = []
        for d in cleaned:
            if d["smiles"] not in seen:
                seen.add(d["smiles"])
                unique.append(d)

        context["data"] = unique
        context["rows"] = len(unique)
        context["quality_score"] = round(len(unique) / max(1, original_count) * 100, 1)
        return context

    async def _step_feature_engineering(
        self, context: dict, pipeline_id: str, callback: Optional[ProgressCallback]
    ) -> dict:
        """Step 4: Extract molecular descriptors using RDKit."""
        from data.processors.feature_engineer import FeatureEngineer

        engineer = FeatureEngineer()
        data = context["data"]
        total = len(data)

        enriched = []
        for i, record in enumerate(data):
            features = engineer.compute_features(record.get("smiles", ""))
            if features:
                record.update(features)
                enriched.append(record)

            if callback and i % max(1, total // 10) == 0:
                await callback(pipeline_id, "engineer", {
                    "status": "running",
                    "progress": round((i + 1) / total * 100, 1),
                    "rows_processed": i + 1,
                    "total_rows": total,
                })

        context["data"] = enriched
        context["rows"] = len(enriched)
        return context

    async def _step_generate_embeddings(
        self, context: dict, pipeline_id: str, callback: Optional[ProgressCallback]
    ) -> dict:
        """Step 5: Generate Uni-Mol embeddings for all molecules."""
        from app.services.embedding import EmbeddingService

        embedding_service = EmbeddingService()
        data = context["data"]
        total = len(data)
        batch_size = 32

        for i in range(0, total, batch_size):
            batch = data[i : i + batch_size]
            smiles_batch = [d["smiles"] for d in batch]
            embeddings = await embedding_service.get_embeddings_batch(smiles_batch)

            for j, emb in enumerate(embeddings):
                data[i + j]["embedding"] = emb

            if callback:
                processed = min(i + batch_size, total)
                await callback(pipeline_id, "embed", {
                    "status": "running",
                    "progress": round(processed / total * 100, 1),
                    "rows_processed": processed,
                    "total_rows": total,
                })

        context["data"] = data
        return context

    async def _step_create_dataset(
        self, context: dict, pipeline_id: str, callback: Optional[ProgressCallback]
    ) -> dict:
        """Step 6: Insert processed data into the database."""
        # In a real implementation, this would batch-insert into PostgreSQL
        data = context["data"]
        context["rows"] = len(data)
        logger.info(f"Dataset created with {len(data)} records")
        return context
