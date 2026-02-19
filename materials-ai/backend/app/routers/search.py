"""
Search API — AI-powered material discovery search.
Vector similarity + multi-objective ranking.
"""

import time
from fastapi import APIRouter, Depends, HTTPException
from sqlalchemy.ext.asyncio import AsyncSession
from sqlalchemy import text, func, select, and_

from app.core.database import get_db
from app.models.schemas import (
    APIResponse,
    SearchQuery,
    SearchResponse,
    MaterialSearchResult,
    PropertyPrediction,
)
from app.services.embedding import EmbeddingService
from app.services.predictor import PropertyPredictor
from app.services.ranker import MultiObjectiveRanker
from app.models.material import Material

router = APIRouter()

embedding_service = EmbeddingService()
predictor = PropertyPredictor()
ranker = MultiObjectiveRanker()


@router.post("", response_model=APIResponse[SearchResponse])
async def search_materials(
    query: SearchQuery,
    db: AsyncSession = Depends(get_db),
):
    """
    AI 기반 소재 탐색 — 벡터 유사도 + 물성 필터링 + 다목적 최적화 랭킹.
    """
    start_time = time.time()

    try:
        # Build dynamic WHERE clause
        conditions = []
        params: dict = {
            "limit": query.limit,
            "offset": query.offset,
        }

        if query.min_thermal_stability is not None:
            conditions.append("m.thermal_stability >= :min_thermal")
            params["min_thermal"] = query.min_thermal_stability
        if query.max_thermal_stability is not None:
            conditions.append("m.thermal_stability <= :max_thermal")
            params["max_thermal"] = query.max_thermal_stability
        if query.min_dielectric is not None:
            conditions.append("m.dielectric_constant >= :min_dielectric")
            params["min_dielectric"] = query.min_dielectric
        if query.max_dielectric is not None:
            conditions.append("m.dielectric_constant <= :max_dielectric")
            params["max_dielectric"] = query.max_dielectric
        if query.min_bandgap is not None:
            conditions.append("m.bandgap >= :min_bandgap")
            params["min_bandgap"] = query.min_bandgap
        if query.max_bandgap is not None:
            conditions.append("m.bandgap <= :max_bandgap")
            params["max_bandgap"] = query.max_bandgap
        if query.min_solubility is not None:
            conditions.append("m.solubility >= :min_solubility")
            params["min_solubility"] = query.min_solubility
        if query.max_solubility is not None:
            conditions.append("m.solubility <= :max_solubility")
            params["max_solubility"] = query.max_solubility
        if query.category:
            conditions.append("m.category = :category")
            params["category"] = query.category

        where_clause = ""
        if conditions:
            where_clause = "WHERE " + " AND ".join(conditions)

        # If SMILES provided, use vector similarity search
        if query.smiles:
            query_embedding = await embedding_service.get_embedding(query.smiles)
            params["query_vec"] = str(query_embedding)

            sql = f"""
                SELECT m.*,
                       1 - (m.embedding <=> :query_vec::vector) as similarity
                FROM materials m
                {where_clause}
                ORDER BY m.embedding <=> :query_vec::vector
                LIMIT :limit OFFSET :offset
            """
        else:
            # No SMILES query — sort by thermal stability descending as default
            sql = f"""
                SELECT m.*, 0.0 as similarity
                FROM materials m
                {where_clause}
                ORDER BY m.thermal_stability DESC NULLS LAST
                LIMIT :limit OFFSET :offset
            """

        result = await db.execute(text(sql), params)
        rows = result.mappings().all()

        # Count total matches
        count_sql = f"SELECT COUNT(*) FROM materials m {where_clause}"
        count_params = {k: v for k, v in params.items() if k not in ("limit", "offset", "query_vec")}
        total_result = await db.execute(text(count_sql), count_params)
        total = total_result.scalar()

        # Build response with multi-objective scores
        results = []
        for row in rows:
            match_score = ranker.compute_score(
                similarity=row["similarity"],
                thermal_stability=row.get("thermal_stability"),
                dielectric_constant=row.get("dielectric_constant"),
                bandgap=row.get("bandgap"),
                weights={
                    "similarity": query.weight_similarity,
                    "thermal": query.weight_thermal,
                    "dielectric": query.weight_dielectric,
                    "bandgap": query.weight_bandgap,
                },
            )

            results.append(
                MaterialSearchResult(
                    id=row["id"],
                    name=row.get("name"),
                    smiles=row["smiles"],
                    category=row.get("category"),
                    molecular_weight=row.get("molecular_weight"),
                    logp=row.get("logp"),
                    hbd=row.get("hbd"),
                    hba=row.get("hba"),
                    tpsa=row.get("tpsa"),
                    rotatable_bonds=row.get("rotatable_bonds"),
                    aromatic_rings=row.get("aromatic_rings"),
                    thermal_stability=row.get("thermal_stability"),
                    dielectric_constant=row.get("dielectric_constant"),
                    bandgap=row.get("bandgap"),
                    solubility=row.get("solubility"),
                    density=row.get("density"),
                    source=row.get("source"),
                    is_verified=row.get("is_verified", False),
                    created_at=row["created_at"],
                    updated_at=row["updated_at"],
                    similarity=row["similarity"],
                    match_score=match_score,
                )
            )

        # Sort by match_score descending
        results.sort(key=lambda r: r.match_score, reverse=True)

        elapsed = (time.time() - start_time) * 1000

        return APIResponse(
            success=True,
            data=SearchResponse(
                results=results,
                total=total,
                query_time_ms=round(elapsed, 2),
            ),
        )
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))


@router.post("/predict", response_model=APIResponse[PropertyPrediction])
async def predict_properties(smiles: str):
    """
    주어진 SMILES 분자 구조에 대한 물성 예측.
    """
    try:
        prediction = await predictor.predict(smiles)
        return APIResponse(success=True, data=prediction)
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))
