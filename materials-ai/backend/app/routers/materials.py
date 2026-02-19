"""
Materials CRUD API — manage materials in the database.
Includes experiment feedback endpoint for active learning.
"""

from fastapi import APIRouter, Depends, HTTPException, Query
from sqlalchemy.ext.asyncio import AsyncSession
from sqlalchemy import select, func

from app.core.database import get_db
from app.models.material import Material, Experiment
from app.models.schemas import (
    APIResponse,
    MaterialCreate,
    MaterialResponse,
    ExperimentResult,
    FeedbackAnalysis,
)
from app.services.embedding import EmbeddingService

router = APIRouter()
embedding_service = EmbeddingService()


@router.get("", response_model=APIResponse[list[MaterialResponse]])
async def list_materials(
    category: str | None = None,
    limit: int = Query(default=20, ge=1, le=100),
    offset: int = Query(default=0, ge=0),
    db: AsyncSession = Depends(get_db),
):
    """소재 목록 조회 (카테고리 필터 지원)."""
    query = select(Material)
    if category:
        query = query.where(Material.category == category)
    query = query.order_by(Material.created_at.desc()).limit(limit).offset(offset)

    result = await db.execute(query)
    materials = result.scalars().all()

    return APIResponse(
        success=True,
        data=[MaterialResponse.model_validate(m) for m in materials],
    )


@router.get("/stats")
async def material_stats(db: AsyncSession = Depends(get_db)):
    """소재 DB 통계 요약."""
    total = await db.scalar(select(func.count(Material.id)))
    verified = await db.scalar(
        select(func.count(Material.id)).where(Material.is_verified.is_(True))
    )
    categories_result = await db.execute(
        select(Material.category, func.count(Material.id))
        .group_by(Material.category)
        .order_by(func.count(Material.id).desc())
    )
    categories = {row[0] or "미분류": row[1] for row in categories_result.all()}

    return APIResponse(
        success=True,
        data={
            "total_materials": total,
            "verified_materials": verified,
            "categories": categories,
        },
    )


@router.get("/{material_id}", response_model=APIResponse[MaterialResponse])
async def get_material(
    material_id: int,
    db: AsyncSession = Depends(get_db),
):
    """특정 소재 상세 조회."""
    material = await db.get(Material, material_id)
    if not material:
        raise HTTPException(status_code=404, detail="Material not found")
    return APIResponse(
        success=True,
        data=MaterialResponse.model_validate(material),
    )


@router.post("", response_model=APIResponse[MaterialResponse])
async def create_material(
    data: MaterialCreate,
    db: AsyncSession = Depends(get_db),
):
    """새 소재 등록 (자동으로 임베딩 생성)."""
    try:
        # Generate embedding
        embedding = await embedding_service.get_embedding(data.smiles)

        material = Material(
            **data.model_dump(),
            embedding=embedding,
        )
        db.add(material)
        await db.flush()
        await db.refresh(material)

        return APIResponse(
            success=True,
            data=MaterialResponse.model_validate(material),
        )
    except Exception as e:
        raise HTTPException(status_code=400, detail=str(e))


@router.delete("/{material_id}")
async def delete_material(
    material_id: int,
    db: AsyncSession = Depends(get_db),
):
    """소재 삭제."""
    material = await db.get(Material, material_id)
    if not material:
        raise HTTPException(status_code=404, detail="Material not found")
    await db.delete(material)
    return APIResponse(success=True, data={"deleted": material_id})


@router.post("/feedback", response_model=APIResponse[FeedbackAnalysis])
async def submit_experiment_feedback(
    result: ExperimentResult,
    db: AsyncSession = Depends(get_db),
):
    """
    실험 결과 피드백 제출 — Active Learning 루프.
    예측값 vs 실측값을 비교하고, 오차가 크면 재학습 트리거.
    """
    material = await db.get(Material, result.material_id)
    if not material:
        raise HTTPException(status_code=404, detail="Material not found")

    # Save experiment
    experiment = Experiment(
        material_id=result.material_id,
        actual_thermal_stability=result.actual_thermal_stability,
        actual_dielectric_constant=result.actual_dielectric_constant,
        actual_bandgap=result.actual_bandgap,
        actual_solubility=result.actual_solubility,
        actual_density=result.actual_density,
        success=result.success,
        notes=result.notes,
        researcher=result.researcher,
    )
    db.add(experiment)
    await db.flush()

    # Calculate prediction errors
    errors = {}
    retrain_threshold = 0.15  # 15% error threshold

    property_pairs = [
        ("thermal_stability", material.thermal_stability, result.actual_thermal_stability),
        ("dielectric_constant", material.dielectric_constant, result.actual_dielectric_constant),
        ("bandgap", material.bandgap, result.actual_bandgap),
        ("solubility", material.solubility, result.actual_solubility),
        ("density", material.density, result.actual_density),
    ]

    max_error = 0.0
    for prop_name, predicted, actual in property_pairs:
        if predicted is not None and actual is not None and predicted != 0:
            error = abs(predicted - actual) / abs(predicted)
            errors[prop_name] = round(error, 4)
            max_error = max(max_error, error)

    retrain_triggered = max_error > retrain_threshold

    return APIResponse(
        success=True,
        data=FeedbackAnalysis(
            material_id=result.material_id,
            prediction_errors=errors,
            model_improvement="재학습 예정" if retrain_triggered else "정상 범위",
            retrain_triggered=retrain_triggered,
        ),
    )


@router.get("/{material_id}/sdf")
async def get_material_sdf(
    material_id: int,
    db: AsyncSession = Depends(get_db),
):
    """소재의 3D 구조를 SDF 형식으로 반환 (3Dmol.js 뷰어용)."""
    material = await db.get(Material, material_id)
    if not material:
        raise HTTPException(status_code=404, detail="Material not found")

    try:
        from rdkit import Chem
        from rdkit.Chem import AllChem

        mol = Chem.MolFromSmiles(material.smiles)
        if mol is None:
            raise ValueError("Invalid SMILES")

        mol = Chem.AddHs(mol)
        AllChem.EmbedMolecule(mol, randomSeed=42)
        AllChem.MMFFOptimizeMolecule(mol)

        from fastapi.responses import PlainTextResponse

        return PlainTextResponse(
            content=Chem.MolToMolBlock(mol),
            media_type="chemical/x-mdl-molfile",
        )
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"SDF generation failed: {e}")
