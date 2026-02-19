"""
Pydantic schemas for request/response validation.
All API responses follow the { success, data, error? } pattern.
"""

from datetime import datetime
from typing import Generic, Optional, TypeVar
from pydantic import BaseModel, Field

T = TypeVar("T")


# ---------- Generic API Response ----------


class APIResponse(BaseModel, Generic[T]):
    """Standard API response wrapper."""
    success: bool
    data: Optional[T] = None
    error: Optional[str] = None


# ---------- Material Schemas ----------


class MaterialBase(BaseModel):
    name: Optional[str] = None
    smiles: str
    category: Optional[str] = None
    molecular_weight: Optional[float] = None
    logp: Optional[float] = None
    hbd: Optional[int] = None
    hba: Optional[int] = None
    tpsa: Optional[float] = None
    rotatable_bonds: Optional[int] = None
    aromatic_rings: Optional[int] = None
    thermal_stability: Optional[float] = None
    dielectric_constant: Optional[float] = None
    bandgap: Optional[float] = None
    solubility: Optional[float] = None
    density: Optional[float] = None
    source: Optional[str] = None


class MaterialCreate(MaterialBase):
    """Schema for creating a new material."""
    pass


class MaterialResponse(MaterialBase):
    """Schema for material in API responses."""
    id: int
    is_verified: bool = False
    created_at: datetime
    updated_at: datetime

    model_config = {"from_attributes": True}


class MaterialSearchResult(MaterialResponse):
    """Material search result with similarity and match scores."""
    similarity: float = Field(default=0.0, description="Vector cosine similarity")
    match_score: float = Field(default=0.0, description="Multi-objective match score")


# ---------- Search Schemas ----------


class SearchQuery(BaseModel):
    """Search query with property filters."""
    smiles: Optional[str] = Field(None, description="Query molecule SMILES")
    category: Optional[str] = Field(None, description="Material category filter")

    # Property range filters
    min_thermal_stability: Optional[float] = Field(None, ge=0, le=500)
    max_thermal_stability: Optional[float] = Field(None, ge=0, le=500)
    min_dielectric: Optional[float] = Field(None, ge=0, le=25)
    max_dielectric: Optional[float] = Field(None, ge=0, le=25)
    min_bandgap: Optional[float] = Field(None, ge=0, le=8)
    max_bandgap: Optional[float] = Field(None, ge=0, le=8)
    min_solubility: Optional[float] = Field(None, ge=0, le=1)
    max_solubility: Optional[float] = Field(None, ge=0, le=1)

    # Weights for multi-objective ranking
    weight_thermal: float = Field(default=0.25, ge=0, le=1)
    weight_dielectric: float = Field(default=0.25, ge=0, le=1)
    weight_bandgap: float = Field(default=0.25, ge=0, le=1)
    weight_similarity: float = Field(default=0.25, ge=0, le=1)

    limit: int = Field(default=20, ge=1, le=100)
    offset: int = Field(default=0, ge=0)


class SearchResponse(BaseModel):
    """Search results."""
    results: list[MaterialSearchResult]
    total: int
    query_time_ms: float


# ---------- Pipeline Schemas ----------


class PipelineStepStatus(BaseModel):
    """Status of a single pipeline step."""
    name: str
    label: str  # Korean UI label
    status: str = Field(
        default="pending",
        description="pending | running | completed | failed",
    )
    progress: float = Field(default=0.0, ge=0, le=100)
    rows_processed: int = 0
    total_rows: int = 0
    quality_score: Optional[float] = None
    time_elapsed_sec: float = 0.0
    error: Optional[str] = None


class PipelineStatus(BaseModel):
    """Overall pipeline status."""
    pipeline_id: str
    status: str  # pending | running | completed | failed
    steps: list[PipelineStepStatus]
    started_at: Optional[datetime] = None
    completed_at: Optional[datetime] = None


class PipelineStartRequest(BaseModel):
    """Request to start a new pipeline run."""
    source_type: str = Field(description="csv | sdf | mol | json")
    file_path: str
    options: dict = Field(default_factory=dict)


# ---------- Property Prediction Schemas ----------


class PropertyPrediction(BaseModel):
    """Predicted properties for a molecule."""
    smiles: str
    thermal_stability: Optional[float] = None
    dielectric_constant: Optional[float] = None
    bandgap: Optional[float] = None
    solubility: Optional[float] = None
    density: Optional[float] = None
    confidence: Optional[float] = None  # model confidence score


# ---------- Experiment / Feedback Schemas ----------


class ExperimentResult(BaseModel):
    """Experiment result submitted by a researcher."""
    material_id: int
    actual_thermal_stability: Optional[float] = None
    actual_dielectric_constant: Optional[float] = None
    actual_bandgap: Optional[float] = None
    actual_solubility: Optional[float] = None
    actual_density: Optional[float] = None
    success: bool = False
    notes: Optional[str] = None
    researcher: Optional[str] = None


class FeedbackAnalysis(BaseModel):
    """Analysis of prediction vs experiment."""
    material_id: int
    prediction_errors: dict[str, float]  # property -> error %
    model_improvement: str  # "재학습 예정" or "정상 범위"
    retrain_triggered: bool = False
