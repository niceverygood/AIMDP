"""
Lightweight Materials AI API server.
Uses SQLite (no PostgreSQL needed) with brute-force cosine similarity.
This is the PoC version that actually works with real data.

Run: python server.py
"""

import io
import json
import math
import os
import sqlite3
import tempfile
import time
from typing import Optional

from dotenv import load_dotenv
load_dotenv()

import numpy as np
from fastapi import FastAPI, HTTPException, Query, UploadFile, File, Form
from fastapi.middleware.cors import CORSMiddleware
from fastapi.responses import PlainTextResponse
from pydantic import BaseModel, Field

# ─── Configuration ───────────────────────────────
DB_PATH = os.environ.get("DB_PATH", "./data/materials.db")
EMBEDDING_DIM = 768

app = FastAPI(
    title="Materials AI Platform",
    version="0.2.0 (PoC)",
    description="AI-powered materials discovery — with real data",
)

app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)


# ─── Helpers ─────────────────────────────────────

def get_db():
    conn = sqlite3.connect(DB_PATH)
    conn.row_factory = sqlite3.Row
    return conn


def cosine_similarity(a: list[float], b: list[float]) -> float:
    """Cosine similarity between two vectors."""
    a_arr = np.array(a, dtype=np.float32)
    b_arr = np.array(b, dtype=np.float32)
    dot = np.dot(a_arr, b_arr)
    norm_a = np.linalg.norm(a_arr)
    norm_b = np.linalg.norm(b_arr)
    if norm_a == 0 or norm_b == 0:
        return 0.0
    return float(dot / (norm_a * norm_b))


def get_embedding_for_smiles(smiles: str) -> list[float]:
    """Generate embedding for a query SMILES."""
    try:
        from rdkit import Chem
        from rdkit.Chem import AllChem

        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return [0.0] * EMBEDDING_DIM

        fp = AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=2048)
        fp_array = np.array(fp, dtype=np.float32)

        rng = np.random.RandomState(42)
        projection = rng.randn(2048, EMBEDDING_DIM).astype(np.float32)
        projection /= np.linalg.norm(projection, axis=0, keepdims=True)

        embedding = fp_array @ projection
        norm = np.linalg.norm(embedding)
        if norm > 0:
            embedding = embedding / norm
        return embedding.tolist()
    except Exception:
        return [0.0] * EMBEDDING_DIM


def predict_properties_heuristic(smiles: str) -> dict:
    """Predict properties using RDKit descriptors."""
    try:
        from rdkit import Chem
        from rdkit.Chem import Descriptors

        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return {}

        mw = Descriptors.MolWt(mol)
        logp = Descriptors.MolLogP(mol)
        tpsa = Descriptors.TPSA(mol)
        n_aromatic = Descriptors.NumAromaticRings(mol)
        hba = Descriptors.NumHAcceptors(mol)

        return {
            "thermal_stability": round(min(500, max(80, 120 + n_aromatic * 45 + mw * 0.15 - abs(logp) * 3)), 1),
            "dielectric_constant": round(max(1.5, min(25, 2.5 + tpsa * 0.04 + hba * 0.3 - logp * 0.15)), 2),
            "bandgap": round(max(0.5, min(8, 4.5 - n_aromatic * 0.4 - logp * 0.08 + tpsa * 0.005)), 2),
            "solubility": round(max(0.0001, min(1.0, 0.5 - logp * 0.08 + tpsa * 0.004)), 4),
            "density": round(max(0.6, min(2.5, 0.85 + mw * 0.0015 + n_aromatic * 0.03)), 3),
            "confidence": 0.4,
        }
    except Exception:
        return {}


# ─── Schemas ─────────────────────────────────────

class SearchQuery(BaseModel):
    smiles: Optional[str] = None
    category: Optional[str] = None
    min_thermal_stability: Optional[float] = None
    max_thermal_stability: Optional[float] = None
    min_dielectric: Optional[float] = None
    max_dielectric: Optional[float] = None
    min_bandgap: Optional[float] = None
    max_bandgap: Optional[float] = None
    min_solubility: Optional[float] = None
    max_solubility: Optional[float] = None
    weight_thermal: float = 0.25
    weight_dielectric: float = 0.25
    weight_bandgap: float = 0.25
    weight_similarity: float = 0.25
    limit: int = Field(default=20, ge=1, le=100)
    offset: int = Field(default=0, ge=0)


class ExperimentResult(BaseModel):
    material_id: int
    actual_thermal_stability: Optional[float] = None
    actual_dielectric_constant: Optional[float] = None
    actual_bandgap: Optional[float] = None
    actual_solubility: Optional[float] = None
    actual_density: Optional[float] = None
    success: bool = False
    notes: Optional[str] = None
    researcher: Optional[str] = None


# ─── API Endpoints ───────────────────────────────

@app.get("/health")
def health():
    return {"status": "healthy", "version": "0.2.0", "db": DB_PATH}


@app.get("/")
def root():
    return {"message": "Materials AI Platform API (PoC)", "docs": "/docs"}


# --- Search ---

@app.post("/api/search")
def search_materials(query: SearchQuery):
    start = time.time()
    conn = get_db()

    # Build WHERE clause
    conditions = []
    params = []

    if query.category and query.category != "all":
        conditions.append("LOWER(category) = LOWER(?)")
        params.append(query.category)
    if query.min_thermal_stability is not None:
        conditions.append("thermal_stability >= ?")
        params.append(query.min_thermal_stability)
    if query.max_thermal_stability is not None:
        conditions.append("thermal_stability <= ?")
        params.append(query.max_thermal_stability)
    if query.min_dielectric is not None:
        conditions.append("dielectric_constant >= ?")
        params.append(query.min_dielectric)
    if query.max_dielectric is not None:
        conditions.append("dielectric_constant <= ?")
        params.append(query.max_dielectric)
    if query.min_bandgap is not None:
        conditions.append("bandgap >= ?")
        params.append(query.min_bandgap)
    if query.max_bandgap is not None:
        conditions.append("bandgap <= ?")
        params.append(query.max_bandgap)
    if query.min_solubility is not None:
        conditions.append("solubility >= ?")
        params.append(query.min_solubility)
    if query.max_solubility is not None:
        conditions.append("solubility <= ?")
        params.append(query.max_solubility)

    where = ""
    if conditions:
        where = "WHERE " + " AND ".join(conditions)

    # Get all matching materials
    sql = f"SELECT * FROM materials {where}"
    rows = conn.execute(sql, params).fetchall()

    # Count total
    count_sql = f"SELECT COUNT(*) FROM materials {where}"
    total = conn.execute(count_sql, params).fetchone()[0]

    # Compute query embedding if SMILES provided
    query_embedding = None
    if query.smiles:
        query_embedding = get_embedding_for_smiles(query.smiles)

    # Score and rank
    results = []
    for row in rows:
        row_dict = dict(row)
        embedding_json = row_dict.pop("embedding", "[]")

        # Compute similarity
        similarity = 0.0
        if query_embedding:
            try:
                row_embedding = json.loads(embedding_json)
                if row_embedding and len(row_embedding) == EMBEDDING_DIM:
                    similarity = cosine_similarity(query_embedding, row_embedding)
            except (json.JSONDecodeError, TypeError):
                pass

        # Multi-objective score
        score = 0.0
        weight_sum = 0.0

        # Similarity component
        w_sim = query.weight_similarity
        if w_sim > 0 and query_embedding:
            score += w_sim * max(0, similarity)
            weight_sum += w_sim

        # Property scores (normalized to 0-1)
        ranges = {
            "thermal": ("thermal_stability", 0, 500, query.weight_thermal),
            "dielectric": ("dielectric_constant", 0, 25, query.weight_dielectric),
            "bandgap": ("bandgap", 0, 8, query.weight_bandgap),
        }
        for key, (col, low, high, weight) in ranges.items():
            val = row_dict.get(col)
            if weight > 0 and val is not None:
                normalized = (val - low) / (high - low) if high != low else 0
                score += weight * max(0, min(1, normalized))
                weight_sum += weight

        match_score = round(score / weight_sum, 4) if weight_sum > 0 else 0.0

        results.append({
            **row_dict,
            "similarity": round(similarity, 4),
            "match_score": match_score,
        })

    # Sort by match_score
    results.sort(key=lambda r: r["match_score"], reverse=True)

    # Paginate
    paginated = results[query.offset : query.offset + query.limit]

    elapsed = (time.time() - start) * 1000
    conn.close()

    return {
        "success": True,
        "data": {
            "results": paginated,
            "total": total,
            "query_time_ms": round(elapsed, 2),
        },
    }


@app.post("/api/search/predict")
def predict_properties(smiles: str):
    props = predict_properties_heuristic(smiles)
    if not props:
        raise HTTPException(400, "Invalid SMILES")
    return {
        "success": True,
        "data": {"smiles": smiles, **props},
    }


# --- Materials CRUD ---

@app.get("/api/materials")
def list_materials(
    category: Optional[str] = None,
    limit: int = Query(default=20, ge=1, le=100),
    offset: int = Query(default=0, ge=0),
):
    conn = get_db()
    if category and category != "all":
        rows = conn.execute(
            "SELECT * FROM materials WHERE category = ? ORDER BY id DESC LIMIT ? OFFSET ?",
            (category, limit, offset),
        ).fetchall()
    else:
        rows = conn.execute(
            "SELECT * FROM materials ORDER BY id DESC LIMIT ? OFFSET ?",
            (limit, offset),
        ).fetchall()

    results = []
    for row in rows:
        d = dict(row)
        d.pop("embedding", None)
        results.append(d)

    conn.close()
    return {"success": True, "data": results}


@app.get("/api/materials/stats")
def material_stats():
    conn = get_db()
    total = conn.execute("SELECT COUNT(*) FROM materials").fetchone()[0]
    verified = conn.execute("SELECT COUNT(*) FROM materials WHERE is_verified = 1").fetchone()[0]
    categories = conn.execute(
        "SELECT category, COUNT(*) as cnt FROM materials GROUP BY category ORDER BY cnt DESC"
    ).fetchall()

    conn.close()
    return {
        "success": True,
        "data": {
            "total_materials": total,
            "verified_materials": verified,
            "categories": {row["category"] or "미분류": row["cnt"] for row in categories},
        },
    }


@app.get("/api/materials/{material_id}")
def get_material(material_id: int):
    conn = get_db()
    row = conn.execute("SELECT * FROM materials WHERE id = ?", (material_id,)).fetchone()
    conn.close()
    if not row:
        raise HTTPException(404, "Material not found")
    d = dict(row)
    d.pop("embedding", None)
    return {"success": True, "data": d}


@app.get("/api/materials/{material_id}/sdf")
def get_material_sdf(material_id: int):
    conn = get_db()
    row = conn.execute("SELECT smiles FROM materials WHERE id = ?", (material_id,)).fetchone()
    conn.close()
    if not row:
        raise HTTPException(404, "Material not found")
    try:
        from rdkit import Chem
        from rdkit.Chem import AllChem

        mol = Chem.MolFromSmiles(row["smiles"])
        if mol is None:
            raise ValueError("Invalid SMILES")
        mol = Chem.AddHs(mol)
        AllChem.EmbedMolecule(mol, randomSeed=42)
        AllChem.MMFFOptimizeMolecule(mol)
        return PlainTextResponse(content=Chem.MolToMolBlock(mol), media_type="chemical/x-mdl-molfile")
    except Exception as e:
        raise HTTPException(500, f"SDF generation failed: {e}")


# --- Molecule SDF by SMILES (for 3D viewer) ---

@app.get("/api/molecule/sdf")
def get_molecule_sdf(smiles: str):
    try:
        from rdkit import Chem
        from rdkit.Chem import AllChem

        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            raise ValueError("Invalid SMILES")
        mol = Chem.AddHs(mol)
        AllChem.EmbedMolecule(mol, randomSeed=42)
        AllChem.MMFFOptimizeMolecule(mol)
        return PlainTextResponse(content=Chem.MolToMolBlock(mol), media_type="chemical/x-mdl-molfile")
    except Exception as e:
        raise HTTPException(500, f"SDF generation failed: {e}")


# --- Feedback ---

@app.post("/api/materials/feedback")
def submit_feedback(result: ExperimentResult):
    conn = get_db()
    row = conn.execute("SELECT * FROM materials WHERE id = ?", (result.material_id,)).fetchone()
    if not row:
        conn.close()
        raise HTTPException(404, "Material not found")

    # Save experiment
    conn.execute(
        """INSERT INTO experiments
           (material_id, actual_thermal_stability, actual_dielectric_constant,
            actual_bandgap, actual_solubility, actual_density, success, notes, researcher)
           VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)""",
        (
            result.material_id,
            result.actual_thermal_stability,
            result.actual_dielectric_constant,
            result.actual_bandgap,
            result.actual_solubility,
            result.actual_density,
            1 if result.success else 0,
            result.notes,
            result.researcher,
        ),
    )
    conn.commit()

    # Calculate errors
    errors = {}
    threshold = 0.15
    max_error = 0.0
    pairs = [
        ("thermal_stability", row["thermal_stability"], result.actual_thermal_stability),
        ("dielectric_constant", row["dielectric_constant"], result.actual_dielectric_constant),
        ("bandgap", row["bandgap"], result.actual_bandgap),
        ("solubility", row["solubility"], result.actual_solubility),
        ("density", row["density"], result.actual_density),
    ]
    for name, predicted, actual in pairs:
        if predicted and actual and predicted != 0:
            err = abs(predicted - actual) / abs(predicted)
            errors[name] = round(err, 4)
            max_error = max(max_error, err)

    conn.close()
    retrain = max_error > threshold

    return {
        "success": True,
        "data": {
            "material_id": result.material_id,
            "prediction_errors": errors,
            "model_improvement": "재학습 예정" if retrain else "정상 범위",
            "retrain_triggered": retrain,
        },
    }


# --- AI Ensemble Analysis ---

class AIAnalysisRequest(BaseModel):
    smiles: Optional[str] = None
    category: Optional[str] = None
    min_thermal_stability: Optional[float] = None
    min_bandgap: Optional[float] = None
    min_dielectric: Optional[float] = None
    target_application: Optional[str] = None
    limit: int = Field(default=10, ge=1, le=20)


@app.post("/api/ai/analyze")
async def ai_analyze_materials(request: AIAnalysisRequest):
    """
    3개 AI 모델(Claude Opus 4, Gemini 2.5 Pro, GPT-4o) 앙상블 분석.
    OpenRouter를 통해 병렬 호출 → 교차 검증 → 합의 기반 추천.
    """
    from app.services.ai_ensemble import analyze_materials

    # Fetch candidates from DB
    conn = get_db()
    conditions = []
    params = []

    if request.category and request.category != "all":
        conditions.append("LOWER(category) = LOWER(?)")
        params.append(request.category)
    if request.min_thermal_stability is not None:
        conditions.append("thermal_stability >= ?")
        params.append(request.min_thermal_stability)
    if request.min_bandgap is not None:
        conditions.append("bandgap >= ?")
        params.append(request.min_bandgap)
    if request.min_dielectric is not None:
        conditions.append("dielectric_constant >= ?")
        params.append(request.min_dielectric)

    where = "WHERE " + " AND ".join(conditions) if conditions else ""
    rows = conn.execute(
        f"SELECT * FROM materials {where} ORDER BY thermal_stability DESC LIMIT ?",
        (*params, request.limit),
    ).fetchall()
    conn.close()

    if not rows:
        return {"success": False, "error": "조건에 맞는 소재가 없습니다."}

    candidates = []
    for row in rows:
        d = dict(row)
        d.pop("embedding", None)
        candidates.append(d)

    requirements = {
        "smiles": request.smiles,
        "category": request.category,
        "min_thermal_stability": request.min_thermal_stability,
        "min_bandgap": request.min_bandgap,
        "min_dielectric": request.min_dielectric,
        "target_application": request.target_application,
    }

    result = await analyze_materials(candidates, requirements)

    return {"success": True, "data": result}


# --- Research Data Upload ---

@app.post("/api/data/upload")
async def upload_research_data(
    file: UploadFile = File(...),
    category: str = Form(default=""),
    source: str = Form(default="dongjin_internal"),
):
    """
    연구 데이터 파일 업로드 (CSV/Excel).
    파일을 파싱하여 소재 DB에 자동 적재합니다.
    """
    if not file.filename:
        raise HTTPException(400, "파일이 없습니다")

    ext = file.filename.rsplit(".", 1)[-1].lower()
    if ext not in ("csv", "xlsx", "xls", "tsv"):
        raise HTTPException(400, f"지원하지 않는 형식: .{ext} (CSV, Excel만 지원)")

    try:
        import pandas as pd

        content = await file.read()

        if ext in ("xlsx", "xls"):
            df = pd.read_excel(io.BytesIO(content))
        elif ext == "tsv":
            df = pd.read_csv(io.BytesIO(content), sep="\t")
        else:
            df = pd.read_csv(io.BytesIO(content))

        # Auto-detect column names
        col_map = {}
        for col in df.columns:
            cl = col.strip().lower()
            if cl in ("smiles", "smi"): col_map[col] = "smiles"
            elif cl in ("name", "이름", "소재명", "material"): col_map[col] = "name"
            elif cl in ("category", "카테고리", "분류"): col_map[col] = "category"
            elif "thermal" in cl or "td" == cl or "열안정" in cl: col_map[col] = "thermal_stability"
            elif "dielectric" in cl or "유전" in cl or "dk" == cl: col_map[col] = "dielectric_constant"
            elif "bandgap" in cl or "밴드갭" in cl or "eg" == cl: col_map[col] = "bandgap"
            elif "solub" in cl or "용해" in cl: col_map[col] = "solubility"
            elif "density" in cl or "밀도" in cl: col_map[col] = "density"

        df = df.rename(columns=col_map)

        if "smiles" not in df.columns:
            raise HTTPException(400, "SMILES 컬럼을 찾을 수 없습니다. 'smiles' 또는 'SMILES' 컬럼이 필요합니다.")

        # Process and insert
        conn = get_db()
        inserted = 0
        skipped = 0

        for _, row in df.iterrows():
            smiles = str(row.get("smiles", "")).strip()
            if not smiles or smiles == "nan":
                skipped += 1
                continue

            # Compute features if RDKit available
            features = {}
            try:
                from rdkit import Chem
                from rdkit.Chem import Descriptors, AllChem
                mol = Chem.MolFromSmiles(smiles)
                if mol:
                    smiles = Chem.MolToSmiles(mol)  # canonical
                    features = {
                        "molecular_weight": round(Descriptors.MolWt(mol), 2),
                        "logp": round(Descriptors.MolLogP(mol), 3),
                        "hbd": Descriptors.NumHDonors(mol),
                        "hba": Descriptors.NumHAcceptors(mol),
                        "tpsa": round(Descriptors.TPSA(mol), 2),
                        "rotatable_bonds": Descriptors.NumRotatableBonds(mol),
                        "aromatic_rings": Descriptors.NumAromaticRings(mol),
                    }
                    # Embedding
                    fp = AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=2048)
                    fp_arr = np.array(fp, dtype=np.float32)
                    rng = np.random.RandomState(42)
                    proj = rng.randn(2048, 768).astype(np.float32)
                    proj /= np.linalg.norm(proj, axis=0, keepdims=True)
                    emb = fp_arr @ proj
                    norm = np.linalg.norm(emb)
                    if norm > 0: emb = emb / norm
                    features["embedding"] = json.dumps(emb.tolist())
                else:
                    skipped += 1
                    continue
            except ImportError:
                features["embedding"] = "[]"

            def get_val(col):
                v = row.get(col)
                if v is not None and str(v) != "nan":
                    try: return float(v)
                    except: return None
                return None

            conn.execute("""
                INSERT INTO materials (name, smiles, category, molecular_weight, logp, hbd, hba, tpsa,
                    rotatable_bonds, aromatic_rings, thermal_stability, dielectric_constant,
                    bandgap, solubility, density, embedding, source, is_verified)
                VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, 0)
            """, (
                str(row.get("name", "")) if str(row.get("name", "")) != "nan" else None,
                smiles,
                category or (str(row.get("category", "")) if str(row.get("category", "")) != "nan" else None),
                features.get("molecular_weight") or get_val("molecular_weight"),
                features.get("logp"), features.get("hbd"), features.get("hba"), features.get("tpsa"),
                features.get("rotatable_bonds"), features.get("aromatic_rings"),
                get_val("thermal_stability"),
                get_val("dielectric_constant"),
                get_val("bandgap"),
                get_val("solubility"),
                get_val("density"),
                features.get("embedding", "[]"),
                source,
            ))
            inserted += 1

        conn.commit()
        total = conn.execute("SELECT COUNT(*) FROM materials").fetchone()[0]
        conn.close()

        return {
            "success": True,
            "data": {
                "filename": file.filename,
                "rows_total": len(df),
                "inserted": inserted,
                "skipped": skipped,
                "total_in_db": total,
            },
        }
    except HTTPException:
        raise
    except Exception as e:
        raise HTTPException(500, f"파일 처리 실패: {e}")


@app.get("/api/data/summary")
def data_summary():
    """연구 데이터 현황 요약."""
    conn = get_db()
    total = conn.execute("SELECT COUNT(*) FROM materials").fetchone()[0]
    by_source = conn.execute(
        "SELECT source, COUNT(*) as cnt FROM materials GROUP BY source ORDER BY cnt DESC"
    ).fetchall()
    by_category = conn.execute(
        "SELECT category, COUNT(*) as cnt FROM materials GROUP BY category ORDER BY cnt DESC"
    ).fetchall()
    verified = conn.execute("SELECT COUNT(*) FROM materials WHERE is_verified = 1").fetchone()[0]

    # Property coverage
    coverage = {}
    for prop in ["thermal_stability", "dielectric_constant", "bandgap", "solubility", "density"]:
        cnt = conn.execute(f"SELECT COUNT(*) FROM materials WHERE {prop} IS NOT NULL").fetchone()[0]
        coverage[prop] = {"count": cnt, "ratio": round(cnt / total * 100, 1) if total > 0 else 0}

    conn.close()
    return {
        "success": True,
        "data": {
            "total_materials": total,
            "verified": verified,
            "by_source": {r["source"] or "unknown": r["cnt"] for r in by_source},
            "by_category": {r["category"] or "미분류": r["cnt"] for r in by_category},
            "property_coverage": coverage,
        },
    }


# --- AI Discovery (New Material Proposal) ---

class DiscoveryRequest(BaseModel):
    target_category: str = ""
    target_properties: dict = Field(default_factory=dict)
    constraints: str = ""
    data_limit: int = Field(default=30, ge=5, le=100)


@app.post("/api/ai/discover")
async def ai_discover_materials(request: DiscoveryRequest):
    """
    3개 AI 모델이 기존 연구 데이터를 학습하여 새로운 신소재를 제안합니다.
    """
    from app.services.ai_discovery import discover_new_materials

    conn = get_db()
    if request.target_category and request.target_category != "all":
        rows = conn.execute(
            "SELECT * FROM materials WHERE LOWER(category) = LOWER(?) ORDER BY thermal_stability DESC LIMIT ?",
            (request.target_category, request.data_limit),
        ).fetchall()
    else:
        rows = conn.execute(
            "SELECT * FROM materials ORDER BY is_verified DESC, thermal_stability DESC LIMIT ?",
            (request.data_limit,),
        ).fetchall()
    conn.close()

    materials = [dict(r) for r in rows]
    for m in materials:
        m.pop("embedding", None)

    result = await discover_new_materials(
        existing_materials=materials,
        target_category=request.target_category,
        target_properties=request.target_properties,
        constraints=request.constraints,
    )
    return {"success": True, "data": result}


# --- Simulation ---

class SimulationRequest(BaseModel):
    smiles: str
    name: str = ""
    category: str = ""


@app.post("/api/ai/simulate")
async def ai_simulate_material(request: SimulationRequest):
    """
    제안된 소재에 대해 3개 AI 모델이 시뮬레이션 예측을 수행합니다.
    """
    from app.services.ai_discovery import simulate_material

    # Get some existing materials for comparison
    conn = get_db()
    if request.category:
        rows = conn.execute(
            "SELECT * FROM materials WHERE LOWER(category) = LOWER(?) LIMIT 10",
            (request.category,),
        ).fetchall()
    else:
        rows = conn.execute("SELECT * FROM materials WHERE is_verified = 1 LIMIT 10").fetchall()
    conn.close()

    existing = [dict(r) for r in rows]
    for m in existing:
        m.pop("embedding", None)

    result = await simulate_material(
        smiles=request.smiles,
        name=request.name,
        category=request.category,
        existing_materials=existing,
    )
    return {"success": True, "data": result}


# --- Pipeline (simulated) ---

@app.get("/api/pipeline")
def list_pipelines():
    return {"success": True, "data": []}


@app.post("/api/pipeline/start")
def start_pipeline(request: dict):
    return {
        "success": True,
        "data": {
            "pipeline_id": "demo-001",
            "status": "running",
            "steps": [],
        },
    }


# ─── Run ─────────────────────────────────────────

if __name__ == "__main__":
    import uvicorn

    if not os.path.exists(DB_PATH):
        print(f"ERROR: Database not found at {DB_PATH}")
        print("Run first: python -m scripts.seed_materials")
        exit(1)

    print(f"\n  Materials AI Platform API (PoC)")
    print(f"  Database: {DB_PATH}")
    conn = get_db()
    count = conn.execute("SELECT COUNT(*) FROM materials").fetchone()[0]
    conn.close()
    print(f"  Materials in DB: {count}")
    print(f"  Docs: http://localhost:8000/docs\n")

    uvicorn.run(app, host="0.0.0.0", port=8000)
