"""
AI Discovery & Simulation Engine.

1. 기존 연구 데이터를 기반으로 새로운 신소재 분자 구조 제안
2. 제안된 소재의 물성 예측 + 합성 가능성 + 비용 + 위험도 시뮬레이션
3. 3개 AI 모델(Claude, Gemini, GPT)이 각각 독립 분석 후 합의

OpenRouter를 통해 3개 모델에 동진쎄미켐의 연구 데이터 컨텍스트를 전달하고,
새로운 분자 구조를 역설계합니다.
"""

import asyncio
import json
import logging
import os
import time
from typing import Optional

import httpx

logger = logging.getLogger(__name__)

OPENROUTER_API_KEY = os.environ.get("OPENROUTER_API_KEY", "")
OPENROUTER_URL = "https://openrouter.ai/api/v1/chat/completions"

MODELS = {
    "claude": "anthropic/claude-sonnet-4",
    "gemini": "google/gemini-2.5-pro-preview-06-05",
    "gpt": "openai/gpt-4o",
}

MODEL_DISPLAY = {
    "claude": {"name": "Claude Opus 4", "color": "#D97706"},
    "gemini": {"name": "Gemini 2.5 Pro", "color": "#4285F4"},
    "gpt": {"name": "GPT-4o", "color": "#10A37F"},
}

DISCOVERY_SYSTEM_PROMPT = """당신은 세계 최고 수준의 소재과학자이자 계산화학 전문가입니다.
동진쎄미켐 R&D 연구소의 기존 연구 데이터를 분석하고, 새로운 신소재 분자 구조를 제안해야 합니다.

기존 데이터에서 구조-물성 관계(Structure-Property Relationship)를 파악하고,
아직 탐색되지 않은 화학 공간에서 목표 물성을 달성할 수 있는 새로운 분자를 설계하세요.

반드시 아래 JSON 형식으로만 응답하세요. markdown 없이 순수 JSON만 출력하세요.

{
  "proposed_materials": [
    {
      "name": "제안 소재명",
      "smiles": "유효한 SMILES 문자열",
      "design_rationale": "이 분자를 설계한 과학적 근거 (3-4문장)",
      "predicted_properties": {
        "thermal_stability": 숫자(°C),
        "dielectric_constant": 숫자,
        "bandgap": 숫자(eV),
        "solubility": 숫자(0-1),
        "density": 숫자(g/cm³)
      },
      "novelty_score": 0.0~1.0,
      "synthesizability": "easy|moderate|difficult|very_difficult",
      "estimated_cost_level": "low|medium|high|very_high",
      "key_building_blocks": ["핵심 출발물질1", "출발물질2"],
      "target_applications": ["응용분야1", "응용분야2"]
    }
  ],
  "data_insights": "기존 데이터에서 발견한 핵심 패턴 (3-4문장)",
  "research_direction": "향후 연구 방향 제안 (2-3문장)"
}

proposed_materials는 3~5개를 제안하세요.
각 분자의 SMILES는 반드시 유효한 화학 구조여야 합니다.
기존 데이터에 없는 새로운 구조를 제안하되, 합성 가능성을 반드시 고려하세요."""

SIMULATION_SYSTEM_PROMPT = """당신은 계산화학 및 공정 시뮬레이션 전문가입니다.
주어진 분자 구조에 대해 상세한 시뮬레이션 결과를 예측해야 합니다.

반드시 아래 JSON 형식으로만 응답하세요. markdown 없이 순수 JSON만 출력하세요.

{
  "simulation_results": {
    "molecular_dynamics": {
      "glass_transition_temp": 숫자(°C),
      "decomposition_temp": 숫자(°C),
      "thermal_expansion_coeff": 숫자(ppm/°C),
      "stability_rating": "excellent|good|moderate|poor"
    },
    "electronic_properties": {
      "homo_energy": 숫자(eV),
      "lumo_energy": 숫자(eV),
      "bandgap": 숫자(eV),
      "electron_mobility": 숫자(cm²/V·s),
      "hole_mobility": 숫자(cm²/V·s)
    },
    "synthesis_simulation": {
      "recommended_route": "합성 경로 설명 (2-3문장)",
      "key_reactions": ["반응1", "반응2", "반응3"],
      "estimated_yield": 숫자(%),
      "estimated_steps": 숫자,
      "difficulty": "easy|moderate|difficult|very_difficult",
      "estimated_cost_per_gram": "숫자(USD)",
      "scale_up_feasibility": "excellent|good|moderate|poor"
    },
    "environmental_assessment": {
      "toxicity_risk": "low|medium|high",
      "biodegradability": "good|moderate|poor",
      "regulatory_concerns": "우려사항 또는 '없음'"
    },
    "performance_prediction": {
      "oled_efficiency": "해당 시 예측값, 비해당 시 null",
      "battery_capacity": "해당 시 예측값, 비해당 시 null",
      "coating_hardness": "해당 시 예측값, 비해당 시 null",
      "overall_score": 0.0~1.0,
      "confidence": 0.0~1.0
    }
  },
  "comparison_with_existing": "기존 소재 대비 장단점 비교 (3-4문장)",
  "recommendation": "최종 추천 의견 (2-3문장)"
}"""


async def _call_openrouter(model_key: str, system: str, user: str) -> dict:
    """Call a model via OpenRouter."""
    api_key = OPENROUTER_API_KEY
    if not api_key:
        return {"success": False, "error": "API key not set", "model": model_key}

    start = time.time()
    try:
        async with httpx.AsyncClient(timeout=60.0) as client:
            resp = await client.post(
                OPENROUTER_URL,
                headers={
                    "Authorization": f"Bearer {api_key}",
                    "Content-Type": "application/json",
                    "HTTP-Referer": "https://aimdp.vercel.app",
                    "X-Title": "Materials AI Platform",
                },
                json={
                    "model": MODELS[model_key],
                    "messages": [
                        {"role": "system", "content": system},
                        {"role": "user", "content": user},
                    ],
                    "temperature": 0.4,
                    "max_tokens": 3000,
                },
            )

        elapsed = time.time() - start
        if resp.status_code != 200:
            return {"success": False, "error": f"HTTP {resp.status_code}", "model": model_key, "elapsed": elapsed}

        content = resp.json()["choices"][0]["message"]["content"]
        json_str = content
        if "```json" in json_str:
            json_str = json_str.split("```json")[1].split("```")[0]
        elif "```" in json_str:
            json_str = json_str.split("```")[1].split("```")[0]

        parsed = json.loads(json_str.strip())
        return {"success": True, "model": model_key, "data": parsed, "elapsed": round(elapsed, 2)}

    except json.JSONDecodeError:
        return {"success": False, "error": "JSON parse failed", "model": model_key, "elapsed": round(time.time() - start, 2)}
    except Exception as e:
        return {"success": False, "error": str(e)[:100], "model": model_key, "elapsed": round(time.time() - start, 2)}


def _format_research_data(materials: list[dict]) -> str:
    """Format existing research data for AI context."""
    lines = []
    for m in materials:
        line = (
            f"  {m.get('name', 'N/A')} | SMILES: {m.get('smiles', '')} | "
            f"카테고리: {m.get('category', 'N/A')} | "
            f"열안정성: {m.get('thermal_stability', 'N/A')}°C | "
            f"유전율: {m.get('dielectric_constant', 'N/A')} | "
            f"밴드갭: {m.get('bandgap', 'N/A')} eV | "
            f"용해도: {m.get('solubility', 'N/A')} | "
            f"밀도: {m.get('density', 'N/A')} g/cm³ | "
            f"분자량: {m.get('molecular_weight', 'N/A')} | "
            f"검증: {'O' if m.get('is_verified') else 'X'}"
        )
        lines.append(line)
    return "\n".join(lines)


async def discover_new_materials(
    existing_materials: list[dict],
    target_category: str = "",
    target_properties: dict = None,
    constraints: str = "",
) -> dict:
    """
    3개 AI 모델이 기존 연구 데이터를 분석하여 새로운 신소재를 제안합니다.
    """
    data_context = _format_research_data(existing_materials)

    target_parts = []
    if target_category:
        target_parts.append(f"- 목표 카테고리: {target_category}")
    if target_properties:
        for k, v in target_properties.items():
            target_parts.append(f"- 목표 {k}: {v}")
    if constraints:
        target_parts.append(f"- 추가 조건: {constraints}")

    target_text = "\n".join(target_parts) if target_parts else "- 기존 소재보다 성능이 우수한 새로운 분자 구조"

    user_prompt = f"""## 동진쎄미켐 기존 연구 데이터 ({len(existing_materials)}건)
{data_context}

## 목표 조건
{target_text}

위 연구 데이터를 분석하여:
1. 구조-물성 관계(SPR) 패턴을 파악하세요
2. 기존 데이터에 없는 새로운 분자 구조를 3~5개 제안하세요
3. 각 제안 분자의 예상 물성, 합성 가능성, 비용을 예측하세요
4. 기존 소재 대비 어떤 점이 개선되는지 설명하세요"""

    # Call 3 models in parallel
    tasks = [_call_openrouter(k, DISCOVERY_SYSTEM_PROMPT, user_prompt) for k in MODELS]
    results = await asyncio.gather(*tasks)

    # Merge proposals from all models
    all_proposals = []
    model_summaries = []
    for r in results:
        display = MODEL_DISPLAY.get(r["model"], {})
        entry = {
            "model": r["model"],
            "display_name": display.get("name", r["model"]),
            "color": display.get("color", "#888"),
            "success": r["success"],
            "error": r.get("error"),
            "elapsed": r.get("elapsed", 0),
        }
        if r["success"]:
            data = r["data"]
            proposals = data.get("proposed_materials", [])
            for p in proposals:
                p["proposed_by"] = display.get("name", r["model"])
                p["model_color"] = display.get("color", "#888")
            all_proposals.extend(proposals)
            entry["data_insights"] = data.get("data_insights", "")
            entry["research_direction"] = data.get("research_direction", "")
            entry["proposal_count"] = len(proposals)
            model_summaries.append(entry)
        else:
            entry["proposal_count"] = 0
            model_summaries.append(entry)

    # Deduplicate by SMILES and merge votes
    unique = {}
    for p in all_proposals:
        smiles = p.get("smiles", "").strip()
        if not smiles:
            continue
        key = smiles.lower()
        if key not in unique:
            unique[key] = {**p, "votes": 1, "proposed_by_list": [p.get("proposed_by", "")]}
        else:
            unique[key]["votes"] += 1
            unique[key]["proposed_by_list"].append(p.get("proposed_by", ""))
            # Keep higher novelty score
            if p.get("novelty_score", 0) > unique[key].get("novelty_score", 0):
                unique[key]["novelty_score"] = p["novelty_score"]

    merged = sorted(unique.values(), key=lambda x: (-x["votes"], -x.get("novelty_score", 0)))

    return {
        "success": True,
        "proposals": merged[:8],
        "total_proposals": len(all_proposals),
        "unique_proposals": len(merged),
        "model_summaries": model_summaries,
        "data_used": len(existing_materials),
    }


async def simulate_material(
    smiles: str,
    name: str = "",
    category: str = "",
    existing_materials: list[dict] = None,
) -> dict:
    """
    제안된 소재에 대해 3개 AI 모델이 시뮬레이션 예측을 수행합니다.
    """
    context = ""
    if existing_materials:
        context = f"\n\n## 비교용 기존 소재 데이터\n{_format_research_data(existing_materials[:10])}"

    user_prompt = f"""## 시뮬레이션 대상 분자
- 이름: {name or '신규 제안 소재'}
- SMILES: {smiles}
- 카테고리: {category or '미정'}
{context}

위 분자에 대해 다음을 시뮬레이션하세요:
1. 분자동역학(MD) 시뮬레이션: 열안정성, 유리전이온도, 분해온도
2. 전자구조 계산: HOMO, LUMO, 밴드갭, 전하 이동도
3. 합성 경로 시뮬레이션: 추천 합성 경로, 예상 수율, 단계 수, 비용
4. 환경 평가: 독성, 생분해성, 규제 우려
5. 성능 예측: 해당 카테고리에서의 예상 성능"""

    tasks = [_call_openrouter(k, SIMULATION_SYSTEM_PROMPT, user_prompt) for k in MODELS]
    results = await asyncio.gather(*tasks)

    # Merge simulation results (average numeric values)
    sim_results = []
    model_details = []
    for r in results:
        display = MODEL_DISPLAY.get(r["model"], {})
        entry = {
            "model": r["model"],
            "display_name": display.get("name", r["model"]),
            "color": display.get("color", "#888"),
            "success": r["success"],
            "error": r.get("error"),
            "elapsed": r.get("elapsed", 0),
        }
        if r["success"]:
            entry["simulation"] = r["data"].get("simulation_results", {})
            entry["comparison"] = r["data"].get("comparison_with_existing", "")
            entry["recommendation"] = r["data"].get("recommendation", "")
            sim_results.append(r["data"].get("simulation_results", {}))
        model_details.append(entry)

    # Average numeric predictions across models
    averaged = _average_simulations(sim_results)

    return {
        "success": True,
        "smiles": smiles,
        "name": name,
        "category": category,
        "averaged_prediction": averaged,
        "model_details": model_details,
        "models_succeeded": sum(1 for r in results if r["success"]),
        "total_models": len(results),
    }


def _average_simulations(sims: list[dict]) -> dict:
    """Average numeric values across multiple simulation results."""
    if not sims:
        return {}

    def avg_field(path: list[str]) -> Optional[float]:
        values = []
        for s in sims:
            obj = s
            for key in path:
                if isinstance(obj, dict):
                    obj = obj.get(key)
                else:
                    obj = None
                    break
            if obj is not None and isinstance(obj, (int, float)):
                values.append(float(obj))
        return round(sum(values) / len(values), 3) if values else None

    def mode_field(path: list[str]) -> Optional[str]:
        values = []
        for s in sims:
            obj = s
            for key in path:
                if isinstance(obj, dict):
                    obj = obj.get(key)
                else:
                    obj = None
                    break
            if obj and isinstance(obj, str):
                values.append(obj)
        if not values:
            return None
        return max(set(values), key=values.count)

    return {
        "molecular_dynamics": {
            "glass_transition_temp": avg_field(["molecular_dynamics", "glass_transition_temp"]),
            "decomposition_temp": avg_field(["molecular_dynamics", "decomposition_temp"]),
            "stability_rating": mode_field(["molecular_dynamics", "stability_rating"]),
        },
        "electronic_properties": {
            "homo_energy": avg_field(["electronic_properties", "homo_energy"]),
            "lumo_energy": avg_field(["electronic_properties", "lumo_energy"]),
            "bandgap": avg_field(["electronic_properties", "bandgap"]),
            "electron_mobility": avg_field(["electronic_properties", "electron_mobility"]),
        },
        "synthesis": {
            "estimated_yield": avg_field(["synthesis_simulation", "estimated_yield"]),
            "estimated_steps": avg_field(["synthesis_simulation", "estimated_steps"]),
            "difficulty": mode_field(["synthesis_simulation", "difficulty"]),
            "scale_up_feasibility": mode_field(["synthesis_simulation", "scale_up_feasibility"]),
        },
        "environment": {
            "toxicity_risk": mode_field(["environmental_assessment", "toxicity_risk"]),
            "biodegradability": mode_field(["environmental_assessment", "biodegradability"]),
        },
        "performance": {
            "overall_score": avg_field(["performance_prediction", "overall_score"]),
            "confidence": avg_field(["performance_prediction", "confidence"]),
        },
    }
