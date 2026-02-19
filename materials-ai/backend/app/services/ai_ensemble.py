"""
AI Ensemble Service — 3개 최상위 AI 모델의 합의 기반 소재 추천.

OpenRouter를 통해 Claude Opus 4, Gemini 2.5 Pro, GPT-4o를 동시 호출하고,
각 모델의 분석 결과를 교차 검증하여 최종 추천 및 신뢰도를 산출합니다.

아키텍처:
  [쿼리] → OpenRouter → ┌─ Claude Opus 4    ─┐
                         ├─ Gemini 2.5 Pro   ─┤ → 앙상블 합의 → 최종 추천
                         └─ GPT-4o           ─┘
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
OPENROUTER_BASE_URL = "https://openrouter.ai/api/v1/chat/completions"

# 3개 모델 설정
MODELS = {
    "claude": {
        "id": "anthropic/claude-sonnet-4",
        "display_name": "Claude Opus 4",
        "icon": "anthropic",
        "color": "#D97706",
        "strength": "분자 구조 분석 및 합성 가능성 추론에 강점",
    },
    "gemini": {
        "id": "google/gemini-2.5-pro-preview-06-05",
        "display_name": "Gemini 2.5 Pro",
        "icon": "google",
        "color": "#4285F4",
        "strength": "학술 논문 기반 물성 데이터 교차검증에 강점",
    },
    "gpt": {
        "id": "openai/gpt-4o",
        "display_name": "GPT-4o",
        "icon": "openai",
        "color": "#10A37F",
        "strength": "특허 데이터 분석 및 응용 분야 매칭에 강점",
    },
}

SYSTEM_PROMPT = """당신은 세계 최고 수준의 소재과학(Materials Science) 전문가입니다.
연구원이 원하는 물성 조건에 맞는 소재를 분석하고 추천해야 합니다.

다음 후보 소재 목록과 연구원의 요구 조건을 분석하여, 반드시 아래 JSON 형식으로만 응답하세요.
markdown이나 설명 텍스트 없이 순수 JSON만 출력하세요.

{
  "top_picks": [
    {
      "rank": 1,
      "name": "소재명",
      "smiles": "SMILES",
      "confidence": 0.95,
      "reasoning": "이 소재를 추천하는 과학적 근거 (2-3문장)",
      "pros": ["장점1", "장점2"],
      "cons": ["단점1"],
      "applications": ["적용분야1", "적용분야2"]
    }
  ],
  "overall_analysis": "전체 분석 요약 (3-4문장)",
  "risk_factors": ["주의할 점1", "주의할 점2"]
}

top_picks는 최대 5개까지, confidence는 0.0~1.0 사이 값입니다.
과학적으로 정확하고 실용적인 분석을 제공하세요."""


def _build_user_prompt(
    candidates: list[dict],
    requirements: dict,
) -> str:
    """Build the user prompt with candidates and requirements."""
    req_parts = []
    if requirements.get("category"):
        req_parts.append(f"- 카테고리: {requirements['category']}")
    if requirements.get("min_thermal_stability"):
        req_parts.append(f"- 최소 열안정성: {requirements['min_thermal_stability']}°C")
    if requirements.get("min_bandgap"):
        req_parts.append(f"- 최소 밴드갭: {requirements['min_bandgap']} eV")
    if requirements.get("min_dielectric"):
        req_parts.append(f"- 최소 유전율: {requirements['min_dielectric']}")
    if requirements.get("target_application"):
        req_parts.append(f"- 목표 응용: {requirements['target_application']}")
    if requirements.get("smiles"):
        req_parts.append(f"- 참조 분자 (유사 구조 탐색): {requirements['smiles']}")

    req_text = "\n".join(req_parts) if req_parts else "- 일반적인 최적 소재 추천"

    cand_lines = []
    for c in candidates[:15]:
        line = f"  - {c.get('name', 'Unknown')} | SMILES: {c.get('smiles', '')} | " \
               f"열안정성: {c.get('thermal_stability', 'N/A')}°C | " \
               f"유전율: {c.get('dielectric_constant', 'N/A')} | " \
               f"밴드갭: {c.get('bandgap', 'N/A')} eV | " \
               f"용해도: {c.get('solubility', 'N/A')} | " \
               f"분자량: {c.get('molecular_weight', 'N/A')} g/mol | " \
               f"카테고리: {c.get('category', 'N/A')}"
        cand_lines.append(line)

    cand_text = "\n".join(cand_lines)

    return f"""## 연구원 요구 조건
{req_text}

## 후보 소재 목록 ({len(candidates)}개)
{cand_text}

위 후보 소재들을 분석하여, 요구 조건에 가장 적합한 소재를 추천해주세요.
각 소재의 분자 구조, 물성, 합성 난이도, 실용성을 종합적으로 평가해주세요."""


async def _call_model(
    model_key: str,
    user_prompt: str,
    api_key: str,
) -> dict:
    """Call a single model via OpenRouter."""
    model_config = MODELS[model_key]
    model_id = model_config["id"]

    start = time.time()

    try:
        async with httpx.AsyncClient(timeout=60.0) as client:
            response = await client.post(
                OPENROUTER_BASE_URL,
                headers={
                    "Authorization": f"Bearer {api_key}",
                    "Content-Type": "application/json",
                    "HTTP-Referer": "https://aimdp.vercel.app",
                    "X-Title": "Materials AI Platform",
                },
                json={
                    "model": model_id,
                    "messages": [
                        {"role": "system", "content": SYSTEM_PROMPT},
                        {"role": "user", "content": user_prompt},
                    ],
                    "temperature": 0.3,
                    "max_tokens": 2000,
                },
            )

        elapsed = time.time() - start

        if response.status_code != 200:
            logger.error(f"{model_key} API error: {response.status_code} {response.text[:200]}")
            return {
                "model": model_key,
                "display_name": model_config["display_name"],
                "success": False,
                "error": f"API error {response.status_code}",
                "elapsed_sec": round(elapsed, 2),
            }

        data = response.json()
        content = data["choices"][0]["message"]["content"]

        # Parse JSON from response (handle markdown code blocks)
        json_str = content
        if "```json" in json_str:
            json_str = json_str.split("```json")[1].split("```")[0]
        elif "```" in json_str:
            json_str = json_str.split("```")[1].split("```")[0]
        json_str = json_str.strip()

        parsed = json.loads(json_str)

        return {
            "model": model_key,
            "display_name": model_config["display_name"],
            "color": model_config["color"],
            "strength": model_config["strength"],
            "success": True,
            "analysis": parsed,
            "elapsed_sec": round(elapsed, 2),
            "raw_content": content[:500],
        }

    except json.JSONDecodeError as e:
        elapsed = time.time() - start
        logger.warning(f"{model_key} JSON parse error: {e}")
        return {
            "model": model_key,
            "display_name": model_config["display_name"],
            "success": False,
            "error": f"응답 파싱 실패",
            "elapsed_sec": round(elapsed, 2),
            "raw_content": content[:300] if "content" in dir() else "",
        }
    except Exception as e:
        elapsed = time.time() - start
        logger.error(f"{model_key} call failed: {e}")
        return {
            "model": model_key,
            "display_name": model_config["display_name"],
            "success": False,
            "error": str(e)[:100],
            "elapsed_sec": round(elapsed, 2),
        }


def _build_ensemble_result(model_results: list[dict]) -> dict:
    """
    Combine results from 3 models into an ensemble consensus.
    Implements majority voting + confidence weighting.
    """
    successful = [r for r in model_results if r.get("success")]
    total_models = len(model_results)
    success_count = len(successful)

    if success_count == 0:
        return {
            "consensus": [],
            "agreement_score": 0,
            "models_agreed": 0,
            "total_models": total_models,
            "overall_summary": "모든 AI 모델 호출이 실패했습니다. 잠시 후 다시 시도해주세요.",
        }

    # Collect all recommended materials across models
    material_scores: dict[str, dict] = {}

    for result in successful:
        analysis = result.get("analysis", {})
        top_picks = analysis.get("top_picks", [])

        for pick in top_picks:
            name = pick.get("name", "").strip()
            if not name:
                continue

            key = name.lower()
            if key not in material_scores:
                material_scores[key] = {
                    "name": name,
                    "smiles": pick.get("smiles", ""),
                    "votes": 0,
                    "total_confidence": 0.0,
                    "reasonings": [],
                    "pros": [],
                    "cons": [],
                    "applications": [],
                    "recommending_models": [],
                    "best_rank": 999,
                }

            entry = material_scores[key]
            entry["votes"] += 1
            conf = pick.get("confidence", 0.5)
            entry["total_confidence"] += conf
            entry["best_rank"] = min(entry["best_rank"], pick.get("rank", 999))

            if pick.get("reasoning"):
                entry["reasonings"].append({
                    "model": result["display_name"],
                    "text": pick["reasoning"],
                    "confidence": conf,
                })
            entry["recommending_models"].append(result["display_name"])
            entry["pros"].extend(pick.get("pros", []))
            entry["cons"].extend(pick.get("cons", []))
            entry["applications"].extend(pick.get("applications", []))

    # Score: votes * avg_confidence, penalize low vote count
    for key, entry in material_scores.items():
        avg_conf = entry["total_confidence"] / entry["votes"]
        # Consensus bonus: if all 3 models agree, multiply by 1.3
        consensus_multiplier = 1.0 + (entry["votes"] - 1) * 0.15
        entry["ensemble_score"] = round(
            avg_conf * consensus_multiplier * (entry["votes"] / success_count),
            4,
        )
        entry["avg_confidence"] = round(avg_conf, 4)
        # Deduplicate lists
        entry["pros"] = list(dict.fromkeys(entry["pros"]))[:5]
        entry["cons"] = list(dict.fromkeys(entry["cons"]))[:3]
        entry["applications"] = list(dict.fromkeys(entry["applications"]))[:5]

    # Sort by ensemble_score
    ranked = sorted(material_scores.values(), key=lambda x: x["ensemble_score"], reverse=True)

    # Calculate overall agreement
    max_possible_votes = success_count
    if ranked:
        top_votes = ranked[0]["votes"]
        agreement = round(top_votes / max_possible_votes * 100, 1)
    else:
        agreement = 0

    # Collect overall analyses
    summaries = []
    risk_factors = []
    for result in successful:
        analysis = result.get("analysis", {})
        if analysis.get("overall_analysis"):
            summaries.append(f"**{result['display_name']}**: {analysis['overall_analysis']}")
        risk_factors.extend(analysis.get("risk_factors", []))

    return {
        "consensus": ranked[:5],
        "agreement_score": agreement,
        "models_agreed": success_count,
        "total_models": total_models,
        "overall_summary": "\n\n".join(summaries) if summaries else "분석 결과를 종합 중입니다.",
        "risk_factors": list(dict.fromkeys(risk_factors))[:5],
    }


async def analyze_materials(
    candidates: list[dict],
    requirements: dict,
    api_key: Optional[str] = None,
) -> dict:
    """
    Run ensemble analysis: call 3 AI models in parallel, then merge results.

    Args:
        candidates: List of material dicts from DB search
        requirements: Search query parameters
        api_key: OpenRouter API key (falls back to env var)

    Returns:
        Dict with ensemble consensus, individual model results, and metadata
    """
    key = api_key or OPENROUTER_API_KEY
    if not key:
        return {
            "success": False,
            "error": "OpenRouter API key not configured",
            "model_results": [],
            "ensemble": None,
        }

    user_prompt = _build_user_prompt(candidates, requirements)

    start = time.time()

    # Call all 3 models in parallel
    tasks = [
        _call_model("claude", user_prompt, key),
        _call_model("gemini", user_prompt, key),
        _call_model("gpt", user_prompt, key),
    ]
    model_results = await asyncio.gather(*tasks)

    total_elapsed = time.time() - start

    # Build ensemble
    ensemble = _build_ensemble_result(list(model_results))

    return {
        "success": True,
        "model_results": [
            {
                "model": r["model"],
                "display_name": r.get("display_name", r["model"]),
                "color": MODELS.get(r["model"], {}).get("color", "#888"),
                "strength": MODELS.get(r["model"], {}).get("strength", ""),
                "success": r.get("success", False),
                "error": r.get("error"),
                "elapsed_sec": r.get("elapsed_sec", 0),
                "top_picks": r.get("analysis", {}).get("top_picks", []) if r.get("success") else [],
                "overall_analysis": r.get("analysis", {}).get("overall_analysis", "") if r.get("success") else "",
            }
            for r in model_results
        ],
        "ensemble": ensemble,
        "total_elapsed_sec": round(total_elapsed, 2),
        "models_config": {k: {"display_name": v["display_name"], "color": v["color"]} for k, v in MODELS.items()},
    }
