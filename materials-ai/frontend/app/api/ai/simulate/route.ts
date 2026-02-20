import { NextRequest, NextResponse } from "next/server";

const OPENROUTER_URL = "https://openrouter.ai/api/v1/chat/completions";
const API_KEY = process.env.OPENROUTER_API_KEY || "";

const MODELS = [
  { key: "claude", id: "anthropic/claude-sonnet-4", name: "Claude Opus 4", color: "#D97706" },
  { key: "gemini", id: "google/gemini-2.5-pro-preview-06-05", name: "Gemini 2.5 Pro", color: "#4285F4" },
  { key: "gpt", id: "openai/gpt-4o", name: "GPT-4o", color: "#10A37F" },
];

const SYSTEM_PROMPT = `당신은 계산화학 및 공정 시뮬레이션 전문가입니다.
주어진 분자 구조에 대해 상세한 시뮬레이션 결과를 예측해야 합니다.

반드시 아래 JSON 형식으로만 응답하세요. markdown 없이 순수 JSON만 출력하세요.

{
  "simulation_results": {
    "molecular_dynamics": {
      "glass_transition_temp": 숫자,
      "decomposition_temp": 숫자,
      "stability_rating": "excellent|good|moderate|poor"
    },
    "electronic_properties": {
      "homo_energy": 숫자,
      "lumo_energy": 숫자,
      "bandgap": 숫자,
      "electron_mobility": 숫자
    },
    "synthesis_simulation": {
      "estimated_yield": 숫자,
      "estimated_steps": 숫자,
      "difficulty": "easy|moderate|difficult|very_difficult",
      "scale_up_feasibility": "excellent|good|moderate|poor"
    },
    "environmental_assessment": {
      "toxicity_risk": "low|medium|high",
      "biodegradability": "good|moderate|poor"
    },
    "performance_prediction": {
      "overall_score": 0.0~1.0,
      "confidence": 0.0~1.0
    }
  },
  "comparison_with_existing": "기존 소재 대비 비교 (2-3문장)",
  "recommendation": "최종 추천 의견 (2-3문장)"
}`;

async function callModel(modelId: string, userPrompt: string): Promise<Record<string, unknown>> {
  const start = Date.now();
  try {
    const res = await fetch(OPENROUTER_URL, {
      method: "POST",
      headers: {
        Authorization: `Bearer ${API_KEY}`,
        "Content-Type": "application/json",
        "HTTP-Referer": "https://aimdp.vercel.app",
        "X-Title": "Materials AI Platform",
      },
      body: JSON.stringify({
        model: modelId,
        messages: [
          { role: "system", content: SYSTEM_PROMPT },
          { role: "user", content: userPrompt },
        ],
        temperature: 0.3,
        max_tokens: 2500,
      }),
    });

    const elapsed = (Date.now() - start) / 1000;
    if (!res.ok) return { success: false, error: `HTTP ${res.status}`, elapsed };

    const data = await res.json();
    let content = data.choices?.[0]?.message?.content || "";
    if (content.includes("```json")) content = content.split("```json")[1].split("```")[0];
    else if (content.includes("```")) content = content.split("```")[1].split("```")[0];

    return { success: true, data: JSON.parse(content.trim()), elapsed: Math.round(elapsed * 100) / 100 };
  } catch (e) {
    return { success: false, error: String(e).slice(0, 100), elapsed: (Date.now() - start) / 1000 };
  }
}

function avgField(sims: Record<string, unknown>[], path: string[]): number | null {
  const values: number[] = [];
  for (const s of sims) {
    let obj: unknown = s;
    for (const k of path) {
      if (obj && typeof obj === "object") obj = (obj as Record<string, unknown>)[k];
      else { obj = null; break; }
    }
    if (typeof obj === "number") values.push(obj);
  }
  return values.length > 0 ? Math.round((values.reduce((a, b) => a + b, 0) / values.length) * 1000) / 1000 : null;
}

function modeField(sims: Record<string, unknown>[], path: string[]): string | null {
  const values: string[] = [];
  for (const s of sims) {
    let obj: unknown = s;
    for (const k of path) {
      if (obj && typeof obj === "object") obj = (obj as Record<string, unknown>)[k];
      else { obj = null; break; }
    }
    if (typeof obj === "string") values.push(obj);
  }
  if (!values.length) return null;
  const counts: Record<string, number> = {};
  values.forEach((v) => (counts[v] = (counts[v] || 0) + 1));
  return Object.entries(counts).sort((a, b) => b[1] - a[1])[0][0];
}

export async function POST(request: NextRequest) {
  if (!API_KEY) {
    return NextResponse.json({ success: false, error: "OPENROUTER_API_KEY not configured" }, { status: 500 });
  }

  const body = await request.json();
  const { smiles, name, category } = body;

  const userPrompt = `## 시뮬레이션 대상 분자
- 이름: ${name || "신규 제안 소재"}
- SMILES: ${smiles}
- 카테고리: ${category || "미정"}

위 분자에 대해 다음을 시뮬레이션하세요:
1. 분자동역학: 유리전이온도, 분해온도, 안정성 등급
2. 전자구조: HOMO, LUMO, 밴드갭, 전자이동도
3. 합성 경로: 수율, 단계, 난이도, 스케일업 가능성
4. 환경 평가: 독성, 생분해성
5. 종합 성능 점수`;

  const rawResults = await Promise.all(
    MODELS.map(async (m) => {
      const r = await callModel(m.id, userPrompt);
      return { model: m.key, display_name: m.name, color: m.color, result: r };
    })
  );

  const sims = rawResults
    .filter((r) => r.result.success && r.result.data)
    .map((r) => (r.result.data as Record<string, unknown>).simulation_results as Record<string, unknown>)
    .filter(Boolean);

  const averaged = {
    molecular_dynamics: {
      glass_transition_temp: avgField(sims, ["molecular_dynamics", "glass_transition_temp"]),
      decomposition_temp: avgField(sims, ["molecular_dynamics", "decomposition_temp"]),
      stability_rating: modeField(sims, ["molecular_dynamics", "stability_rating"]),
    },
    electronic_properties: {
      homo_energy: avgField(sims, ["electronic_properties", "homo_energy"]),
      lumo_energy: avgField(sims, ["electronic_properties", "lumo_energy"]),
      bandgap: avgField(sims, ["electronic_properties", "bandgap"]),
      electron_mobility: avgField(sims, ["electronic_properties", "electron_mobility"]),
    },
    synthesis: {
      estimated_yield: avgField(sims, ["synthesis_simulation", "estimated_yield"]),
      estimated_steps: avgField(sims, ["synthesis_simulation", "estimated_steps"]),
      difficulty: modeField(sims, ["synthesis_simulation", "difficulty"]),
      scale_up_feasibility: modeField(sims, ["synthesis_simulation", "scale_up_feasibility"]),
    },
    environment: {
      toxicity_risk: modeField(sims, ["environmental_assessment", "toxicity_risk"]),
      biodegradability: modeField(sims, ["environmental_assessment", "biodegradability"]),
    },
    performance: {
      overall_score: avgField(sims, ["performance_prediction", "overall_score"]),
      confidence: avgField(sims, ["performance_prediction", "confidence"]),
    },
  };

  const modelDetails = rawResults.map((r) => ({
    model: r.model,
    display_name: r.display_name,
    color: r.color,
    success: r.result.success,
    error: r.result.error,
    elapsed: r.result.elapsed,
    simulation: r.result.success ? (r.result.data as Record<string, unknown>).simulation_results : undefined,
    comparison: r.result.success ? (r.result.data as Record<string, unknown>).comparison_with_existing : undefined,
    recommendation: r.result.success ? (r.result.data as Record<string, unknown>).recommendation : undefined,
  }));

  return NextResponse.json({
    success: true,
    data: {
      smiles, name, category,
      averaged_prediction: averaged,
      model_details: modelDetails,
      models_succeeded: rawResults.filter((r) => r.result.success).length,
      total_models: rawResults.length,
    },
  });
}
