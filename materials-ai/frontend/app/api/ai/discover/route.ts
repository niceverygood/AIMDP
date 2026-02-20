import { NextRequest, NextResponse } from "next/server";

const OPENROUTER_URL = "https://openrouter.ai/api/v1/chat/completions";
const API_KEY = process.env.OPENROUTER_API_KEY || "";

const MODELS = [
  { key: "claude", id: "anthropic/claude-sonnet-4", name: "Claude Opus 4", color: "#D97706" },
  { key: "gemini", id: "google/gemini-2.5-pro-preview-06-05", name: "Gemini 2.5 Pro", color: "#4285F4" },
  { key: "gpt", id: "openai/gpt-4o", name: "GPT-4o", color: "#10A37F" },
];

const SYSTEM_PROMPT = `당신은 세계 최고 수준의 소재과학자이자 계산화학 전문가입니다.
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
        "thermal_stability": 숫자,
        "dielectric_constant": 숫자,
        "bandgap": 숫자,
        "solubility": 숫자,
        "density": 숫자
      },
      "novelty_score": 0.0~1.0,
      "synthesizability": "easy|moderate|difficult|very_difficult",
      "estimated_cost_level": "low|medium|high|very_high",
      "key_building_blocks": ["출발물질1", "출발물질2"],
      "target_applications": ["응용분야1", "응용분야2"]
    }
  ],
  "data_insights": "기존 데이터에서 발견한 핵심 패턴 (3-4문장)",
  "research_direction": "향후 연구 방향 제안 (2-3문장)"
}

proposed_materials는 3~5개를 제안하세요. 각 분자의 SMILES는 반드시 유효한 화학 구조여야 합니다.`;

// Curated OLED/battery/semiconductor reference data
const REFERENCE_DATA: Record<string, string> = {
  OLED: `BPhen | c1ccc(-c2ccnc3c2ccc2c(-c4ccccc4)ccnc23)cc1 | 열안정성: 380°C | 밴드갭: 3.5eV | 유전율: 3.2
TPBi | c1ccc(-n2c(-c3cc(-c4nc5ccccc5[nH]4)cc(-c4nc5ccccc5[nH]4)c3)nc3ccccc32)cc1 | 열안정성: 420°C | 밴드갭: 3.2eV
Ir(ppy)3 | c1ccc(-c2ccccn2)cc1 | 열안정성: 450°C | 밴드갭: 2.4eV | OLED 인광 발광
TCTA | c1ccc(-n2c3ccccc3c3ccccc32)cc1 | 열안정성: 395°C | 밴드갭: 3.4eV
CBP | c1ccc2c(c1)[nH]c1ccccc12 | 열안정성: 365°C | 밴드갭: 3.6eV
NPB | c1ccc(N(c2ccccc2)c2ccc(-c3ccc(N(c4ccccc4)c4ccc5ccccc5c4)cc3)cc2)cc1 | 열안정성: 380°C
DPEPO | O=P(c1ccccc1)(c1ccccc1)c1ccccc1 | 열안정성: 350°C | 밴드갭: 4.1eV
Alq3 | Oc1cccc2cccnc12.[Al] | 열안정성: 412°C | 밴드갭: 2.8eV`,
  battery: `EC | O=C1OCCO1 | 열안정성: 248°C | 유전율: 89.8 | 배터리 전해질
DMC | COC(=O)OC | 열안정성: 90°C | 유전율: 3.1
FEC | O=C1OCC(F)O1 | 열안정성: 212°C | 유전율: 78.4 | SEI 형성 첨가제
VC | O=c1occo1 | 열안정성: 162°C | 유전율: 126.0
LiPF6 | F[P-](F)(F)(F)(F)F.[Li+] | 열안정성: 200°C | 리튬염
EMC | CCOC(=O)OC | 열안정성: 107°C | 유전율: 2.9`,
  semiconductor: `Polyimide | O=c1oc(=O)c2cc3c(=O)oc(=O)c3cc12 | 열안정성: 480°C | 유전율: 3.4
SU-8 | c1cc(CC2CO2)c(CC2CO2)cc1Cc1cc(CC2CO2)c(CC2CO2)cc1 | 열안정성: 300°C
BCB | C1=CC2=CC=CC2=C1 | 열안정성: 350°C | 유전율: 2.65
PTFE monomer | FC(F)=C(F)F | 열안정성: 327°C | 유전율: 2.1`,
  hard_coating: `TEOS | CCO[Si](OCC)(OCC)OCC | 열안정성: 165°C | 유전율: 3.9
GPTMS | CO[Si](OC)(OC)CCCOCC1CO1 | 열안정성: 200°C
TMSPM | CO[Si](OC)(OC)CCCOC(=O)C(=C)C | 열안정성: 220°C`,
  display: `5CB | CCCCCC1=CC=C(C#N)C=C1 | 열안정성: 168°C | 유전율: 11.0 | 액정
7CB | CCCCCCCC1=CC=C(C#N)C=C1 | 열안정성: 185°C | 유전율: 12.5`,
};

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
        temperature: 0.4,
        max_tokens: 3000,
      }),
    });

    const elapsed = (Date.now() - start) / 1000;
    if (!res.ok) return { success: false, error: `HTTP ${res.status}`, elapsed };

    const data = await res.json();
    let content = data.choices?.[0]?.message?.content || "";
    if (content.includes("```json")) content = content.split("```json")[1].split("```")[0];
    else if (content.includes("```")) content = content.split("```")[1].split("```")[0];

    const parsed = JSON.parse(content.trim());
    return { success: true, data: parsed, elapsed: Math.round(elapsed * 100) / 100 };
  } catch (e) {
    return { success: false, error: String(e).slice(0, 100), elapsed: (Date.now() - start) / 1000 };
  }
}

export async function POST(request: NextRequest) {
  if (!API_KEY) {
    return NextResponse.json({ success: false, error: "OPENROUTER_API_KEY not configured" }, { status: 500 });
  }

  const body = await request.json();
  const category = body.target_category || "OLED";
  const constraints = body.constraints || "";

  const refData = REFERENCE_DATA[category] || REFERENCE_DATA["OLED"];
  const userPrompt = `## 동진쎄미켐 기존 연구 데이터 (${category} 카테고리)
${refData}

## 목표 조건
- 카테고리: ${category}
${constraints ? `- 추가 조건: ${constraints}` : "- 기존 소재보다 성능이 우수한 새로운 분자 구조"}

위 연구 데이터를 분석하여:
1. 구조-물성 관계(SPR) 패턴을 파악하세요
2. 기존 데이터에 없는 새로운 분자 구조를 3~5개 제안하세요
3. 각 제안 분자의 예상 물성, 합성 가능성, 비용을 예측하세요`;

  // Call 3 models in parallel
  const rawResults = await Promise.all(
    MODELS.map(async (m) => {
      const r = await callModel(m.id, userPrompt);
      return { model: m.key, display_name: m.name, color: m.color, result: r };
    })
  );

  // Merge proposals
  const allProposals: Record<string, unknown>[] = [];
  const modelSummaries = rawResults.map((r) => {
    const entry: Record<string, unknown> = {
      model: r.model, display_name: r.display_name, color: r.color,
      success: r.result.success, error: r.result.error, elapsed: r.result.elapsed,
    };
    if (r.result.success && r.result.data) {
      const d = r.result.data as Record<string, unknown>;
      const proposals = (d.proposed_materials || []) as Record<string, unknown>[];
      proposals.forEach((p) => {
        p.proposed_by = r.display_name;
        p.model_color = r.color;
      });
      allProposals.push(...proposals);
      entry.data_insights = d.data_insights;
      entry.research_direction = d.research_direction;
      entry.proposal_count = proposals.length;
    } else {
      entry.proposal_count = 0;
    }
    return entry;
  });

  // Deduplicate by SMILES
  const unique: Record<string, Record<string, unknown>> = {};
  for (const p of allProposals) {
    const smiles = String(p.smiles || "").trim();
    if (!smiles) continue;
    const key = smiles.toLowerCase();
    if (!unique[key]) {
      unique[key] = { ...p, votes: 1, proposed_by_list: [p.proposed_by] };
    } else {
      (unique[key].votes as number)++;
      (unique[key].proposed_by_list as string[]).push(String(p.proposed_by));
    }
  }

  const merged = Object.values(unique).sort(
    (a, b) => (b.votes as number) - (a.votes as number) || (b.novelty_score as number || 0) - (a.novelty_score as number || 0)
  );

  return NextResponse.json({
    success: true,
    data: {
      proposals: merged.slice(0, 8),
      total_proposals: allProposals.length,
      unique_proposals: merged.length,
      model_summaries: modelSummaries,
    },
  });
}
