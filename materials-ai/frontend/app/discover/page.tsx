"use client";

import { useState } from "react";
import { motion, AnimatePresence } from "framer-motion";
import { Button } from "@/components/ui/button";
import { Card, CardContent, CardHeader, CardTitle } from "@/components/ui/card";
import { Badge } from "@/components/ui/badge";
import { Input } from "@/components/ui/input";
import { Label } from "@/components/ui/label";
import { Select, SelectContent, SelectItem, SelectTrigger, SelectValue } from "@/components/ui/select";
import { HelpTooltip } from "@/components/HelpTooltip";
import { Tabs, TabsContent, TabsList, TabsTrigger } from "@/components/ui/tabs";

const API = process.env.NEXT_PUBLIC_API_URL || "http://localhost:8000";

interface Proposal {
  name: string;
  smiles: string;
  design_rationale: string;
  predicted_properties: Record<string, number>;
  novelty_score: number;
  synthesizability: string;
  estimated_cost_level: string;
  key_building_blocks: string[];
  target_applications: string[];
  proposed_by: string;
  model_color: string;
  votes?: number;
  proposed_by_list?: string[];
}

interface SimResult {
  averaged_prediction: {
    molecular_dynamics: Record<string, unknown>;
    electronic_properties: Record<string, unknown>;
    synthesis: Record<string, unknown>;
    environment: Record<string, unknown>;
    performance: Record<string, unknown>;
  };
  model_details: Array<{
    display_name: string;
    color: string;
    success: boolean;
    elapsed: number;
    simulation?: Record<string, unknown>;
    comparison?: string;
    recommendation?: string;
  }>;
  models_succeeded: number;
}

const SYNTH_LABELS: Record<string, { label: string; color: string }> = {
  easy: { label: "쉬움", color: "#00E676" },
  moderate: { label: "보통", color: "#00B4D8" },
  difficult: { label: "어려움", color: "#FF9100" },
  very_difficult: { label: "매우 어려움", color: "#FF5252" },
};

const COST_LABELS: Record<string, { label: string; color: string }> = {
  low: { label: "저비용", color: "#00E676" },
  medium: { label: "중간", color: "#00B4D8" },
  high: { label: "고비용", color: "#FF9100" },
  very_high: { label: "매우 높음", color: "#FF5252" },
};

export default function DiscoverPage() {
  const [category, setCategory] = useState("OLED");
  const [constraints, setConstraints] = useState("");
  const [discovering, setDiscovering] = useState(false);
  const [proposals, setProposals] = useState<Proposal[]>([]);
  const [modelSummaries, setModelSummaries] = useState<unknown[]>([]);

  const [simulating, setSimulating] = useState<string | null>(null);
  const [simResults, setSimResults] = useState<Record<string, SimResult>>({});

  const handleDiscover = async () => {
    setDiscovering(true);
    setProposals([]);
    try {
      const r = await fetch(`${API}/api/ai/discover`, {
        method: "POST",
        headers: { "Content-Type": "application/json" },
        body: JSON.stringify({
          target_category: category,
          target_properties: {},
          constraints,
          data_limit: 30,
        }),
      });
      const d = await r.json();
      if (d.success && d.data) {
        setProposals(d.data.proposals || []);
        setModelSummaries(d.data.model_summaries || []);
      }
    } catch (err) {
      console.error(err);
    } finally {
      setDiscovering(false);
    }
  };

  const handleSimulate = async (p: Proposal) => {
    const key = p.smiles;
    setSimulating(key);
    try {
      const r = await fetch(`${API}/api/ai/simulate`, {
        method: "POST",
        headers: { "Content-Type": "application/json" },
        body: JSON.stringify({ smiles: p.smiles, name: p.name, category }),
      });
      const d = await r.json();
      if (d.success && d.data) {
        setSimResults((prev) => ({ ...prev, [key]: d.data }));
      }
    } catch (err) {
      console.error(err);
    } finally {
      setSimulating(null);
    }
  };

  return (
    <div className="p-8 space-y-8 grid-pattern min-h-screen">
      {/* Header */}
      <motion.div initial={{ opacity: 0, y: -20 }} animate={{ opacity: 1, y: 0 }}>
        <div className="flex items-center gap-2">
          <h1 className="text-2xl font-bold text-white">AI 신소재 발굴</h1>
          <HelpTooltip
            title="AI 신소재 발굴 + 시뮬레이션"
            description="3개 AI 모델이 기존 연구 데이터를 학습하여 새로운 분자 구조를 제안하고, 시뮬레이션으로 성능을 예측합니다."
            details={[
              "1단계 — 발굴: 기존 연구 데이터의 구조-물성 관계를 분석하여 새로운 분자 구조 제안",
              "2단계 — 시뮬레이션: 제안 소재의 열역학/전자구조/합성경로/환경영향 시뮬레이션",
              "3개 모델이 독립 분석 → 교차검증 → 합의 결과 도출",
            ]}
            size="md"
          />
        </div>
        <p className="text-sm text-[#8892B0] mt-1">기존 연구 데이터 학습 → AI 신소재 제안 → 시뮬레이션 예측</p>
      </motion.div>

      {/* Discovery Config */}
      <Card className="bg-[#0D1B2A] border-white/[0.06]">
        <CardContent className="p-6">
          <div className="flex flex-wrap items-end gap-4">
            <div className="space-y-1.5">
              <Label className="text-xs text-[#8892B0]">목표 카테고리</Label>
              <Select value={category} onValueChange={setCategory}>
                <SelectTrigger className="w-40 bg-[#112240] border-white/[0.08] text-white">
                  <SelectValue />
                </SelectTrigger>
                <SelectContent className="bg-[#112240] border-white/[0.08]">
                  <SelectItem value="OLED" className="text-white">OLED</SelectItem>
                  <SelectItem value="battery" className="text-white">배터리</SelectItem>
                  <SelectItem value="semiconductor" className="text-white">반도체</SelectItem>
                  <SelectItem value="hard_coating" className="text-white">하드코팅</SelectItem>
                  <SelectItem value="display" className="text-white">디스플레이</SelectItem>
                  <SelectItem value="all" className="text-white">전체</SelectItem>
                </SelectContent>
              </Select>
            </div>

            <div className="space-y-1.5 flex-1 min-w-[200px]">
              <Label className="text-xs text-[#8892B0]">추가 조건 (선택)</Label>
              <Input
                value={constraints}
                onChange={(e) => setConstraints(e.target.value)}
                placeholder="예: 열안정성 400°C 이상, 합성 난이도 낮을 것"
                className="bg-[#112240] border-white/[0.08] text-white text-sm placeholder:text-[#8892B0]/40"
              />
            </div>

            <Button
              onClick={handleDiscover}
              disabled={discovering}
              className="bg-gradient-to-r from-[#D97706] via-[#4285F4] to-[#10A37F] text-white"
            >
              {discovering ? (
                <div className="flex items-center gap-2">
                  <div className="w-4 h-4 border-2 border-white border-t-transparent rounded-full animate-spin" />
                  3개 AI 분석 중...
                </div>
              ) : (
                <div className="flex items-center gap-2">
                  <svg className="w-4 h-4" fill="none" viewBox="0 0 24 24" stroke="currentColor" strokeWidth={2}>
                    <path strokeLinecap="round" strokeLinejoin="round" d="M9.813 15.904 9 18.75l-.813-2.846a4.5 4.5 0 0 0-3.09-3.09L2.25 12l2.846-.813a4.5 4.5 0 0 0 3.09-3.09L9 5.25l.813 2.846a4.5 4.5 0 0 0 3.09 3.09L15.75 12l-2.846.813a4.5 4.5 0 0 0-3.09 3.09Z" />
                  </svg>
                  신소재 발굴 시작
                </div>
              )}
            </Button>
          </div>
        </CardContent>
      </Card>

      {/* Loading */}
      <AnimatePresence>
        {discovering && (
          <motion.div initial={{ opacity: 0 }} animate={{ opacity: 1 }} exit={{ opacity: 0 }} className="space-y-3">
            {["Claude Opus 4", "Gemini 2.5 Pro", "GPT-4o"].map((name, i) => (
              <motion.div key={name} initial={{ x: -20, opacity: 0 }} animate={{ x: 0, opacity: 1 }} transition={{ delay: i * 0.2 }}
                className="flex items-center gap-3 p-3 rounded-lg bg-[#0D1B2A] border border-white/[0.06]"
              >
                <span>{["🟠", "🔵", "🟢"][i]}</span>
                <span className="text-xs text-white">{name}</span>
                <div className="flex-1 text-[10px] text-[#8892B0]">기존 데이터 학습 → 구조-물성 패턴 분석 → 신소재 설계 중...</div>
                <div className="w-4 h-4 border-2 border-[#00B4D8] border-t-transparent rounded-full animate-spin" />
              </motion.div>
            ))}
          </motion.div>
        )}
      </AnimatePresence>

      {/* Proposals */}
      {proposals.length > 0 && (
        <motion.div initial={{ opacity: 0, y: 20 }} animate={{ opacity: 1, y: 0 }} className="space-y-4">
          <div className="flex items-center gap-2">
            <h2 className="text-lg font-bold text-white">AI 제안 신소재</h2>
            <Badge variant="outline" className="text-[10px] border-[#00E676]/30 text-[#00E676]">
              {proposals.length}개 제안
            </Badge>
          </div>

          <div className="grid grid-cols-1 lg:grid-cols-2 gap-4">
            {proposals.map((p, idx) => {
              const sim = simResults[p.smiles];
              const synthInfo = SYNTH_LABELS[p.synthesizability] || { label: p.synthesizability, color: "#8892B0" };
              const costInfo = COST_LABELS[p.estimated_cost_level] || { label: p.estimated_cost_level, color: "#8892B0" };

              return (
                <Card key={idx} className="bg-[#0D1B2A] border-white/[0.06] hover:border-[#00B4D8]/20 transition-all">
                  <CardContent className="p-5 space-y-3">
                    {/* Header */}
                    <div className="flex items-start justify-between">
                      <div>
                        <div className="flex items-center gap-2">
                          <div className="w-7 h-7 rounded-lg bg-[#00B4D8]/10 flex items-center justify-center text-xs font-bold text-[#00B4D8]">
                            {idx + 1}
                          </div>
                          <h3 className="text-sm font-bold text-white">{p.name}</h3>
                        </div>
                        <code className="text-[10px] text-[#8892B0] font-mono mt-1 block">
                          {p.smiles.length > 50 ? p.smiles.slice(0, 50) + "..." : p.smiles}
                        </code>
                      </div>
                      <div className="text-right">
                        <div className="text-lg font-bold font-mono text-[#00E676]">
                          {((p.novelty_score || 0) * 100).toFixed(0)}%
                        </div>
                        <p className="text-[9px] text-[#8892B0]">신규성</p>
                      </div>
                    </div>

                    {/* Proposing model */}
                    <div className="flex flex-wrap gap-1">
                      {(p.proposed_by_list || [p.proposed_by]).map((m, i) => (
                        <Badge key={i} variant="outline" className="text-[9px]"
                          style={{ borderColor: m.includes("Claude") ? "#D9770640" : m.includes("Gemini") ? "#4285F440" : "#10A37F40",
                                   color: m.includes("Claude") ? "#D97706" : m.includes("Gemini") ? "#4285F4" : "#10A37F" }}>
                          {m}
                        </Badge>
                      ))}
                      <Badge variant="outline" className="text-[9px]" style={{ borderColor: synthInfo.color + "40", color: synthInfo.color }}>
                        합성: {synthInfo.label}
                      </Badge>
                      <Badge variant="outline" className="text-[9px]" style={{ borderColor: costInfo.color + "40", color: costInfo.color }}>
                        비용: {costInfo.label}
                      </Badge>
                    </div>

                    {/* Rationale */}
                    <p className="text-[11px] text-[#8892B0] leading-relaxed">{p.design_rationale}</p>

                    {/* Predicted properties */}
                    {p.predicted_properties && (
                      <div className="grid grid-cols-3 gap-1.5">
                        {Object.entries(p.predicted_properties).map(([k, v]) => (
                          <div key={k} className="p-1.5 rounded bg-[#112240]/50 text-center">
                            <p className="text-[9px] text-[#8892B0]">{k.replace(/_/g, " ")}</p>
                            <p className="text-[10px] font-mono text-white">{typeof v === "number" ? v.toFixed(2) : v}</p>
                          </div>
                        ))}
                      </div>
                    )}

                    {/* Applications */}
                    {p.target_applications?.length > 0 && (
                      <div className="flex flex-wrap gap-1">
                        {p.target_applications.map((a, i) => (
                          <Badge key={i} variant="outline" className="text-[9px] border-[#7C4DFF]/20 text-[#7C4DFF]/80">{a}</Badge>
                        ))}
                      </div>
                    )}

                    {/* Simulate button */}
                    <Button
                      size="sm"
                      variant="outline"
                      className="w-full border-[#00B4D8]/20 text-[#00B4D8] hover:bg-[#00B4D8]/10 text-xs"
                      onClick={() => handleSimulate(p)}
                      disabled={simulating === p.smiles}
                    >
                      {simulating === p.smiles ? (
                        <div className="flex items-center gap-2">
                          <div className="w-3 h-3 border-2 border-[#00B4D8] border-t-transparent rounded-full animate-spin" />
                          시뮬레이션 중...
                        </div>
                      ) : sim ? "시뮬레이션 결과 보기" : "시뮬레이션 실행"}
                    </Button>

                    {/* Simulation Results */}
                    <AnimatePresence>
                      {sim && (
                        <motion.div initial={{ opacity: 0, height: 0 }} animate={{ opacity: 1, height: "auto" }} className="space-y-2 pt-2 border-t border-white/[0.04]">
                          <Tabs defaultValue="thermo" className="w-full">
                            <TabsList className="bg-[#112240]/50 w-full">
                              <TabsTrigger value="thermo" className="text-[10px] flex-1">열역학</TabsTrigger>
                              <TabsTrigger value="electronic" className="text-[10px] flex-1">전자구조</TabsTrigger>
                              <TabsTrigger value="synth" className="text-[10px] flex-1">합성</TabsTrigger>
                              <TabsTrigger value="env" className="text-[10px] flex-1">환경</TabsTrigger>
                            </TabsList>

                            <TabsContent value="thermo" className="space-y-1 mt-2">
                              <SimRow label="유리전이온도" value={sim.averaged_prediction.molecular_dynamics.glass_transition_temp} unit="°C" />
                              <SimRow label="분해온도" value={sim.averaged_prediction.molecular_dynamics.decomposition_temp} unit="°C" />
                              <SimRow label="안정성 등급" value={sim.averaged_prediction.molecular_dynamics.stability_rating} />
                            </TabsContent>

                            <TabsContent value="electronic" className="space-y-1 mt-2">
                              <SimRow label="HOMO" value={sim.averaged_prediction.electronic_properties.homo_energy} unit="eV" />
                              <SimRow label="LUMO" value={sim.averaged_prediction.electronic_properties.lumo_energy} unit="eV" />
                              <SimRow label="밴드갭" value={sim.averaged_prediction.electronic_properties.bandgap} unit="eV" />
                              <SimRow label="전자 이동도" value={sim.averaged_prediction.electronic_properties.electron_mobility} unit="cm²/V·s" />
                            </TabsContent>

                            <TabsContent value="synth" className="space-y-1 mt-2">
                              <SimRow label="예상 수율" value={sim.averaged_prediction.synthesis.estimated_yield} unit="%" />
                              <SimRow label="합성 단계" value={sim.averaged_prediction.synthesis.estimated_steps} unit="단계" />
                              <SimRow label="합성 난이도" value={sim.averaged_prediction.synthesis.difficulty} />
                              <SimRow label="스케일업" value={sim.averaged_prediction.synthesis.scale_up_feasibility} />
                            </TabsContent>

                            <TabsContent value="env" className="space-y-1 mt-2">
                              <SimRow label="독성 위험" value={sim.averaged_prediction.environment.toxicity_risk} />
                              <SimRow label="생분해성" value={sim.averaged_prediction.environment.biodegradability} />
                            </TabsContent>
                          </Tabs>

                          {/* Model recommendations */}
                          <div className="space-y-1.5 pt-2">
                            {sim.model_details.filter((m) => m.success).map((m, i) => (
                              <div key={i} className="p-2 rounded bg-[#112240]/30 text-[10px]">
                                <span style={{ color: m.color }} className="font-medium">{m.display_name}:</span>{" "}
                                <span className="text-[#8892B0]">{String((m as Record<string, unknown>).recommendation || (m as Record<string, unknown>).comparison || "")}</span>
                              </div>
                            ))}
                          </div>

                          {/* Confidence */}
                          <div className="flex items-center justify-between p-2 rounded bg-[#00E676]/5 border border-[#00E676]/10">
                            <span className="text-[10px] text-[#8892B0]">종합 신뢰도 ({sim.models_succeeded}/3 모델)</span>
                            <span className="text-sm font-bold font-mono text-[#00E676]">
                              {((sim.averaged_prediction.performance.overall_score as number || 0) * 100).toFixed(1)}%
                            </span>
                          </div>
                        </motion.div>
                      )}
                    </AnimatePresence>
                  </CardContent>
                </Card>
              );
            })}
          </div>
        </motion.div>
      )}
    </div>
  );
}

function SimRow({ label, value, unit }: { label: string; value: unknown; unit?: string }) {
  const display = value !== null && value !== undefined ? `${value}${unit ? ` ${unit}` : ""}` : "—";
  return (
    <div className="flex justify-between px-2 py-1 rounded bg-[#0A1628]/50">
      <span className="text-[10px] text-[#8892B0]">{label}</span>
      <span className="text-[10px] font-mono text-white">{display}</span>
    </div>
  );
}
