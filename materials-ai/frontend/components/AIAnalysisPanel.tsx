"use client";

import { useState } from "react";
import { motion, AnimatePresence } from "framer-motion";
import { Button } from "@/components/ui/button";
import { Card, CardContent } from "@/components/ui/card";
import { Badge } from "@/components/ui/badge";
import { HelpTooltip } from "@/components/HelpTooltip";

const API_URL = process.env.NEXT_PUBLIC_API_URL || "http://localhost:8000";

interface ModelResult {
  model: string;
  display_name: string;
  color: string;
  strength: string;
  success: boolean;
  error?: string;
  elapsed_sec: number;
  top_picks: Array<{
    rank: number;
    name: string;
    smiles: string;
    confidence: number;
    reasoning: string;
    pros: string[];
    cons: string[];
    applications: string[];
  }>;
  overall_analysis: string;
}

interface ConsensusItem {
  name: string;
  smiles: string;
  votes: number;
  avg_confidence: number;
  ensemble_score: number;
  reasonings: Array<{ model: string; text: string; confidence: number }>;
  recommending_models: string[];
  pros: string[];
  cons: string[];
  applications: string[];
}

interface EnsembleResult {
  consensus: ConsensusItem[];
  agreement_score: number;
  models_agreed: number;
  total_models: number;
  overall_summary: string;
  risk_factors: string[];
}

interface AIAnalysisData {
  model_results: ModelResult[];
  ensemble: EnsembleResult;
  total_elapsed_sec: number;
}

interface AIAnalysisPanelProps {
  searchQuery?: {
    category?: string;
    min_thermal_stability?: number;
    min_bandgap?: number;
    min_dielectric?: number;
    smiles?: string;
  };
}

const MODEL_ICONS: Record<string, string> = {
  claude: "🟠",
  gemini: "🔵",
  gpt: "🟢",
};

export function AIAnalysisPanel({ searchQuery }: AIAnalysisPanelProps) {
  const [analysisData, setAnalysisData] = useState<AIAnalysisData | null>(null);
  const [loading, setLoading] = useState(false);
  const [error, setError] = useState<string | null>(null);

  const handleAnalyze = async () => {
    setLoading(true);
    setError(null);

    try {
      const res = await fetch(`${API_URL}/api/ai/analyze`, {
        method: "POST",
        headers: { "Content-Type": "application/json" },
        body: JSON.stringify({
          category: searchQuery?.category,
          min_thermal_stability: searchQuery?.min_thermal_stability,
          min_bandgap: searchQuery?.min_bandgap,
          min_dielectric: searchQuery?.min_dielectric,
          smiles: searchQuery?.smiles,
          limit: 10,
        }),
      });

      if (!res.ok) throw new Error(`API error ${res.status}`);

      const data = await res.json();
      if (data.success) {
        setAnalysisData(data.data);
      } else {
        setError(data.error || "분석 실패");
      }
    } catch (err) {
      setError(err instanceof Error ? err.message : "네트워크 오류 — 백엔드 API가 실행 중인지 확인하세요");
    } finally {
      setLoading(false);
    }
  };

  return (
    <div className="space-y-4">
      {/* Header + Trigger Button */}
      <div className="flex items-center justify-between">
        <div className="flex items-center gap-2">
          <div className="flex -space-x-1">
            <span className="text-base">🟠</span>
            <span className="text-base">🔵</span>
            <span className="text-base">🟢</span>
          </div>
          <h3 className="text-sm font-bold text-white">AI 트리플 앙상블 분석</h3>
          <HelpTooltip
            title="트리플 AI 앙상블 분석"
            description="3개의 최상위 AI 모델이 동시에 후보 소재를 분석하고, 교차 검증을 통해 합의된 추천을 제공합니다."
            details={[
              "Claude Opus 4: 분자 구조 분석 및 합성 가능성 추론",
              "Gemini 2.5 Pro: 학술 논문 기반 물성 데이터 교차검증",
              "GPT-4o: 특허 데이터 분석 및 응용 분야 매칭",
              "3개 모델이 모두 동의한 소재일수록 신뢰도가 높습니다.",
            ]}
          />
        </div>
        <Button
          onClick={handleAnalyze}
          disabled={loading}
          size="sm"
          className="bg-gradient-to-r from-[#D97706] via-[#4285F4] to-[#10A37F] hover:opacity-90 text-white text-xs"
        >
          {loading ? (
            <div className="flex items-center gap-2">
              <div className="w-3 h-3 border-2 border-white border-t-transparent rounded-full animate-spin" />
              <span>3개 AI 분석 중...</span>
            </div>
          ) : (
            "AI 앙상블 분석 시작"
          )}
        </Button>
      </div>

      {error && (
        <div className="p-3 rounded-lg bg-[#FF5252]/10 border border-[#FF5252]/20 text-xs text-[#FF5252]">
          {error}
        </div>
      )}

      {/* Loading Animation */}
      <AnimatePresence>
        {loading && (
          <motion.div
            initial={{ opacity: 0, height: 0 }}
            animate={{ opacity: 1, height: "auto" }}
            exit={{ opacity: 0, height: 0 }}
            className="space-y-3"
          >
            {[
              { name: "Claude Opus 4", color: "#D97706", icon: "🟠" },
              { name: "Gemini 2.5 Pro", color: "#4285F4", icon: "🔵" },
              { name: "GPT-4o", color: "#10A37F", icon: "🟢" },
            ].map((model, i) => (
              <motion.div
                key={model.name}
                initial={{ opacity: 0, x: -20 }}
                animate={{ opacity: 1, x: 0 }}
                transition={{ delay: i * 0.3 }}
                className="flex items-center gap-3 p-3 rounded-lg border border-white/[0.06] bg-[#112240]/50"
              >
                <span>{model.icon}</span>
                <span className="text-xs text-white font-medium">{model.name}</span>
                <div className="flex-1 h-1 bg-white/[0.04] rounded-full overflow-hidden">
                  <motion.div
                    className="h-full rounded-full"
                    style={{ backgroundColor: model.color }}
                    initial={{ width: "0%" }}
                    animate={{ width: "100%" }}
                    transition={{ duration: 8 + i * 3, ease: "easeInOut" }}
                  />
                </div>
                <span className="text-[10px] text-[#8892B0]">분석 중...</span>
              </motion.div>
            ))}
          </motion.div>
        )}
      </AnimatePresence>

      {/* Results */}
      <AnimatePresence>
        {analysisData && !loading && (
          <motion.div
            initial={{ opacity: 0, y: 20 }}
            animate={{ opacity: 1, y: 0 }}
            exit={{ opacity: 0 }}
            className="space-y-4"
          >
            {/* Model Status Bar */}
            <div className="flex items-center gap-2 flex-wrap">
              {analysisData.model_results.map((r) => (
                <Badge
                  key={r.model}
                  variant="outline"
                  className="text-[10px] gap-1"
                  style={{
                    borderColor: r.success ? r.color + "40" : "#FF525240",
                    color: r.success ? r.color : "#FF5252",
                  }}
                >
                  {MODEL_ICONS[r.model]}
                  {r.display_name}
                  {r.success ? ` ✓ ${r.elapsed_sec}s` : " ✗ 실패"}
                </Badge>
              ))}
              <Badge
                variant="outline"
                className="text-[10px] border-white/10 text-[#8892B0]"
              >
                총 {analysisData.total_elapsed_sec}초
              </Badge>
            </div>

            {/* Consensus Agreement Score */}
            {analysisData.ensemble && (
              <Card className="bg-gradient-to-r from-[#D97706]/5 via-[#4285F4]/5 to-[#10A37F]/5 border-white/[0.08]">
                <CardContent className="p-4">
                  <div className="flex items-center justify-between mb-3">
                    <div className="flex items-center gap-2">
                      <div className="w-8 h-8 rounded-lg bg-gradient-to-br from-[#D97706] via-[#4285F4] to-[#10A37F] flex items-center justify-center text-white text-xs font-bold">
                        AI
                      </div>
                      <div>
                        <p className="text-xs font-bold text-white">
                          모델 합의율 {analysisData.ensemble.agreement_score}%
                        </p>
                        <p className="text-[10px] text-[#8892B0]">
                          {analysisData.ensemble.models_agreed}/{analysisData.ensemble.total_models} 모델 응답 성공
                        </p>
                      </div>
                    </div>
                    <div
                      className="text-2xl font-bold font-mono"
                      style={{
                        color:
                          analysisData.ensemble.agreement_score >= 80
                            ? "#00E676"
                            : analysisData.ensemble.agreement_score >= 50
                              ? "#00B4D8"
                              : "#FF9100",
                      }}
                    >
                      {analysisData.ensemble.agreement_score}%
                    </div>
                  </div>

                  {/* Agreement bar */}
                  <div className="w-full h-2 bg-white/[0.04] rounded-full overflow-hidden mb-3">
                    <div
                      className="h-full rounded-full bg-gradient-to-r from-[#D97706] via-[#4285F4] to-[#10A37F] transition-all duration-1000"
                      style={{ width: `${analysisData.ensemble.agreement_score}%` }}
                    />
                  </div>

                  {/* Top Consensus Picks */}
                  <div className="space-y-3">
                    {analysisData.ensemble.consensus.slice(0, 3).map((item, idx) => (
                      <div
                        key={item.name}
                        className="p-3 rounded-lg bg-[#0A1628]/50 border border-white/[0.04]"
                      >
                        <div className="flex items-start justify-between mb-2">
                          <div className="flex items-center gap-2">
                            <div className="w-6 h-6 rounded bg-[#00B4D8]/10 flex items-center justify-center text-xs font-bold text-[#00B4D8]">
                              {idx + 1}
                            </div>
                            <div>
                              <p className="text-xs font-bold text-white">{item.name}</p>
                              <code className="text-[9px] text-[#8892B0] font-mono">
                                {item.smiles?.slice(0, 35)}
                                {(item.smiles?.length || 0) > 35 ? "..." : ""}
                              </code>
                            </div>
                          </div>
                          <div className="text-right">
                            <p className="text-sm font-bold font-mono text-[#00E676]">
                              {(item.ensemble_score * 100).toFixed(1)}%
                            </p>
                            <p className="text-[9px] text-[#8892B0]">
                              {item.votes}개 모델 추천
                            </p>
                          </div>
                        </div>

                        {/* Model voting indicators */}
                        <div className="flex gap-1 mb-2">
                          {item.recommending_models.map((m) => (
                            <Badge
                              key={m}
                              className="text-[9px] px-1.5 py-0 bg-white/[0.04]"
                              variant="outline"
                              style={{
                                borderColor:
                                  m.includes("Claude") ? "#D9770640" :
                                  m.includes("Gemini") ? "#4285F440" : "#10A37F40",
                                color:
                                  m.includes("Claude") ? "#D97706" :
                                  m.includes("Gemini") ? "#4285F4" : "#10A37F",
                              }}
                            >
                              {m.includes("Claude") ? "🟠" : m.includes("Gemini") ? "🔵" : "🟢"} {m}
                            </Badge>
                          ))}
                        </div>

                        {/* Reasoning from each model */}
                        {item.reasonings.map((r, ri) => (
                          <div key={ri} className="flex gap-2 mb-1.5">
                            <span className="text-[10px] text-[#8892B0] shrink-0 mt-0.5">
                              {r.model.includes("Claude") ? "🟠" : r.model.includes("Gemini") ? "🔵" : "🟢"}
                            </span>
                            <p className="text-[10px] text-[#8892B0]/80 leading-relaxed">
                              {r.text}
                            </p>
                          </div>
                        ))}

                        {/* Pros / Cons / Applications */}
                        <div className="flex flex-wrap gap-1 mt-2">
                          {item.pros.slice(0, 3).map((p, i) => (
                            <Badge key={`p${i}`} variant="outline" className="text-[9px] border-[#00E676]/20 text-[#00E676]/80">
                              ✓ {p}
                            </Badge>
                          ))}
                          {item.cons.slice(0, 2).map((c, i) => (
                            <Badge key={`c${i}`} variant="outline" className="text-[9px] border-[#FF9100]/20 text-[#FF9100]/80">
                              ⚠ {c}
                            </Badge>
                          ))}
                        </div>

                        {item.applications.length > 0 && (
                          <div className="flex flex-wrap gap-1 mt-1.5">
                            {item.applications.slice(0, 3).map((a, i) => (
                              <Badge key={`a${i}`} variant="outline" className="text-[9px] border-[#7C4DFF]/20 text-[#7C4DFF]/80">
                                {a}
                              </Badge>
                            ))}
                          </div>
                        )}
                      </div>
                    ))}
                  </div>

                  {/* Risk Factors */}
                  {analysisData.ensemble.risk_factors.length > 0 && (
                    <div className="mt-3 p-2 rounded bg-[#FF9100]/5 border border-[#FF9100]/10">
                      <p className="text-[10px] font-medium text-[#FF9100] mb-1">주의사항</p>
                      {analysisData.ensemble.risk_factors.map((r, i) => (
                        <p key={i} className="text-[10px] text-[#8892B0] leading-relaxed">• {r}</p>
                      ))}
                    </div>
                  )}
                </CardContent>
              </Card>
            )}
          </motion.div>
        )}
      </AnimatePresence>
    </div>
  );
}
