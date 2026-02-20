"use client";

import { useState, useEffect, useCallback } from "react";
import { motion, AnimatePresence } from "framer-motion";
import { Button } from "@/components/ui/button";
import { Card, CardContent, CardHeader, CardTitle } from "@/components/ui/card";
import { Badge } from "@/components/ui/badge";
import { Input } from "@/components/ui/input";
import { Label } from "@/components/ui/label";
import { HelpTooltip } from "@/components/HelpTooltip";
import { supabase } from "@/lib/supabase";

interface FeedbackPanelProps {
  materialName?: string;
  materialSmiles?: string;
  category?: string;
  predictedProps?: {
    thermal_stability?: number | null;
    dielectric_constant?: number | null;
    bandgap?: number | null;
    solubility?: number | null;
    density?: number | null;
  };
  onSubmitted?: () => void;
}

interface FeedbackEntry {
  id: number;
  material_name: string;
  material_smiles: string;
  category: string;
  rating: number;
  experiment_success: boolean;
  actual_thermal_stability: number | null;
  actual_bandgap: number | null;
  pros: string;
  cons: string;
  notes: string;
  researcher: string;
  created_at: string;
}

export function FeedbackPanel({
  materialName,
  materialSmiles,
  category,
  predictedProps,
  onSubmitted,
}: FeedbackPanelProps) {
  const [open, setOpen] = useState(false);
  const [submitting, setSubmitting] = useState(false);
  const [submitted, setSubmitted] = useState(false);

  // Form fields
  const [rating, setRating] = useState(0);
  const [success, setSuccess] = useState<boolean | null>(null);
  const [actualThermal, setActualThermal] = useState("");
  const [actualDielectric, setActualDielectric] = useState("");
  const [actualBandgap, setActualBandgap] = useState("");
  const [actualSolubility, setActualSolubility] = useState("");
  const [actualDensity, setActualDensity] = useState("");
  const [pros, setPros] = useState("");
  const [cons, setCons] = useState("");
  const [notes, setNotes] = useState("");
  const [researcher, setResearcher] = useState("");

  const handleSubmit = async () => {
    if (rating === 0) return;
    setSubmitting(true);

    const feedback = {
      material_name: materialName || "",
      material_smiles: materialSmiles || "",
      category: category || "",
      rating,
      experiment_success: success ?? false,
      actual_thermal_stability: actualThermal ? parseFloat(actualThermal) : null,
      actual_dielectric_constant: actualDielectric ? parseFloat(actualDielectric) : null,
      actual_bandgap: actualBandgap ? parseFloat(actualBandgap) : null,
      actual_solubility: actualSolubility ? parseFloat(actualSolubility) : null,
      actual_density: actualDensity ? parseFloat(actualDensity) : null,
      pros,
      cons,
      notes,
      researcher,
      feedback_type: "evaluation",
    };

    const { error } = await supabase.from("researcher_feedback").insert([feedback]);

    if (error) {
      alert(`저장 실패: ${error.message}\n\nSupabase에 researcher_feedback 테이블을 먼저 생성하세요.\n/api/setup-feedback 에서 SQL을 확인할 수 있습니다.`);
    } else {
      setSubmitted(true);
      onSubmitted?.();
    }

    setSubmitting(false);
  };

  if (submitted) {
    return (
      <motion.div initial={{ opacity: 0 }} animate={{ opacity: 1 }}
        className="p-3 rounded-lg bg-[#00E676]/10 border border-[#00E676]/20 text-xs text-[#00E676] flex items-center gap-2"
      >
        <svg className="w-4 h-4" fill="none" viewBox="0 0 24 24" stroke="currentColor" strokeWidth={2}>
          <path strokeLinecap="round" strokeLinejoin="round" d="M9 12.75 11.25 15 15 9.75M21 12a9 9 0 1 1-18 0 9 9 0 0 1 18 0Z" />
        </svg>
        피드백이 저장되었습니다. 다음 AI 분석에 반영됩니다.
      </motion.div>
    );
  }

  return (
    <div>
      <Button
        onClick={() => setOpen(!open)}
        variant="outline"
        size="sm"
        className="w-full text-[10px] border-[#FF9100]/20 text-[#FF9100] hover:bg-[#FF9100]/10"
      >
        {open ? "접기" : "연구원 평가 입력"}
      </Button>

      <AnimatePresence>
        {open && (
          <motion.div
            initial={{ opacity: 0, height: 0 }}
            animate={{ opacity: 1, height: "auto" }}
            exit={{ opacity: 0, height: 0 }}
            className="overflow-hidden"
          >
            <div className="mt-3 p-4 rounded-xl bg-[#112240]/50 border border-white/[0.06] space-y-3">
              <div className="flex items-center gap-2">
                <span className="text-xs font-medium text-white">연구원 평가</span>
                <HelpTooltip
                  title="연구원 평가"
                  description="이 소재에 대한 실험 결과와 평가를 입력하면, 다음 AI 분석 시 이 데이터가 반영됩니다."
                  details={[
                    "별점: 소재의 전반적인 적합성 평가 (1~5)",
                    "실측 물성: 실험으로 측정한 실제 값 입력",
                    "장단점: 이 소재의 실질적인 장점/단점",
                    "피드백은 다음 AI 신소재 발굴 시 학습 데이터로 활용됩니다.",
                  ]}
                />
              </div>

              {/* Star Rating */}
              <div className="space-y-1">
                <Label className="text-[10px] text-[#8892B0]">적합성 평점</Label>
                <div className="flex gap-1">
                  {[1, 2, 3, 4, 5].map((star) => (
                    <button
                      key={star}
                      onClick={() => setRating(star)}
                      className={`text-lg transition-colors ${
                        star <= rating ? "text-[#FFD600]" : "text-[#8892B0]/30"
                      } hover:text-[#FFD600]`}
                    >
                      ★
                    </button>
                  ))}
                  {rating > 0 && <span className="text-[10px] text-[#8892B0] ml-2 self-center">{rating}/5</span>}
                </div>
              </div>

              {/* Experiment Success */}
              <div className="space-y-1">
                <Label className="text-[10px] text-[#8892B0]">실험 결과</Label>
                <div className="flex gap-2">
                  <button
                    onClick={() => setSuccess(true)}
                    className={`px-3 py-1.5 rounded text-[10px] border transition-all ${
                      success === true
                        ? "bg-[#00E676]/10 border-[#00E676]/30 text-[#00E676]"
                        : "border-white/[0.06] text-[#8892B0] hover:border-white/20"
                    }`}
                  >
                    성공
                  </button>
                  <button
                    onClick={() => setSuccess(false)}
                    className={`px-3 py-1.5 rounded text-[10px] border transition-all ${
                      success === false
                        ? "bg-[#FF5252]/10 border-[#FF5252]/30 text-[#FF5252]"
                        : "border-white/[0.06] text-[#8892B0] hover:border-white/20"
                    }`}
                  >
                    실패
                  </button>
                  <button
                    onClick={() => setSuccess(null)}
                    className={`px-3 py-1.5 rounded text-[10px] border transition-all ${
                      success === null
                        ? "bg-white/[0.04] border-white/20 text-white"
                        : "border-white/[0.06] text-[#8892B0] hover:border-white/20"
                    }`}
                  >
                    미실험
                  </button>
                </div>
              </div>

              {/* Actual Properties */}
              <div className="space-y-1">
                <Label className="text-[10px] text-[#8892B0]">실측 물성 (선택)</Label>
                <div className="grid grid-cols-2 gap-2">
                  <div>
                    <span className="text-[9px] text-[#8892B0]">열안정성 (°C) {predictedProps?.thermal_stability ? `← 예측: ${predictedProps.thermal_stability}` : ""}</span>
                    <Input value={actualThermal} onChange={(e) => setActualThermal(e.target.value)} placeholder="실측값"
                      className="h-7 text-[10px] bg-[#0A1628] border-white/[0.06] text-white" />
                  </div>
                  <div>
                    <span className="text-[9px] text-[#8892B0]">밴드갭 (eV) {predictedProps?.bandgap ? `← 예측: ${predictedProps.bandgap}` : ""}</span>
                    <Input value={actualBandgap} onChange={(e) => setActualBandgap(e.target.value)} placeholder="실측값"
                      className="h-7 text-[10px] bg-[#0A1628] border-white/[0.06] text-white" />
                  </div>
                  <div>
                    <span className="text-[9px] text-[#8892B0]">유전율</span>
                    <Input value={actualDielectric} onChange={(e) => setActualDielectric(e.target.value)} placeholder="실측값"
                      className="h-7 text-[10px] bg-[#0A1628] border-white/[0.06] text-white" />
                  </div>
                  <div>
                    <span className="text-[9px] text-[#8892B0]">용해도</span>
                    <Input value={actualSolubility} onChange={(e) => setActualSolubility(e.target.value)} placeholder="실측값"
                      className="h-7 text-[10px] bg-[#0A1628] border-white/[0.06] text-white" />
                  </div>
                </div>
              </div>

              {/* Pros / Cons */}
              <div className="grid grid-cols-2 gap-2">
                <div>
                  <Label className="text-[10px] text-[#00E676]">장점</Label>
                  <Input value={pros} onChange={(e) => setPros(e.target.value)} placeholder="예: 합성 용이, 높은 수율"
                    className="h-7 text-[10px] bg-[#0A1628] border-white/[0.06] text-white" />
                </div>
                <div>
                  <Label className="text-[10px] text-[#FF5252]">단점</Label>
                  <Input value={cons} onChange={(e) => setCons(e.target.value)} placeholder="예: 결정화 어려움"
                    className="h-7 text-[10px] bg-[#0A1628] border-white/[0.06] text-white" />
                </div>
              </div>

              {/* Notes + Researcher */}
              <div>
                <Label className="text-[10px] text-[#8892B0]">추가 메모</Label>
                <Input value={notes} onChange={(e) => setNotes(e.target.value)} placeholder="자유 형식 메모"
                  className="h-7 text-[10px] bg-[#0A1628] border-white/[0.06] text-white" />
              </div>
              <div>
                <Label className="text-[10px] text-[#8892B0]">연구원 이름</Label>
                <Input value={researcher} onChange={(e) => setResearcher(e.target.value)} placeholder="이름"
                  className="h-7 text-[10px] bg-[#0A1628] border-white/[0.06] text-white" />
              </div>

              {/* Submit */}
              <Button
                onClick={handleSubmit}
                disabled={rating === 0 || submitting}
                size="sm"
                className="w-full bg-[#FF9100] hover:bg-[#FF9100]/80 text-white text-[10px]"
              >
                {submitting ? "저장 중..." : "평가 저장"}
              </Button>
            </div>
          </motion.div>
        )}
      </AnimatePresence>
    </div>
  );
}

// ─── Feedback History Component ───

interface FeedbackHistoryProps {
  category?: string;
  limit?: number;
}

export function FeedbackHistory({ category, limit = 10 }: FeedbackHistoryProps) {
  const [entries, setEntries] = useState<FeedbackEntry[]>([]);
  const [loading, setLoading] = useState(true);

  useEffect(() => {
    async function load() {
      let query = supabase
        .from("researcher_feedback")
        .select("*")
        .order("created_at", { ascending: false })
        .limit(limit);

      if (category) {
        query = query.eq("category", category);
      }

      const { data } = await query;
      setEntries((data as FeedbackEntry[]) || []);
      setLoading(false);
    }
    load();
  }, [category, limit]);

  if (loading) return null;
  if (entries.length === 0) return null;

  return (
    <div className="space-y-2">
      <div className="flex items-center gap-2">
        <span className="text-xs font-medium text-white">최근 연구원 피드백</span>
        <Badge variant="outline" className="text-[9px] border-[#FF9100]/30 text-[#FF9100]">
          {entries.length}건
        </Badge>
      </div>
      <div className="space-y-1.5">
        {entries.map((e) => (
          <div key={e.id} className="p-2.5 rounded-lg bg-[#112240]/30 border border-white/[0.04]">
            <div className="flex items-center justify-between mb-1">
              <div className="flex items-center gap-2">
                <span className="text-[10px] font-medium text-white">{e.material_name || "Unknown"}</span>
                <span className="text-[9px] text-[#FFD600]">{"★".repeat(e.rating)}{"☆".repeat(5 - e.rating)}</span>
              </div>
              <Badge variant="outline" className={`text-[8px] ${
                e.experiment_success ? "border-[#00E676]/20 text-[#00E676]" : "border-[#FF5252]/20 text-[#FF5252]"
              }`}>
                {e.experiment_success ? "성공" : "실패/미실험"}
              </Badge>
            </div>
            {e.notes && <p className="text-[9px] text-[#8892B0] leading-relaxed">{e.notes}</p>}
            {(e.pros || e.cons) && (
              <div className="flex gap-2 mt-1">
                {e.pros && <span className="text-[8px] text-[#00E676]">✓ {e.pros}</span>}
                {e.cons && <span className="text-[8px] text-[#FF5252]">✗ {e.cons}</span>}
              </div>
            )}
            <div className="flex justify-between mt-1">
              <span className="text-[8px] text-[#8892B0]">{e.researcher || "익명"}</span>
              <span className="text-[8px] text-[#8892B0]">{new Date(e.created_at).toLocaleDateString("ko")}</span>
            </div>
          </div>
        ))}
      </div>
    </div>
  );
}

// ─── Helper: Get feedback summary for AI prompts ───

export async function getFeedbackForAI(category?: string): Promise<string> {
  let query = supabase
    .from("researcher_feedback")
    .select("*")
    .order("created_at", { ascending: false })
    .limit(20);

  if (category) {
    query = query.eq("category", category);
  }

  const { data } = await query;
  if (!data || data.length === 0) return "";

  const lines = (data as FeedbackEntry[]).map((e) => {
    const parts = [`${e.material_name} (${e.material_smiles || "N/A"}): 평점 ${e.rating}/5`];
    if (e.experiment_success) parts.push("실험 성공");
    if (e.actual_thermal_stability) parts.push(`실측 열안정성: ${e.actual_thermal_stability}°C`);
    if (e.actual_bandgap) parts.push(`실측 밴드갭: ${e.actual_bandgap}eV`);
    if (e.pros) parts.push(`장점: ${e.pros}`);
    if (e.cons) parts.push(`단점: ${e.cons}`);
    if (e.notes) parts.push(`메모: ${e.notes}`);
    return parts.join(" | ");
  });

  return `\n\n## 연구원 실험 피드백 (최근 ${data.length}건) — 이 정보를 반드시 참고하세요\n${lines.join("\n")}`;
}
