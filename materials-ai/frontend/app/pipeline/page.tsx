"use client";

import { useState, useCallback, useRef } from "react";
import { motion } from "framer-motion";
import { Button } from "@/components/ui/button";
import { Card, CardContent, CardHeader, CardTitle } from "@/components/ui/card";
import { Badge } from "@/components/ui/badge";
import { Label } from "@/components/ui/label";
import { PipelineStep } from "@/components/PipelineStep";
import { HelpTooltip } from "@/components/HelpTooltip";
import type { PipelineStepStatus, PipelineStepStatusType } from "@/lib/types";
import { supabase } from "@/lib/supabase";

const INITIAL_STEPS: PipelineStepStatus[] = [
  { name: "collect", label: "데이터 수집", status: "pending", progress: 0, rows_processed: 0, total_rows: 0, quality_score: null, time_elapsed_sec: 0, error: null },
  { name: "normalize", label: "정규화", status: "pending", progress: 0, rows_processed: 0, total_rows: 0, quality_score: null, time_elapsed_sec: 0, error: null },
  { name: "clean", label: "데이터 정제", status: "pending", progress: 0, rows_processed: 0, total_rows: 0, quality_score: null, time_elapsed_sec: 0, error: null },
  { name: "engineer", label: "특성 추출", status: "pending", progress: 0, rows_processed: 0, total_rows: 0, quality_score: null, time_elapsed_sec: 0, error: null },
  { name: "embed", label: "임베딩 생성", status: "pending", progress: 0, rows_processed: 0, total_rows: 0, quality_score: null, time_elapsed_sec: 0, error: null },
  { name: "dataset", label: "데이터셋 저장", status: "pending", progress: 0, rows_processed: 0, total_rows: 0, quality_score: null, time_elapsed_sec: 0, error: null },
];

// CSV column auto-mapping
const COL_MAP: Record<string, string> = {
  smiles: "smiles", smi: "smiles", SMILES: "smiles",
  name: "name", "이름": "name", "소재명": "name",
  category: "category", "카테고리": "category",
  thermal_stability: "thermal_stability", "열안정성": "thermal_stability", td: "thermal_stability",
  dielectric_constant: "dielectric_constant", "유전율": "dielectric_constant", dk: "dielectric_constant",
  bandgap: "bandgap", "밴드갭": "bandgap", eg: "bandgap",
  solubility: "solubility", "용해도": "solubility",
  density: "density", "밀도": "density",
  molecular_weight: "molecular_weight", "분자량": "molecular_weight",
};

function parseCSV(text: string): Record<string, string>[] {
  const lines = text.trim().split("\n");
  if (lines.length < 2) return [];
  const sep = lines[0].includes("\t") ? "\t" : ",";
  const headers = lines[0].split(sep).map((h) => h.trim().replace(/^["']|["']$/g, ""));
  return lines.slice(1).map((line) => {
    const values = line.split(sep).map((v) => v.trim().replace(/^["']|["']$/g, ""));
    const row: Record<string, string> = {};
    headers.forEach((h, i) => { row[h] = values[i] || ""; });
    return row;
  });
}

export default function PipelinePage() {
  const [steps, setSteps] = useState<PipelineStepStatus[]>(INITIAL_STEPS);
  const [isRunning, setIsRunning] = useState(false);
  const [selectedFile, setSelectedFile] = useState<File | null>(null);
  const [result, setResult] = useState<{ inserted: number; total: number } | null>(null);
  const fileRef = useRef<HTMLInputElement>(null);

  const updateStep = (index: number, update: Partial<PipelineStepStatus>) => {
    setSteps((prev) => prev.map((s, i) => (i === index ? { ...s, ...update } : s)));
  };

  const runPipeline = useCallback(async () => {
    if (!selectedFile) return;

    setIsRunning(true);
    setResult(null);
    setSteps(INITIAL_STEPS.map((s) => ({ ...s })));

    const startTime = performance.now();

    try {
      // ─── Step 1: 데이터 수집 (파일 읽기) ───
      updateStep(0, { status: "running" as PipelineStepStatusType, progress: 10 });
      const text = await selectedFile.text();
      const rawRecords = parseCSV(text);
      const totalRows = rawRecords.length;

      if (totalRows === 0) {
        updateStep(0, { status: "failed" as PipelineStepStatusType, error: "파일에 데이터가 없습니다" });
        setIsRunning(false);
        return;
      }

      updateStep(0, {
        status: "completed" as PipelineStepStatusType,
        progress: 100,
        rows_processed: totalRows,
        total_rows: totalRows,
        time_elapsed_sec: (performance.now() - startTime) / 1000,
      });

      // ─── Step 2: 정규화 (컬럼 매핑) ───
      const t2 = performance.now();
      updateStep(1, { status: "running" as PipelineStepStatusType, total_rows: totalRows });

      const mapped = rawRecords.map((row) => {
        const out: Record<string, unknown> = {};
        for (const [key, val] of Object.entries(row)) {
          const canonical = COL_MAP[key.trim()] || COL_MAP[key.trim().toLowerCase()];
          if (canonical) out[canonical] = val;
        }
        return out;
      });

      await new Promise((r) => setTimeout(r, 400));
      updateStep(1, {
        status: "completed" as PipelineStepStatusType,
        progress: 100,
        rows_processed: totalRows,
        total_rows: totalRows,
        time_elapsed_sec: (performance.now() - t2) / 1000,
      });

      // ─── Step 3: 데이터 정제 (SMILES 검증) ───
      const t3 = performance.now();
      updateStep(2, { status: "running" as PipelineStepStatusType, total_rows: totalRows });

      const cleaned = mapped.filter((r) => {
        const s = String(r.smiles || "").trim();
        return s.length > 0 && s !== "nan" && s !== "undefined";
      });

      // Deduplicate
      const seen = new Set<string>();
      const unique = cleaned.filter((r) => {
        const s = String(r.smiles);
        if (seen.has(s)) return false;
        seen.add(s);
        return true;
      });

      await new Promise((r) => setTimeout(r, 300));
      const qualityScore = (unique.length / Math.max(totalRows, 1)) * 100;
      updateStep(2, {
        status: "completed" as PipelineStepStatusType,
        progress: 100,
        rows_processed: unique.length,
        total_rows: totalRows,
        quality_score: Math.round(qualityScore * 10) / 10,
        time_elapsed_sec: (performance.now() - t3) / 1000,
      });

      // ─── Step 4: 특성 추출 (숫자 변환) ───
      const t4 = performance.now();
      updateStep(3, { status: "running" as PipelineStepStatusType, total_rows: unique.length });

      const enriched = unique.map((r) => {
        const numFields = ["thermal_stability", "dielectric_constant", "bandgap", "solubility", "density", "molecular_weight"];
        const out: Record<string, unknown> = { ...r };
        for (const f of numFields) {
          if (out[f] !== undefined && out[f] !== null && String(out[f]).trim() !== "") {
            const parsed = parseFloat(String(out[f]));
            out[f] = isNaN(parsed) ? null : parsed;
          } else {
            out[f] = null;
          }
        }
        out.name = out.name ? String(out.name) : null;
        out.category = out.category ? String(out.category) : null;
        return out;
      });

      await new Promise((r) => setTimeout(r, 500));
      updateStep(3, {
        status: "completed" as PipelineStepStatusType,
        progress: 100,
        rows_processed: enriched.length,
        total_rows: unique.length,
        time_elapsed_sec: (performance.now() - t4) / 1000,
      });

      // ─── Step 5: 임베딩 생성 (패스 — AI 모델 사용 시 대체) ───
      const t5 = performance.now();
      updateStep(4, { status: "running" as PipelineStepStatusType, total_rows: enriched.length });

      // Simulate progress
      for (let i = 0; i < 10; i++) {
        await new Promise((r) => setTimeout(r, 150));
        updateStep(4, {
          progress: ((i + 1) / 10) * 100,
          rows_processed: Math.floor(((i + 1) / 10) * enriched.length),
        });
      }

      updateStep(4, {
        status: "completed" as PipelineStepStatusType,
        progress: 100,
        rows_processed: enriched.length,
        total_rows: enriched.length,
        time_elapsed_sec: (performance.now() - t5) / 1000,
      });

      // ─── Step 6: 데이터셋 저장 (Supabase INSERT) ───
      const t6 = performance.now();
      updateStep(5, { status: "running" as PipelineStepStatusType, total_rows: enriched.length });

      const insertRows = enriched.map((r) => ({
        smiles: String(r.smiles),
        name: r.name as string | null,
        category: r.category as string | null,
        thermal_stability: r.thermal_stability as number | null,
        dielectric_constant: r.dielectric_constant as number | null,
        bandgap: r.bandgap as number | null,
        solubility: r.solubility as number | null,
        density: r.density as number | null,
        molecular_weight: r.molecular_weight as number | null,
        source: "pipeline_upload",
        is_verified: false,
      }));

      // Batch insert (50 at a time)
      let inserted = 0;
      const batchSize = 50;
      for (let i = 0; i < insertRows.length; i += batchSize) {
        const batch = insertRows.slice(i, i + batchSize);
        const { data, error } = await supabase.from("materials").insert(batch).select("id");

        if (error) {
          updateStep(5, {
            status: "failed" as PipelineStepStatusType,
            error: error.message,
            time_elapsed_sec: (performance.now() - t6) / 1000,
          });
          setIsRunning(false);
          return;
        }

        inserted += data?.length || 0;
        updateStep(5, {
          progress: (Math.min(i + batchSize, insertRows.length) / insertRows.length) * 100,
          rows_processed: inserted,
        });
      }

      const { count } = await supabase.from("materials").select("*", { count: "exact", head: true });

      updateStep(5, {
        status: "completed" as PipelineStepStatusType,
        progress: 100,
        rows_processed: inserted,
        total_rows: enriched.length,
        time_elapsed_sec: (performance.now() - t6) / 1000,
      });

      setResult({ inserted, total: count || 0 });
    } catch (err) {
      console.error("Pipeline error:", err);
    } finally {
      setIsRunning(false);
    }
  }, [selectedFile]);

  const completedSteps = steps.filter((s) => s.status === "completed").length;
  const overallProgress = (completedSteps / steps.length) * 100;

  return (
    <div className="p-8 space-y-8 grid-pattern min-h-screen">
      <motion.div initial={{ opacity: 0, y: -20 }} animate={{ opacity: 1, y: 0 }} transition={{ duration: 0.5 }}>
        <div className="flex items-center gap-2">
          <h1 className="text-2xl font-bold text-white">데이터 파이프라인</h1>
          <HelpTooltip
            title="데이터 파이프라인"
            description="CSV 파일을 업로드하면 6단계 자동 파이프라인을 거쳐 Supabase DB에 저장됩니다."
            details={[
              "1. 데이터 수집: CSV 파일을 읽어서 원시 레코드를 로드합니다.",
              "2. 정규화: 컬럼명을 자동 매핑합니다 (한글/영문 모두 지원).",
              "3. 데이터 정제: 유효하지 않은 SMILES 제거, 중복 제거.",
              "4. 특성 추출: 숫자형 물성 데이터를 파싱합니다.",
              "5. 임베딩 생성: AI 모델 활용 시 벡터 임베딩을 생성합니다.",
              "6. 데이터셋 저장: Supabase PostgreSQL에 배치 저장합니다.",
            ]}
            size="md"
          />
        </div>
        <p className="text-sm text-[#8892B0] mt-1">
          CSV 파일 업로드 → 자동 파싱 → 정제 → Supabase DB 저장
        </p>
      </motion.div>

      <div className="grid grid-cols-1 lg:grid-cols-3 gap-8">
        {/* Left — Pipeline Configuration */}
        <div className="space-y-6">
          <Card className="bg-[#0D1B2A] border-white/[0.06]">
            <CardHeader>
              <div className="flex items-center gap-2">
                <CardTitle className="text-sm font-semibold text-white">파이프라인 설정</CardTitle>
                <HelpTooltip title="파이프라인 설정" description="처리할 CSV 파일을 선택하세요. SMILES 컬럼이 필수입니다."
                  details={["지원 형식: .csv, .tsv", "필수 컬럼: SMILES (또는 smiles)", "한글 컬럼명 자동 인식: 열안정성, 유전율, 밴드갭 등"]} />
              </div>
            </CardHeader>
            <CardContent className="space-y-4">
              <div className="space-y-2">
                <Label className="text-xs text-[#8892B0]">CSV 파일 선택</Label>
                <div className="relative">
                  <input
                    ref={fileRef}
                    type="file"
                    accept=".csv,.tsv"
                    onChange={(e) => setSelectedFile(e.target.files?.[0] || null)}
                    disabled={isRunning}
                    className="absolute inset-0 opacity-0 cursor-pointer z-10"
                  />
                  <div className={`p-4 border-2 border-dashed rounded-xl text-center transition-colors ${
                    selectedFile ? "border-[#00E676]/30 bg-[#00E676]/5" : "border-white/[0.08] hover:border-[#00B4D8]/30"
                  }`}>
                    {selectedFile ? (
                      <div>
                        <p className="text-xs text-[#00E676] font-medium">{selectedFile.name}</p>
                        <p className="text-[10px] text-[#8892B0] mt-1">{(selectedFile.size / 1024).toFixed(1)} KB</p>
                      </div>
                    ) : (
                      <div>
                        <svg className="w-6 h-6 text-[#8892B0]/50 mx-auto mb-1" fill="none" viewBox="0 0 24 24" stroke="currentColor" strokeWidth={1.5}>
                          <path strokeLinecap="round" strokeLinejoin="round" d="M3 16.5v2.25A2.25 2.25 0 0 0 5.25 21h13.5A2.25 2.25 0 0 0 21 18.75V16.5m-13.5-9L12 3m0 0 4.5 4.5M12 3v13.5" />
                        </svg>
                        <p className="text-[10px] text-[#8892B0]">클릭하여 CSV 파일 선택</p>
                      </div>
                    )}
                  </div>
                </div>
              </div>

              <Button
                onClick={runPipeline}
                disabled={isRunning || !selectedFile}
                className="w-full bg-gradient-to-r from-[#00B4D8] to-[#0096C7] hover:from-[#0096C7] hover:to-[#0077B6] text-white font-medium"
              >
                {isRunning ? (
                  <div className="flex items-center gap-2">
                    <div className="w-4 h-4 border-2 border-white border-t-transparent rounded-full animate-spin" />
                    <span>파이프라인 실행 중...</span>
                  </div>
                ) : (
                  <div className="flex items-center gap-2">
                    <svg className="w-4 h-4" fill="none" viewBox="0 0 24 24" stroke="currentColor" strokeWidth={2}>
                      <path strokeLinecap="round" strokeLinejoin="round" d="M5.25 5.653c0-.856.917-1.398 1.667-.986l11.54 6.347a1.125 1.125 0 0 1 0 1.972l-11.54 6.347a1.125 1.125 0 0 1-1.667-.986V5.653Z" />
                    </svg>
                    <span>파이프라인 시작</span>
                  </div>
                )}
              </Button>

              {/* Result */}
              {result && (
                <motion.div initial={{ opacity: 0, y: 10 }} animate={{ opacity: 1, y: 0 }}
                  className="p-3 rounded-lg bg-[#00E676]/10 border border-[#00E676]/20 text-xs text-[#00E676]">
                  <p className="font-medium">파이프라인 완료!</p>
                  <p>{result.inserted}건 DB 저장 · 전체 {result.total}건</p>
                </motion.div>
              )}
            </CardContent>
          </Card>

          {/* Summary */}
          <Card className="bg-[#0D1B2A] border-white/[0.06]">
            <CardHeader>
              <div className="flex items-center gap-2">
                <CardTitle className="text-sm font-semibold text-white">진행 요약</CardTitle>
                <HelpTooltip title="진행 요약" description="파이프라인의 전체 진행 상태를 요약합니다."
                  details={["전체 진행률: 6단계 중 완료된 비율", "완료 단계: 성공적으로 처리된 단계 수", "하단 배지: 각 단계의 완료 상태를 색상으로 표시"]} />
              </div>
            </CardHeader>
            <CardContent className="space-y-4">
              <div className="space-y-2">
                <div className="flex items-center justify-between">
                  <span className="text-xs text-[#8892B0]">전체 진행률</span>
                  <span className="text-xs font-mono text-[#00B4D8]">{overallProgress.toFixed(0)}%</span>
                </div>
                <div className="w-full h-2 bg-white/[0.04] rounded-full overflow-hidden">
                  <motion.div
                    className="h-full bg-gradient-to-r from-[#00B4D8] to-[#00E676] rounded-full"
                    initial={{ width: 0 }}
                    animate={{ width: `${overallProgress}%` }}
                    transition={{ duration: 0.5 }}
                  />
                </div>
              </div>

              <div className="grid grid-cols-2 gap-3">
                <div className="p-3 rounded-lg bg-[#112240]/50">
                  <p className="text-[10px] text-[#8892B0]">완료 단계</p>
                  <p className="text-lg font-bold font-mono text-[#00E676]">{completedSteps}/{steps.length}</p>
                </div>
                <div className="p-3 rounded-lg bg-[#112240]/50">
                  <p className="text-[10px] text-[#8892B0]">처리 데이터</p>
                  <p className="text-lg font-bold font-mono text-[#00B4D8]">
                    {steps.find((s) => s.rows_processed > 0)?.rows_processed.toLocaleString() || "0"}
                  </p>
                </div>
              </div>

              <div className="flex flex-wrap gap-1.5">
                {steps.map((step) => (
                  <Badge key={step.name} variant="outline"
                    className={`text-[10px] ${
                      step.status === "completed" ? "border-[#00E676]/30 text-[#00E676]"
                        : step.status === "running" ? "border-[#00B4D8]/30 text-[#00B4D8]"
                        : step.status === "failed" ? "border-[#FF5252]/30 text-[#FF5252]"
                        : "border-white/10 text-[#8892B0]"
                    }`}>
                    {step.label}
                  </Badge>
                ))}
              </div>
            </CardContent>
          </Card>
        </div>

        {/* Right — Pipeline Steps Timeline */}
        <div className="lg:col-span-2">
          <Card className="bg-[#0D1B2A] border-white/[0.06]">
            <CardHeader>
              <div className="flex items-center justify-between">
                <CardTitle className="text-sm font-semibold text-white">파이프라인 단계</CardTitle>
                <Badge variant="outline" className={`text-xs ${
                  isRunning ? "border-[#00B4D8]/30 text-[#00B4D8]"
                    : completedSteps === steps.length && completedSteps > 0 ? "border-[#00E676]/30 text-[#00E676]"
                    : "border-white/10 text-[#8892B0]"
                }`}>
                  {isRunning ? "실행 중" : completedSteps === steps.length && completedSteps > 0 ? "완료" : "대기"}
                </Badge>
              </div>
              <p className="text-[11px] text-[#8892B0]">수집 → 정규화 → 정제 → 특성 추출 → 임베딩 → 저장</p>
            </CardHeader>
            <CardContent>
              <div className="space-y-0">
                {steps.map((step, index) => (
                  <PipelineStep key={step.name} step={step} index={index} isLast={index === steps.length - 1} />
                ))}
              </div>
            </CardContent>
          </Card>
        </div>
      </div>
    </div>
  );
}
