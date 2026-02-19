"use client";

import { useState, useCallback } from "react";
import { motion } from "framer-motion";
import { Button } from "@/components/ui/button";
import { Card, CardContent, CardHeader, CardTitle } from "@/components/ui/card";
import { Badge } from "@/components/ui/badge";
import {
  Select,
  SelectContent,
  SelectItem,
  SelectTrigger,
  SelectValue,
} from "@/components/ui/select";
import { Input } from "@/components/ui/input";
import { Label } from "@/components/ui/label";
import { PipelineStep } from "@/components/PipelineStep";
import { HelpTooltip } from "@/components/HelpTooltip";
import type { PipelineStepStatus, PipelineStepStatusType } from "@/lib/types";

const INITIAL_STEPS: PipelineStepStatus[] = [
  { name: "collect", label: "데이터 수집", status: "pending", progress: 0, rows_processed: 0, total_rows: 0, quality_score: null, time_elapsed_sec: 0, error: null },
  { name: "normalize", label: "정규화", status: "pending", progress: 0, rows_processed: 0, total_rows: 0, quality_score: null, time_elapsed_sec: 0, error: null },
  { name: "clean", label: "데이터 정제", status: "pending", progress: 0, rows_processed: 0, total_rows: 0, quality_score: null, time_elapsed_sec: 0, error: null },
  { name: "engineer", label: "특성 추출", status: "pending", progress: 0, rows_processed: 0, total_rows: 0, quality_score: null, time_elapsed_sec: 0, error: null },
  { name: "embed", label: "임베딩 생성", status: "pending", progress: 0, rows_processed: 0, total_rows: 0, quality_score: null, time_elapsed_sec: 0, error: null },
  { name: "dataset", label: "데이터셋 생성", status: "pending", progress: 0, rows_processed: 0, total_rows: 0, quality_score: null, time_elapsed_sec: 0, error: null },
];

export default function PipelinePage() {
  const [steps, setSteps] = useState<PipelineStepStatus[]>(INITIAL_STEPS);
  const [isRunning, setIsRunning] = useState(false);
  const [sourceType, setSourceType] = useState("csv");
  const [filePath, setFilePath] = useState("");

  const simulatePipeline = useCallback(async () => {
    setIsRunning(true);
    setSteps(INITIAL_STEPS.map((s) => ({ ...s })));

    const totalRows = Math.floor(Math.random() * 5000) + 1000;

    for (let stepIndex = 0; stepIndex < INITIAL_STEPS.length; stepIndex++) {
      const stepName = INITIAL_STEPS[stepIndex].name;
      const stepDuration = 1500 + Math.random() * 2000;
      const intervals = 20;
      const intervalTime = stepDuration / intervals;

      // Set step to running
      setSteps((prev) =>
        prev.map((s, i) =>
          i === stepIndex
            ? { ...s, status: "running" as PipelineStepStatusType, total_rows: totalRows }
            : s
        )
      );

      // Simulate progress
      for (let p = 1; p <= intervals; p++) {
        await new Promise((r) => setTimeout(r, intervalTime));
        const progress = (p / intervals) * 100;
        const rowsProcessed = Math.floor((p / intervals) * totalRows);

        setSteps((prev) =>
          prev.map((s, i) =>
            i === stepIndex
              ? {
                  ...s,
                  progress,
                  rows_processed: rowsProcessed,
                  time_elapsed_sec: (p * intervalTime) / 1000,
                }
              : s
          )
        );
      }

      // Complete step
      setSteps((prev) =>
        prev.map((s, i) =>
          i === stepIndex
            ? {
                ...s,
                status: "completed" as PipelineStepStatusType,
                progress: 100,
                rows_processed: totalRows,
                quality_score:
                  stepName === "clean"
                    ? 85 + Math.random() * 12
                    : null,
                time_elapsed_sec: stepDuration / 1000,
              }
            : s
        )
      );
    }

    setIsRunning(false);
  }, []);

  const handleStart = useCallback(async () => {
    // Try real API first
    try {
      const apiUrl = process.env.NEXT_PUBLIC_API_URL || "http://localhost:8000";
      const res = await fetch(`${apiUrl}/api/pipeline/start`, {
        method: "POST",
        headers: { "Content-Type": "application/json" },
        body: JSON.stringify({
          source_type: sourceType,
          file_path: filePath || "/data/sample.csv",
          options: {},
        }),
      });

      if (res.ok) {
        // TODO: Connect to WebSocket for real-time updates
        return;
      }
    } catch {
      // API not available, simulate
    }

    await simulatePipeline();
  }, [sourceType, filePath, simulatePipeline]);

  const completedSteps = steps.filter((s) => s.status === "completed").length;
  const overallProgress = (completedSteps / steps.length) * 100;

  return (
    <div className="p-8 space-y-8 grid-pattern min-h-screen">
      {/* Header */}
      <motion.div
        initial={{ opacity: 0, y: -20 }}
        animate={{ opacity: 1, y: 0 }}
        transition={{ duration: 0.5 }}
      >
        <div className="flex items-center gap-2">
          <h1 className="text-2xl font-bold text-white">데이터 파이프라인</h1>
          <HelpTooltip
            title="데이터 파이프라인"
            description="원시 데이터(CSV, SDF 등)를 AI가 학습할 수 있는 형태로 변환하는 6단계 자동화 파이프라인입니다."
            details={[
              "1. 데이터 수집: 파일을 읽어서 원시 레코드를 로드합니다.",
              "2. 정규화: 단위 통일, SMILES 표준화 등을 수행합니다.",
              "3. 데이터 정제: 중복 제거, 결측치 처리를 합니다.",
              "4. 특성 추출: RDKit로 분자 기술자를 계산합니다.",
              "5. 임베딩 생성: 분자를 768차원 벡터로 변환합니다.",
              "6. 데이터셋 생성: 최종 학습 데이터셋을 DB에 저장합니다.",
            ]}
            size="md"
          />
        </div>
        <p className="text-sm text-[#8892B0] mt-1">
          원시 데이터를 AI 학습용 데이터셋으로 변환하는 파이프라인 관리
        </p>
      </motion.div>

      <div className="grid grid-cols-1 lg:grid-cols-3 gap-8">
        {/* Left — Pipeline Configuration */}
        <div className="space-y-6">
          <Card className="bg-[#0D1B2A] border-white/[0.06]">
            <CardHeader>
              <div className="flex items-center gap-2">
                <CardTitle className="text-sm font-semibold text-white">
                  파이프라인 설정
                </CardTitle>
                <HelpTooltip
                  title="파이프라인 설정"
                  description="처리할 데이터의 형식과 파일 경로를 지정합니다."
                  details={[
                    "CSV: 쉼표 또는 탭으로 구분된 표 형식 데이터",
                    "SDF: 분자 구조 데이터 파일 (3D 좌표 포함)",
                    "MOL: 단일 분자 구조 파일",
                    "JSON: 키-값 형식의 데이터",
                  ]}
                />
              </div>
            </CardHeader>
            <CardContent className="space-y-4">
              <div className="space-y-2">
                <Label className="text-xs text-[#8892B0]">데이터 형식</Label>
                <Select value={sourceType} onValueChange={setSourceType}>
                  <SelectTrigger className="bg-[#112240] border-white/[0.08] text-white">
                    <SelectValue />
                  </SelectTrigger>
                  <SelectContent className="bg-[#112240] border-white/[0.08]">
                    <SelectItem value="csv" className="text-white">CSV</SelectItem>
                    <SelectItem value="sdf" className="text-white">SDF</SelectItem>
                    <SelectItem value="mol" className="text-white">MOL</SelectItem>
                    <SelectItem value="json" className="text-white">JSON</SelectItem>
                  </SelectContent>
                </Select>
              </div>

              <div className="space-y-2">
                <Label className="text-xs text-[#8892B0]">파일 경로</Label>
                <Input
                  value={filePath}
                  onChange={(e) => setFilePath(e.target.value)}
                  placeholder="/data/materials.csv"
                  className="bg-[#112240] border-white/[0.08] text-white placeholder:text-[#8892B0]/50 font-mono text-sm"
                />
              </div>

              <Button
                onClick={handleStart}
                disabled={isRunning}
                className="w-full bg-gradient-to-r from-[#00B4D8] to-[#0096C7] hover:from-[#0096C7] hover:to-[#0077B6] text-white font-medium"
              >
                {isRunning ? (
                  <div className="flex items-center gap-2">
                    <div className="w-4 h-4 border-2 border-white border-t-transparent rounded-full animate-spin" />
                    <span>파이프라인 실행 중...</span>
                  </div>
                ) : (
                  <div className="flex items-center gap-2">
                    <svg
                      className="w-4 h-4"
                      fill="none"
                      viewBox="0 0 24 24"
                      stroke="currentColor"
                      strokeWidth={2}
                    >
                      <path
                        strokeLinecap="round"
                        strokeLinejoin="round"
                        d="M5.25 5.653c0-.856.917-1.398 1.667-.986l11.54 6.347a1.125 1.125 0 0 1 0 1.972l-11.54 6.347a1.125 1.125 0 0 1-1.667-.986V5.653Z"
                      />
                    </svg>
                    <span>파이프라인 시작</span>
                  </div>
                )}
              </Button>
            </CardContent>
          </Card>

          {/* Summary */}
          <Card className="bg-[#0D1B2A] border-white/[0.06]">
            <CardHeader>
              <div className="flex items-center gap-2">
                <CardTitle className="text-sm font-semibold text-white">
                  진행 요약
                </CardTitle>
                <HelpTooltip
                  title="진행 요약"
                  description="파이프라인의 전체 진행 상태를 요약합니다."
                  details={[
                    "전체 진행률: 6단계 중 완료된 비율",
                    "완료 단계: 성공적으로 처리된 단계 수",
                    "처리 데이터: 현재까지 처리된 행(row) 수",
                    "하단 배지: 각 단계의 완료 상태를 색상으로 표시",
                  ]}
                />
              </div>
            </CardHeader>
            <CardContent className="space-y-4">
              <div className="space-y-2">
                <div className="flex items-center justify-between">
                  <span className="text-xs text-[#8892B0]">전체 진행률</span>
                  <span className="text-xs font-mono text-[#00B4D8]">
                    {overallProgress.toFixed(0)}%
                  </span>
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
                  <p className="text-lg font-bold font-mono text-[#00E676]">
                    {completedSteps}/{steps.length}
                  </p>
                </div>
                <div className="p-3 rounded-lg bg-[#112240]/50">
                  <p className="text-[10px] text-[#8892B0]">처리 데이터</p>
                  <p className="text-lg font-bold font-mono text-[#00B4D8]">
                    {steps
                      .find((s) => s.status === "running" || s.status === "completed")
                      ?.rows_processed.toLocaleString() || "0"}
                  </p>
                </div>
              </div>

              {/* Step badges */}
              <div className="flex flex-wrap gap-1.5">
                {steps.map((step) => (
                  <Badge
                    key={step.name}
                    variant="outline"
                    className={`text-[10px] ${
                      step.status === "completed"
                        ? "border-[#00E676]/30 text-[#00E676]"
                        : step.status === "running"
                          ? "border-[#00B4D8]/30 text-[#00B4D8]"
                          : step.status === "failed"
                            ? "border-[#FF5252]/30 text-[#FF5252]"
                            : "border-white/10 text-[#8892B0]"
                    }`}
                  >
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
                <CardTitle className="text-sm font-semibold text-white">
                  파이프라인 단계
                </CardTitle>
                <Badge
                  variant="outline"
                  className={`text-xs ${
                    isRunning
                      ? "border-[#00B4D8]/30 text-[#00B4D8]"
                      : completedSteps === steps.length && completedSteps > 0
                        ? "border-[#00E676]/30 text-[#00E676]"
                        : "border-white/10 text-[#8892B0]"
                  }`}
                >
                  {isRunning
                    ? "실행 중"
                    : completedSteps === steps.length && completedSteps > 0
                      ? "완료"
                      : "대기"}
                </Badge>
              </div>
              <p className="text-[11px] text-[#8892B0]">
                수집 → 정규화 → 정제 → 특성 추출 → 임베딩 → 데이터셋
              </p>
            </CardHeader>
            <CardContent>
              <div className="space-y-0">
                {steps.map((step, index) => (
                  <PipelineStep
                    key={step.name}
                    step={step}
                    index={index}
                    isLast={index === steps.length - 1}
                  />
                ))}
              </div>
            </CardContent>
          </Card>
        </div>
      </div>
    </div>
  );
}
