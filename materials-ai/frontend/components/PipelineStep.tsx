"use client";

import { Progress } from "@/components/ui/progress";
import type { PipelineStepStatus } from "@/lib/types";

interface PipelineStepProps {
  step: PipelineStepStatus;
  index: number;
  isLast: boolean;
}

const STATUS_CONFIG = {
  pending: {
    color: "#8892B0",
    bgColor: "bg-[#8892B0]/10",
    borderColor: "border-[#8892B0]/20",
    icon: "○",
    label: "대기",
  },
  running: {
    color: "#00B4D8",
    bgColor: "bg-[#00B4D8]/10",
    borderColor: "border-[#00B4D8]/20",
    icon: "◉",
    label: "진행",
  },
  completed: {
    color: "#00E676",
    bgColor: "bg-[#00E676]/10",
    borderColor: "border-[#00E676]/20",
    icon: "✓",
    label: "완료",
  },
  failed: {
    color: "#FF5252",
    bgColor: "bg-[#FF5252]/10",
    borderColor: "border-[#FF5252]/20",
    icon: "✗",
    label: "실패",
  },
};

export function PipelineStep({ step, index, isLast }: PipelineStepProps) {
  const config = STATUS_CONFIG[step.status] || STATUS_CONFIG.pending;

  return (
    <div className="flex gap-4">
      {/* Timeline */}
      <div className="flex flex-col items-center">
        <div
          className={`w-10 h-10 rounded-full flex items-center justify-center text-sm font-bold border-2 transition-all duration-500 ${config.bgColor} ${config.borderColor}`}
          style={{ color: config.color }}
        >
          {step.status === "running" ? (
            <div
              className="w-5 h-5 border-2 border-t-transparent rounded-full animate-spin"
              style={{ borderColor: config.color, borderTopColor: "transparent" }}
            />
          ) : (
            config.icon
          )}
        </div>
        {!isLast && (
          <div
            className="w-0.5 flex-1 min-h-[40px] transition-colors duration-500"
            style={{
              backgroundColor:
                step.status === "completed"
                  ? "#00E676"
                  : "rgba(255,255,255,0.06)",
            }}
          />
        )}
      </div>

      {/* Content */}
      <div className="flex-1 pb-8">
        <div
          className={`p-4 rounded-xl border transition-all duration-500 ${config.bgColor} ${config.borderColor}`}
        >
          <div className="flex items-center justify-between mb-2">
            <div>
              <h3 className="text-sm font-semibold text-white">
                {step.label}
              </h3>
              <p className="text-[11px] text-[#8892B0]">
                Step {index + 1} · {config.label}
              </p>
            </div>
            {step.time_elapsed_sec > 0 && (
              <span className="text-xs font-mono text-[#8892B0]">
                {step.time_elapsed_sec.toFixed(1)}s
              </span>
            )}
          </div>

          {/* Progress bar */}
          {(step.status === "running" || step.status === "completed") && (
            <div className="space-y-2">
              <Progress
                value={step.progress}
                className="h-1.5 bg-white/[0.04]"
              />
              <div className="flex items-center justify-between">
                <span className="text-[10px] text-[#8892B0]">
                  {step.rows_processed.toLocaleString()} /{" "}
                  {step.total_rows.toLocaleString()} 행
                </span>
                <span
                  className="text-[10px] font-mono font-medium"
                  style={{ color: config.color }}
                >
                  {step.progress.toFixed(1)}%
                </span>
              </div>
            </div>
          )}

          {/* Quality Score */}
          {step.quality_score !== null && step.quality_score !== undefined && (
            <div className="mt-2 flex items-center gap-2">
              <span className="text-[10px] text-[#8892B0]">품질 점수</span>
              <span className="text-xs font-mono text-[#00E676]">
                {step.quality_score.toFixed(1)}%
              </span>
            </div>
          )}

          {/* Error message */}
          {step.error && (
            <div className="mt-2 p-2 rounded bg-[#FF5252]/10 border border-[#FF5252]/20">
              <p className="text-[11px] text-[#FF5252]">{step.error}</p>
            </div>
          )}
        </div>
      </div>
    </div>
  );
}
