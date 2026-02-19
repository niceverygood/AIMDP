"use client";

import { useState } from "react";
import {
  Popover,
  PopoverContent,
  PopoverTrigger,
} from "@/components/ui/popover";

interface HelpTooltipProps {
  /** Short title for the help popup */
  title: string;
  /** Detailed description */
  description: string;
  /** Optional bullet points for extra detail */
  details?: string[];
  /** Button size variant */
  size?: "sm" | "md";
  /** Popover alignment */
  align?: "start" | "center" | "end";
  /** Popover side */
  side?: "top" | "bottom" | "left" | "right";
}

export function HelpTooltip({
  title,
  description,
  details,
  size = "sm",
  align = "center",
  side = "bottom",
}: HelpTooltipProps) {
  const [open, setOpen] = useState(false);

  const btnSize = size === "sm" ? "w-5 h-5 text-[10px]" : "w-6 h-6 text-xs";

  return (
    <Popover open={open} onOpenChange={setOpen}>
      <PopoverTrigger asChild>
        <button
          className={`
            ${btnSize} rounded-full flex items-center justify-center
            border border-white/10 text-[#8892B0]
            hover:border-[#00B4D8]/40 hover:text-[#00B4D8] hover:bg-[#00B4D8]/5
            active:scale-95
            transition-all duration-200 flex-shrink-0
          `}
          aria-label={`도움말: ${title}`}
        >
          ?
        </button>
      </PopoverTrigger>
      <PopoverContent
        side={side}
        align={align}
        className="w-72 bg-[#112240] border-[#00B4D8]/20 shadow-[0_0_30px_rgba(0,180,216,0.1)] p-0 rounded-xl"
      >
        {/* Header */}
        <div className="px-4 pt-3 pb-2 border-b border-white/[0.06]">
          <div className="flex items-center gap-2">
            <div className="w-5 h-5 rounded-md bg-[#00B4D8]/10 flex items-center justify-center">
              <svg
                className="w-3 h-3 text-[#00B4D8]"
                fill="none"
                viewBox="0 0 24 24"
                stroke="currentColor"
                strokeWidth={2}
              >
                <path
                  strokeLinecap="round"
                  strokeLinejoin="round"
                  d="M9.879 7.519c1.171-1.025 3.071-1.025 4.242 0 1.172 1.025 1.172 2.687 0 3.712-.203.179-.43.326-.67.442-.745.361-1.45.999-1.45 1.827v.75M21 12a9 9 0 1 1-18 0 9 9 0 0 1 18 0Zm-9 5.25h.008v.008H12v-.008Z"
                />
              </svg>
            </div>
            <h4 className="text-xs font-semibold text-white">{title}</h4>
          </div>
        </div>

        {/* Body */}
        <div className="px-4 py-3 space-y-2">
          <p className="text-[11px] text-[#8892B0] leading-relaxed">
            {description}
          </p>

          {details && details.length > 0 && (
            <ul className="space-y-1.5 pt-1">
              {details.map((item, i) => (
                <li key={i} className="flex items-start gap-2 text-[11px]">
                  <span className="text-[#00B4D8] mt-0.5 flex-shrink-0">•</span>
                  <span className="text-[#8892B0]/90 leading-relaxed">{item}</span>
                </li>
              ))}
            </ul>
          )}
        </div>

        {/* Footer hint */}
        <div className="px-4 pb-2.5">
          <p className="text-[9px] text-[#8892B0]/40">
            아무 곳이나 클릭하면 닫힙니다
          </p>
        </div>
      </PopoverContent>
    </Popover>
  );
}
