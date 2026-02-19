"use client";

interface StatsCardProps {
  title: string;
  value: string | number;
  subtitle?: string;
  icon: React.ReactNode;
  trend?: {
    value: number;
    label: string;
  };
  color?: string;
}

export function StatsCard({
  title,
  value,
  subtitle,
  icon,
  trend,
  color = "#00B4D8",
}: StatsCardProps) {
  return (
    <div className="relative p-5 rounded-xl bg-[#0D1B2A] border border-white/[0.06] overflow-hidden group hover:border-white/[0.1] transition-all duration-300">
      {/* Background glow */}
      <div
        className="absolute -top-12 -right-12 w-32 h-32 rounded-full opacity-[0.03] group-hover:opacity-[0.06] transition-opacity"
        style={{ backgroundColor: color }}
      />

      <div className="flex items-start justify-between">
        <div>
          <p className="text-xs text-[#8892B0] mb-1">{title}</p>
          <p className="text-2xl font-bold font-mono text-white">{value}</p>
          {subtitle && (
            <p className="text-[11px] text-[#8892B0] mt-1">{subtitle}</p>
          )}
          {trend && (
            <div className="flex items-center gap-1 mt-2">
              <svg
                className={`w-3 h-3 ${trend.value >= 0 ? "text-[#00E676]" : "text-[#FF5252]"}`}
                fill="none"
                viewBox="0 0 24 24"
                stroke="currentColor"
                strokeWidth={2}
              >
                <path
                  strokeLinecap="round"
                  strokeLinejoin="round"
                  d={
                    trend.value >= 0
                      ? "M4.5 19.5l15-15m0 0H8.25m11.25 0v11.25"
                      : "M4.5 4.5l15 15m0 0V8.25m0 11.25H8.25"
                  }
                />
              </svg>
              <span
                className={`text-[10px] font-mono ${trend.value >= 0 ? "text-[#00E676]" : "text-[#FF5252]"}`}
              >
                {trend.value > 0 ? "+" : ""}
                {trend.value}% {trend.label}
              </span>
            </div>
          )}
        </div>
        <div
          className="w-10 h-10 rounded-xl flex items-center justify-center"
          style={{ backgroundColor: `${color}15` }}
        >
          <div style={{ color }}>{icon}</div>
        </div>
      </div>
    </div>
  );
}
