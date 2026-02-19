"use client";

import {
  Radar,
  RadarChart,
  PolarGrid,
  PolarAngleAxis,
  PolarRadiusAxis,
  ResponsiveContainer,
  Legend,
} from "recharts";

interface PropertyRadarProps {
  properties: {
    thermal_stability: number | null;
    dielectric_constant: number | null;
    bandgap: number | null;
    solubility: number | null;
    density: number | null;
  };
  /** Optional actual (measured) properties for comparison */
  actualProperties?: {
    thermal_stability?: number | null;
    dielectric_constant?: number | null;
    bandgap?: number | null;
    solubility?: number | null;
    density?: number | null;
  };
}

// Normalization ranges for each property
const RANGES: Record<string, { min: number; max: number; label: string }> = {
  thermal_stability: { min: 0, max: 500, label: "열안정성" },
  dielectric_constant: { min: 0, max: 25, label: "유전율" },
  bandgap: { min: 0, max: 8, label: "밴드갭" },
  solubility: { min: 0, max: 1, label: "용해도" },
  density: { min: 0.5, max: 3.0, label: "밀도" },
};

function normalize(key: string, value: number | null | undefined): number {
  if (value === null || value === undefined) return 0;
  const range = RANGES[key];
  if (!range) return 0;
  return Math.max(
    0,
    Math.min(100, ((value - range.min) / (range.max - range.min)) * 100)
  );
}

export function PropertyRadar({
  properties,
  actualProperties,
}: PropertyRadarProps) {
  const data = Object.entries(RANGES).map(([key, { label }]) => ({
    property: label,
    predicted: normalize(
      key,
      properties[key as keyof typeof properties] as number | null
    ),
    ...(actualProperties
      ? {
          actual: normalize(
            key,
            actualProperties[key as keyof typeof actualProperties] as
              | number
              | null
              | undefined
          ),
        }
      : {}),
  }));

  return (
    <ResponsiveContainer width="100%" height={280}>
      <RadarChart cx="50%" cy="50%" outerRadius="70%" data={data}>
        <PolarGrid stroke="rgba(255,255,255,0.06)" />
        <PolarAngleAxis
          dataKey="property"
          tick={{ fill: "#8892B0", fontSize: 11 }}
        />
        <PolarRadiusAxis
          angle={90}
          domain={[0, 100]}
          tick={false}
          axisLine={false}
        />
        <Radar
          name="예측값"
          dataKey="predicted"
          stroke="#00B4D8"
          fill="#00B4D8"
          fillOpacity={0.2}
          strokeWidth={2}
        />
        {actualProperties && (
          <Radar
            name="실측값"
            dataKey="actual"
            stroke="#00E676"
            fill="#00E676"
            fillOpacity={0.1}
            strokeWidth={2}
            strokeDasharray="4 4"
          />
        )}
        <Legend
          wrapperStyle={{ fontSize: 11, color: "#8892B0" }}
        />
      </RadarChart>
    </ResponsiveContainer>
  );
}
