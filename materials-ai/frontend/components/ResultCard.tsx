"use client";

import { Card, CardContent, CardHeader } from "@/components/ui/card";
import { Badge } from "@/components/ui/badge";
import { PropertyRadar } from "./PropertyRadar";
import type { MaterialSearchResult } from "@/lib/types";

interface ResultCardProps {
  material: MaterialSearchResult;
  rank: number;
}

export function ResultCard({ material, rank }: ResultCardProps) {
  const scoreColor =
    material.match_score >= 0.8
      ? "text-[#00E676]"
      : material.match_score >= 0.5
        ? "text-[#00B4D8]"
        : material.match_score >= 0.3
          ? "text-[#FF9100]"
          : "text-[#FF5252]";

  const scoreBarColor =
    material.match_score >= 0.8
      ? "bg-[#00E676]"
      : material.match_score >= 0.5
        ? "bg-[#00B4D8]"
        : material.match_score >= 0.3
          ? "bg-[#FF9100]"
          : "bg-[#FF5252]";

  return (
    <Card className="bg-[#0D1B2A] border-white/[0.06] hover:border-[#00B4D8]/20 transition-all duration-300 hover:shadow-[0_0_30px_rgba(0,180,216,0.05)] group">
      <CardHeader className="pb-3">
        <div className="flex items-start justify-between">
          <div className="flex items-center gap-3">
            <div className="w-8 h-8 rounded-lg bg-[#112240] flex items-center justify-center text-sm font-mono font-bold text-[#00B4D8]">
              {rank}
            </div>
            <div>
              <h3 className="text-sm font-semibold text-white">
                {material.name || `Material #${material.id}`}
              </h3>
              <code className="text-[11px] text-[#8892B0] font-mono">
                {material.smiles.length > 40
                  ? material.smiles.slice(0, 40) + "..."
                  : material.smiles}
              </code>
            </div>
          </div>
          <div className="text-right">
            <div className={`text-lg font-bold font-mono ${scoreColor}`}>
              {(material.match_score * 100).toFixed(1)}%
            </div>
            <p className="text-[10px] text-[#8892B0]">적합도</p>
          </div>
        </div>

        {/* Score Bar */}
        <div className="w-full h-1 bg-white/[0.04] rounded-full mt-2">
          <div
            className={`h-full rounded-full transition-all duration-500 ${scoreBarColor}`}
            style={{ width: `${material.match_score * 100}%` }}
          />
        </div>
      </CardHeader>

      <CardContent className="space-y-4">
        {/* Badges */}
        <div className="flex flex-wrap gap-1.5">
          {material.category && (
            <Badge
              variant="outline"
              className="text-[10px] border-[#7C4DFF]/30 text-[#7C4DFF] bg-[#7C4DFF]/5"
            >
              {material.category}
            </Badge>
          )}
          {material.is_verified && (
            <Badge className="text-[10px] bg-[#00E676]/10 text-[#00E676] border-[#00E676]/30">
              검증됨
            </Badge>
          )}
          {material.similarity > 0.01 && (
            <Badge
              variant="outline"
              className="text-[10px] border-[#00B4D8]/30 text-[#00B4D8]"
            >
              유사도 {(material.similarity * 100).toFixed(1)}%
            </Badge>
          )}
          {material.source && (
            <Badge
              variant="outline"
              className="text-[10px] border-white/10 text-[#8892B0]"
            >
              {material.source}
            </Badge>
          )}
        </div>

        {/* Property Radar Chart */}
        <PropertyRadar
          properties={{
            thermal_stability: material.thermal_stability,
            dielectric_constant: material.dielectric_constant,
            bandgap: material.bandgap,
            solubility: material.solubility,
            density: material.density,
          }}
        />

        {/* Property Values */}
        <div className="grid grid-cols-2 gap-2">
          <PropertyItem
            label="열안정성"
            value={material.thermal_stability}
            unit="°C"
            color="#FF9100"
          />
          <PropertyItem
            label="유전율"
            value={material.dielectric_constant}
            unit=""
            color="#7C4DFF"
          />
          <PropertyItem
            label="밴드갭"
            value={material.bandgap}
            unit="eV"
            color="#00B4D8"
          />
          <PropertyItem
            label="용해도"
            value={material.solubility}
            unit=""
            color="#00E676"
          />
          <PropertyItem
            label="밀도"
            value={material.density}
            unit="g/cm³"
            color="#FF5252"
          />
          <PropertyItem
            label="분자량"
            value={material.molecular_weight}
            unit="g/mol"
            color="#8892B0"
          />
        </div>
      </CardContent>
    </Card>
  );
}

function PropertyItem({
  label,
  value,
  unit,
  color,
}: {
  label: string;
  value: number | null;
  unit: string;
  color: string;
}) {
  return (
    <div className="flex items-center justify-between px-3 py-2 rounded-lg bg-[#112240]/50">
      <span className="text-[11px] text-[#8892B0]">{label}</span>
      <span className="text-xs font-mono font-medium" style={{ color }}>
        {value !== null ? `${value.toFixed(2)} ${unit}` : "—"}
      </span>
    </div>
  );
}
