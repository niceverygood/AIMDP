"use client";

import { useState, useCallback } from "react";
import { Slider } from "@/components/ui/slider";
import { Input } from "@/components/ui/input";
import { Label } from "@/components/ui/label";
import { Button } from "@/components/ui/button";
import { Badge } from "@/components/ui/badge";
import {
  Select,
  SelectContent,
  SelectItem,
  SelectTrigger,
  SelectValue,
} from "@/components/ui/select";
import { HelpTooltip } from "@/components/HelpTooltip";
import { PROPERTY_FILTERS, MATERIAL_CATEGORIES } from "@/lib/types";
import type { SearchQuery } from "@/lib/types";

interface SearchFiltersProps {
  onSearch: (query: SearchQuery) => void;
  loading?: boolean;
}

type FilterRange = {
  min: number;
  max: number;
  enabled: boolean;
};

export function SearchFilters({ onSearch, loading }: SearchFiltersProps) {
  const [smiles, setSmiles] = useState("");
  const [category, setCategory] = useState<string>("");
  const [filters, setFilters] = useState<Record<string, FilterRange>>(() => {
    const init: Record<string, FilterRange> = {};
    for (const f of PROPERTY_FILTERS) {
      init[f.key] = { min: f.min, max: f.max, enabled: false };
    }
    return init;
  });

  const toggleFilter = useCallback((key: string) => {
    setFilters((prev) => ({
      ...prev,
      [key]: { ...prev[key], enabled: !prev[key].enabled },
    }));
  }, []);

  const updateFilter = useCallback((key: string, values: number[]) => {
    setFilters((prev) => ({
      ...prev,
      [key]: { ...prev[key], min: values[0], max: values[1] },
    }));
  }, []);

  const handleSearch = useCallback(() => {
    const query: SearchQuery = {
      limit: 20,
      offset: 0,
    };

    if (smiles.trim()) query.smiles = smiles.trim();
    if (category) query.category = category;

    for (const f of PROPERTY_FILTERS) {
      const range = filters[f.key];
      if (range.enabled) {
        const minKey = `min_${f.key}` as keyof SearchQuery;
        const maxKey = `max_${f.key}` as keyof SearchQuery;
        (query as Record<string, unknown>)[minKey] = range.min;
        (query as Record<string, unknown>)[maxKey] = range.max;
      }
    }

    onSearch(query);
  }, [smiles, category, filters, onSearch]);

  return (
    <div className="space-y-6">
      {/* SMILES Input */}
      <div className="space-y-2">
        <div className="flex items-center gap-1.5">
          <Label className="text-sm text-[#8892B0]">분자 구조 (SMILES)</Label>
          <HelpTooltip
            title="SMILES 입력"
            description="분자의 구조를 텍스트로 표현하는 SMILES 표기법입니다. 입력하면 구조적으로 유사한 소재를 찾습니다."
            details={[
              "예: c1ccccc1 (벤젠), CCO (에탄올)",
              "예: c1ccc2[nH]c3ccccc3c2c1 (카바졸)",
              "빈칸으로 두면 물성 조건만으로 검색합니다.",
              "유사도 검색은 Morgan Fingerprint 기반 코사인 유사도를 사용합니다.",
            ]}
            side="right"
          />
        </div>
        <div className="flex gap-2">
          <Input
            value={smiles}
            onChange={(e) => setSmiles(e.target.value)}
            placeholder="예: CCO, c1ccccc1, CC(=O)O"
            className="bg-[#112240] border-white/[0.08] text-white placeholder:text-[#8892B0]/50 font-mono text-sm"
          />
        </div>
      </div>

      {/* Category Filter */}
      <div className="space-y-2">
        <div className="flex items-center gap-1.5">
          <Label className="text-sm text-[#8892B0]">소재 카테고리</Label>
          <HelpTooltip
            title="카테고리 필터"
            description="특정 소재 유형만 검색합니다. '전체'를 선택하면 모든 카테고리에서 검색합니다."
            details={[
              "OLED: 유기발광 소재 (전자수송층, 호스트, 발광체)",
              "배터리: 전해액, 첨가제 (EC, DMC, FEC 등)",
              "반도체: 절연막, 포토레지스트 소재",
              "하드코팅: 실란 계열 코팅 소재",
            ]}
            side="right"
          />
        </div>
        <Select value={category} onValueChange={setCategory}>
          <SelectTrigger className="bg-[#112240] border-white/[0.08] text-white">
            <SelectValue placeholder="전체 카테고리" />
          </SelectTrigger>
          <SelectContent className="bg-[#112240] border-white/[0.08]">
            <SelectItem value="all" className="text-white">
              전체
            </SelectItem>
            {MATERIAL_CATEGORIES.map((cat) => (
              <SelectItem
                key={cat.value}
                value={cat.value}
                className="text-white"
              >
                {cat.label}
              </SelectItem>
            ))}
          </SelectContent>
        </Select>
      </div>

      {/* Property Filters */}
      <div className="space-y-4">
        <div className="flex items-center gap-1.5">
          <Label className="text-sm text-[#8892B0]">물성 조건</Label>
          <HelpTooltip
            title="물성 조건 필터"
            description="원하는 물성 범위를 설정합니다. 체크박스를 클릭하면 해당 조건이 활성화되고, 슬라이더로 범위를 조절할 수 있습니다."
            details={[
              "열안정성(°C): 소재가 분해되기 시작하는 온도. OLED 소재는 보통 300°C 이상.",
              "유전율: 전기장에 대한 반응 정도. 배터리 전해질은 높을수록 좋음.",
              "밴드갭(eV): HOMO-LUMO 에너지 차이. OLED 호스트는 3.0+ eV 필요.",
              "용해도: 물에 대한 용해도. 0(불용)~1(완전용해).",
            ]}
            side="right"
          />
        </div>
        {PROPERTY_FILTERS.map((f) => {
          const range = filters[f.key];
          return (
            <div
              key={f.key}
              className={`p-3 rounded-lg border transition-all duration-200 ${
                range.enabled
                  ? "bg-[#00B4D8]/5 border-[#00B4D8]/20"
                  : "bg-[#112240]/50 border-white/[0.04]"
              }`}
            >
              <div className="flex items-center justify-between mb-2">
                <button
                  onClick={() => toggleFilter(f.key)}
                  className="flex items-center gap-2 text-sm"
                >
                  <div
                    className={`w-4 h-4 rounded border flex items-center justify-center transition-colors ${
                      range.enabled
                        ? "bg-[#00B4D8] border-[#00B4D8]"
                        : "border-white/20"
                    }`}
                  >
                    {range.enabled && (
                      <svg
                        className="w-3 h-3 text-white"
                        fill="none"
                        viewBox="0 0 24 24"
                        stroke="currentColor"
                        strokeWidth={3}
                      >
                        <path
                          strokeLinecap="round"
                          strokeLinejoin="round"
                          d="m4.5 12.75 6 6 9-13.5"
                        />
                      </svg>
                    )}
                  </div>
                  <span
                    className={
                      range.enabled ? "text-white" : "text-[#8892B0]"
                    }
                  >
                    {f.label}
                  </span>
                </button>
                {range.enabled && (
                  <Badge
                    variant="outline"
                    className="text-[10px] border-[#00B4D8]/30 text-[#00B4D8]"
                  >
                    {range.min}
                    {f.unit} — {range.max}
                    {f.unit}
                  </Badge>
                )}
              </div>

              {range.enabled && (
                <div className="pt-2">
                  <Slider
                    min={f.min}
                    max={f.max}
                    step={f.step}
                    value={[range.min, range.max]}
                    onValueChange={(values) => updateFilter(f.key, values)}
                    className="w-full"
                  />
                  <div className="flex justify-between mt-1">
                    <span className="text-[10px] text-[#8892B0] font-mono">
                      {f.min}
                      {f.unit}
                    </span>
                    <span className="text-[10px] text-[#8892B0] font-mono">
                      {f.max}
                      {f.unit}
                    </span>
                  </div>
                </div>
              )}
            </div>
          );
        })}
      </div>

      {/* Search Button */}
      <Button
        onClick={handleSearch}
        disabled={loading}
        className="w-full bg-gradient-to-r from-[#00B4D8] to-[#0096C7] hover:from-[#0096C7] hover:to-[#0077B6] text-white font-medium"
      >
        {loading ? (
          <div className="flex items-center gap-2">
            <div className="w-4 h-4 border-2 border-white border-t-transparent rounded-full animate-spin" />
            <span>탐색 중...</span>
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
                d="m21 21-5.197-5.197m0 0A7.5 7.5 0 1 0 5.196 5.196a7.5 7.5 0 0 0 10.607 10.607Z"
              />
            </svg>
            <span>AI 소재 탐색</span>
          </div>
        )}
      </Button>
    </div>
  );
}
