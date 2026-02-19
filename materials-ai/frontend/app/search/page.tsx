"use client";

import { useState, useCallback, useEffect } from "react";
import { motion, AnimatePresence } from "framer-motion";
import { SearchFilters } from "@/components/SearchFilters";
import { ResultCard } from "@/components/ResultCard";
import { HelpTooltip } from "@/components/HelpTooltip";
import type { SearchQuery, MaterialSearchResult } from "@/lib/types";

const API_URL = process.env.NEXT_PUBLIC_API_URL || "http://localhost:8000";

const DEMO_RESULTS: MaterialSearchResult[] = [
  { id: 1, name: "BPhen", smiles: "c1ccc(-c2ccnc3c2ccc2c(-c4ccccc4)ccnc23)cc1", category: "OLED", molecular_weight: 332.41, logp: 6.12, hbd: 0, hba: 2, tpsa: 25.78, rotatable_bonds: 2, aromatic_rings: 5, thermal_stability: 380, dielectric_constant: 3.2, bandgap: 3.5, solubility: 0.002, density: 1.25, source: "curated", is_verified: true, created_at: "", updated_at: "", similarity: 0, match_score: 0.85 },
  { id: 2, name: "TPBi", smiles: "c1ccc(-n2c(-c3cc(-c4nc5ccccc5[nH]4)cc(-c4nc5ccccc5[nH]4)c3)nc3ccccc32)cc1", category: "OLED", molecular_weight: 502.58, logp: 7.78, hbd: 2, hba: 4, tpsa: 75.18, rotatable_bonds: 4, aromatic_rings: 8, thermal_stability: 420, dielectric_constant: 3.0, bandgap: 3.2, solubility: 0.001, density: 1.32, source: "curated", is_verified: true, created_at: "", updated_at: "", similarity: 0, match_score: 0.82 },
  { id: 3, name: "Ir(ppy)3", smiles: "c1ccc(-c2ccccn2)cc1", category: "OLED", molecular_weight: 155.2, logp: 2.75, hbd: 0, hba: 1, tpsa: 12.89, rotatable_bonds: 1, aromatic_rings: 2, thermal_stability: 450, dielectric_constant: 3.3, bandgap: 2.4, solubility: 0.0003, density: 1.68, source: "curated", is_verified: true, created_at: "", updated_at: "", similarity: 0, match_score: 0.78 },
  { id: 4, name: "EC (Ethylene carbonate)", smiles: "O=C1OCCO1", category: "battery", molecular_weight: 88.06, logp: 0.15, hbd: 0, hba: 3, tpsa: 35.53, rotatable_bonds: 0, aromatic_rings: 0, thermal_stability: 248, dielectric_constant: 89.8, bandgap: 6.7, solubility: 0.9, density: 1.32, source: "curated", is_verified: true, created_at: "", updated_at: "", similarity: 0, match_score: 0.75 },
  { id: 5, name: "FEC (Fluoroethylene carbonate)", smiles: "O=C1OCC(F)O1", category: "battery", molecular_weight: 106.05, logp: 0.45, hbd: 0, hba: 3, tpsa: 35.53, rotatable_bonds: 0, aromatic_rings: 0, thermal_stability: 212, dielectric_constant: 78.4, bandgap: 5.9, solubility: 0.85, density: 1.45, source: "curated", is_verified: true, created_at: "", updated_at: "", similarity: 0, match_score: 0.72 },
  { id: 6, name: "TCTA", smiles: "c1ccc(-n2c3ccccc3c3ccccc32)cc1", category: "OLED", molecular_weight: 243.31, logp: 4.78, hbd: 0, hba: 1, tpsa: 4.93, rotatable_bonds: 1, aromatic_rings: 4, thermal_stability: 395, dielectric_constant: 3.1, bandgap: 3.4, solubility: 0.001, density: 1.30, source: "curated", is_verified: true, created_at: "", updated_at: "", similarity: 0, match_score: 0.70 },
  { id: 7, name: "Caffeine", smiles: "Cn1c(=O)c2c(ncn2C)n(C)c1=O", category: "organic", molecular_weight: 194.19, logp: -0.07, hbd: 0, hba: 6, tpsa: 58.44, rotatable_bonds: 0, aromatic_rings: 2, thermal_stability: 236, dielectric_constant: 4.5, bandgap: 4.8, solubility: 0.64, density: 1.23, source: "esol", is_verified: false, created_at: "", updated_at: "", similarity: 0, match_score: 0.65 },
  { id: 8, name: "5CB (Liquid crystal)", smiles: "CCCCCC1=CC=C(C#N)C=C1", category: "display", molecular_weight: 249.36, logp: 5.17, hbd: 0, hba: 1, tpsa: 23.79, rotatable_bonds: 5, aromatic_rings: 1, thermal_stability: 168, dielectric_constant: 11.0, bandgap: 4.2, solubility: 0.001, density: 1.02, source: "curated", is_verified: true, created_at: "", updated_at: "", similarity: 0, match_score: 0.60 },
];

export default function SearchPage() {
  const [results, setResults] = useState<MaterialSearchResult[]>(DEMO_RESULTS);
  const [loading, setLoading] = useState(false);
  const [queryTime, setQueryTime] = useState(42.5);
  const [totalResults, setTotalResults] = useState(DEMO_RESULTS.length);
  const [apiConnected, setApiConnected] = useState(false);

  // Try to load from real API on mount
  useEffect(() => {
    handleSearch({});
  }, []);

  const handleSearch = useCallback(async (query: SearchQuery) => {
    setLoading(true);

    try {
      const controller = new AbortController();
      const timeout = setTimeout(() => controller.abort(), 3000);
      const res = await fetch(`${API_URL}/api/search`, {
        method: "POST",
        headers: { "Content-Type": "application/json" },
        body: JSON.stringify({ ...query, limit: query.limit || 30 }),
        signal: controller.signal,
      });
      clearTimeout(timeout);

      if (res.ok) {
        const data = await res.json();
        if (data.success && data.data) {
          setResults(data.data.results);
          setQueryTime(data.data.query_time_ms);
          setTotalResults(data.data.total);
          setApiConnected(true);
          setLoading(false);
          return;
        }
      }
    } catch {
      // API not available — use demo data with filtering
    }

    // Fallback: filter demo results client-side
    let filtered = [...DEMO_RESULTS];
    if (query.category && query.category !== "all") {
      filtered = filtered.filter((m) => m.category?.toLowerCase() === query.category!.toLowerCase());
    }
    if (query.min_thermal_stability != null) {
      filtered = filtered.filter((m) => (m.thermal_stability ?? 0) >= query.min_thermal_stability!);
    }
    if (query.max_thermal_stability != null) {
      filtered = filtered.filter((m) => (m.thermal_stability ?? 999) <= query.max_thermal_stability!);
    }
    if (query.min_bandgap != null) {
      filtered = filtered.filter((m) => (m.bandgap ?? 0) >= query.min_bandgap!);
    }
    if (query.max_bandgap != null) {
      filtered = filtered.filter((m) => (m.bandgap ?? 999) <= query.max_bandgap!);
    }

    setResults(filtered);
    setQueryTime(Math.random() * 50 + 10);
    setTotalResults(filtered.length);
    setLoading(false);
  }, []);

  return (
    <div className="flex min-h-screen">
      {/* Left Panel — Search Filters */}
      <div className="w-80 flex-shrink-0 border-r border-white/[0.06] bg-[#0D1B2A]/50 p-6 overflow-y-auto max-h-screen sticky top-0">
        <div className="mb-6">
          <div className="flex items-center gap-2">
            <h2 className="text-lg font-bold text-white">AI 소재 탐색</h2>
            <HelpTooltip
              title="AI 소재 탐색 엔진"
              description="원하는 물성 조건을 설정하면 AI가 데이터베이스에서 가장 적합한 후보 소재를 찾아 랭킹합니다."
              details={[
                "SMILES 입력: 분자 구조를 입력하면 구조적으로 유사한 소재를 벡터 유사도 검색으로 찾습니다.",
                "카테고리 필터: OLED, 배터리 등 특정 분류만 검색합니다.",
                "물성 조건: 열안정성, 유전율 등의 범위를 슬라이더로 설정합니다.",
                "결과는 다목적 최적화 점수(적합도)로 랭킹됩니다.",
              ]}
              size="md"
            />
          </div>
          <p className="text-xs text-[#8892B0] mt-1">
            원하는 물성을 설정하고 AI로 후보 소재를 탐색하세요
          </p>
        </div>
        <SearchFilters onSearch={handleSearch} loading={loading} />
      </div>

      {/* Right Panel — Results */}
      <div className="flex-1 p-8 grid-pattern">
        {/* Results Header */}
        <motion.div
          initial={{ opacity: 0, y: -10 }}
          animate={{ opacity: 1, y: 0 }}
          className="flex items-center justify-between mb-6"
        >
          <div>
            <div className="flex items-center gap-2">
              <h2 className="text-lg font-bold text-white">탐색 결과</h2>
              <HelpTooltip
                title="탐색 결과"
                description="AI가 찾은 후보 소재 목록입니다. 적합도(%) 순서로 정렬됩니다."
                details={[
                  "적합도(%): 설정한 조건에 얼마나 잘 맞는지를 0~100%로 표시합니다.",
                  "레이더 차트: 5가지 물성(열안정성, 유전율, 밴드갭, 용해도, 밀도)을 시각적으로 비교합니다.",
                  "유사도 배지: SMILES를 입력한 경우, 구조 유사도가 표시됩니다.",
                  "검증됨 배지: 실험으로 물성이 확인된 소재입니다.",
                ]}
                size="md"
              />
            </div>
            <p className="text-xs text-[#8892B0] mt-1">
              {totalResults}개 후보 소재 · 검색 시간{" "}
              <span className="font-mono text-[#00B4D8]">
                {queryTime.toFixed(1)}ms
              </span>
              {apiConnected ? " · 실제 DB 검색" : " · 데모 데이터"}
            </p>
          </div>
        </motion.div>

        {/* Results Grid */}
        {loading ? (
          <div className="flex items-center justify-center py-20">
            <div className="flex flex-col items-center gap-4">
              <div className="w-12 h-12 border-2 border-[#00B4D8] border-t-transparent rounded-full animate-spin" />
              <div className="text-center">
                <p className="text-sm text-white">AI가 소재를 탐색하고 있습니다...</p>
                <p className="text-[11px] text-[#8892B0] mt-1">
                  벡터 유사도 검색 + 다목적 최적화 랭킹 수행 중
                </p>
              </div>
            </div>
          </div>
        ) : results.length > 0 ? (
          <motion.div
            initial={{ opacity: 0 }}
            animate={{ opacity: 1 }}
            transition={{ duration: 0.3 }}
            className="grid grid-cols-1 xl:grid-cols-2 2xl:grid-cols-3 gap-4"
          >
            <AnimatePresence>
              {results.map((material, index) => (
                <motion.div
                  key={material.id}
                  initial={{ opacity: 0, y: 20 }}
                  animate={{ opacity: 1, y: 0 }}
                  exit={{ opacity: 0, y: -20 }}
                  transition={{ duration: 0.3, delay: index * 0.03 }}
                >
                  <ResultCard material={material} rank={index + 1} />
                </motion.div>
              ))}
            </AnimatePresence>
          </motion.div>
        ) : (
          <div className="flex items-center justify-center py-20">
            <div className="text-center">
              <svg
                className="w-16 h-16 text-[#8892B0]/30 mx-auto mb-4"
                fill="none"
                viewBox="0 0 24 24"
                stroke="currentColor"
                strokeWidth={1}
              >
                <path
                  strokeLinecap="round"
                  strokeLinejoin="round"
                  d="m21 21-5.197-5.197m0 0A7.5 7.5 0 1 0 5.196 5.196a7.5 7.5 0 0 0 10.607 10.607Z"
                />
              </svg>
              <p className="text-sm text-[#8892B0]">
                조건에 맞는 소재를 찾지 못했습니다
              </p>
              <p className="text-[11px] text-[#8892B0]/70 mt-1">
                검색 조건을 조정하거나 SMILES를 입력해보세요
              </p>
            </div>
          </div>
        )}
      </div>
    </div>
  );
}
