"use client";

import { useState, useEffect, useCallback } from "react";
import { motion } from "framer-motion";
import { Button } from "@/components/ui/button";
import { Card, CardContent, CardHeader, CardTitle } from "@/components/ui/card";
import { Badge } from "@/components/ui/badge";
import { Input } from "@/components/ui/input";
import { Label } from "@/components/ui/label";
import { Progress } from "@/components/ui/progress";
import { Select, SelectContent, SelectItem, SelectTrigger, SelectValue } from "@/components/ui/select";
import { HelpTooltip } from "@/components/HelpTooltip";
import { getDataSummary } from "@/lib/db";
import { supabase } from "@/lib/supabase";

interface DataSummary {
  total_materials: number;
  verified: number;
  by_source: Record<string, number>;
  by_category: Record<string, number>;
  property_coverage: Record<string, { count: number; ratio: number }>;
}

const SOURCE_LABELS: Record<string, string> = {
  dongjin_internal: "동진쎄미켐 내부",
  curated: "큐레이션 데이터",
  esol: "ESOL 데이터셋",
  generated: "자동 생성",
  unknown: "기타",
};

const PROP_LABELS: Record<string, string> = {
  thermal_stability: "열안정성",
  dielectric_constant: "유전율",
  bandgap: "밴드갭",
  solubility: "용해도",
  density: "밀도",
};

export default function DataPage() {
  const [summary, setSummary] = useState<DataSummary | null>(null);
  const [uploading, setUploading] = useState(false);
  const [demoLoading, setDemoLoading] = useState(false);
  const [uploadResult, setUploadResult] = useState<Record<string, unknown> | null>(null);
  const [category, setCategory] = useState("");

  const fetchSummary = useCallback(async () => {
    try {
      const data = await getDataSummary();
      if (data.total_materials > 0) {
        setSummary(data);
      }
    } catch (err) {
      console.error("Data summary error:", err);
    }
  }, []);

  useEffect(() => { fetchSummary(); }, [fetchSummary]);

  const handleUpload = async (e: React.ChangeEvent<HTMLInputElement>) => {
    const file = e.target.files?.[0];
    if (!file) return;

    setUploading(true);
    setUploadResult(null);

    const formData = new FormData();
    formData.append("file", file);
    formData.append("category", category);
    formData.append("source", "dongjin_internal");

    try {
      const r = await fetch("/api/data/upload", { method: "POST", body: formData });
      const d = await r.json();
      setUploadResult(d.success ? d.data : { error: d.error || "업로드 실패" });
      if (d.success) fetchSummary();
    } catch (err) {
      setUploadResult({ error: "서버 연결 실패" });
    } finally {
      setUploading(false);
      e.target.value = "";
    }
  };

  return (
    <div className="p-8 space-y-8 grid-pattern min-h-screen">
      <motion.div initial={{ opacity: 0, y: -20 }} animate={{ opacity: 1, y: 0 }}>
        <div className="flex items-center gap-2">
          <h1 className="text-2xl font-bold text-white">연구 데이터 관리</h1>
          <HelpTooltip
            title="연구 데이터 관리"
            description="동진쎄미켐의 연구 데이터를 업로드하여 AI가 학습할 수 있도록 합니다."
            details={[
              "CSV/Excel 파일을 업로드하면 자동으로 파싱하여 DB에 저장합니다.",
              "SMILES 컬럼이 필수이며, 물성 데이터가 있으면 함께 저장됩니다.",
              "RDKit로 분자 기술자를 자동 계산하고, 벡터 임베딩을 생성합니다.",
              "업로드된 데이터가 많을수록 AI 신소재 제안의 정확도가 높아집니다.",
            ]}
            size="md"
          />
        </div>
        <p className="text-sm text-[#8892B0] mt-1">AI 학습을 위한 연구 데이터 업로드 및 현황</p>
      </motion.div>

      <div className="grid grid-cols-1 lg:grid-cols-3 gap-6">
        {/* Upload Panel */}
        <Card className="bg-[#0D1B2A] border-white/[0.06]">
          <CardHeader>
            <div className="flex items-center gap-2">
              <CardTitle className="text-sm font-semibold text-white">데이터 업로드</CardTitle>
              <HelpTooltip
                title="파일 업로드"
                description="CSV 또는 Excel 파일을 업로드하세요. SMILES 컬럼이 필수입니다."
                details={[
                  "지원 형식: .csv, .xlsx, .xls, .tsv",
                  "필수 컬럼: smiles (또는 SMILES)",
                  "선택 컬럼: name, category, thermal_stability, bandgap, dielectric_constant, solubility, density",
                  "한글 컬럼명도 자동 인식: 열안정성, 유전율, 밴드갭 등",
                ]}
              />
            </div>
          </CardHeader>
          <CardContent className="space-y-4">
            {/* Demo Data Button */}
            <Button
              onClick={async () => {
                setDemoLoading(true);
                setUploadResult(null);
                try {
                  const demoData = [
                    { name: "BPhen", smiles: "c1ccc(-c2ccnc3c2ccc2c(-c4ccccc4)ccnc23)cc1", category: "OLED", thermal_stability: 380, dielectric_constant: 3.2, bandgap: 3.5, solubility: 0.002, density: 1.25, source: "demo", is_verified: true },
                    { name: "TPBi", smiles: "c1ccc(-n2c(-c3cc(-c4nc5ccccc5[nH]4)cc(-c4nc5ccccc5[nH]4)c3)nc3ccccc32)cc1", category: "OLED", thermal_stability: 420, dielectric_constant: 3.0, bandgap: 3.2, solubility: 0.001, density: 1.32, source: "demo", is_verified: true },
                    { name: "Ir(ppy)3", smiles: "c1ccc(-c2ccccn2)cc1", category: "OLED", thermal_stability: 450, dielectric_constant: 3.3, bandgap: 2.4, solubility: 0.0003, density: 1.68, source: "demo", is_verified: true },
                    { name: "TCTA", smiles: "c1ccc(-n2c3ccccc3c3ccccc32)cc1", category: "OLED", thermal_stability: 395, dielectric_constant: 3.1, bandgap: 3.4, solubility: 0.001, density: 1.30, source: "demo", is_verified: true },
                    { name: "DPEPO", smiles: "O=P(c1ccccc1)(c1ccccc1)c1ccccc1", category: "OLED", thermal_stability: 350, dielectric_constant: 3.8, bandgap: 4.1, solubility: 0.003, density: 1.28, source: "demo", is_verified: true },
                    { name: "EC", smiles: "O=C1OCCO1", category: "battery", thermal_stability: 248, dielectric_constant: 89.8, bandgap: 6.7, solubility: 0.9, density: 1.32, source: "demo", is_verified: true },
                    { name: "FEC", smiles: "O=C1OCC(F)O1", category: "battery", thermal_stability: 212, dielectric_constant: 78.4, bandgap: 5.9, solubility: 0.85, density: 1.45, source: "demo", is_verified: true },
                    { name: "VC", smiles: "O=c1occo1", category: "battery", thermal_stability: 162, dielectric_constant: 126.0, bandgap: 5.5, solubility: 0.9, density: 1.36, source: "demo", is_verified: true },
                    { name: "LiTFSI", smiles: "O=S(=O)(C(F)(F)F)N([Li])S(=O)(=O)C(F)(F)F", category: "battery", thermal_stability: 360, solubility: 0.95, density: 1.33, source: "demo", is_verified: true },
                    { name: "PMDA", smiles: "O=c1oc(=O)c2cc3c(=O)oc(=O)c3cc12", category: "semiconductor", thermal_stability: 480, dielectric_constant: 3.4, bandgap: 3.0, solubility: 0.0001, density: 1.42, source: "demo", is_verified: true },
                    { name: "TEOS", smiles: "CCO[Si](OCC)(OCC)OCC", category: "hard_coating", thermal_stability: 165, dielectric_constant: 3.9, bandgap: 8.9, solubility: 0.01, density: 0.93, source: "demo", is_verified: true },
                    { name: "5CB", smiles: "CCCCCC1=CC=C(C#N)C=C1", category: "display", thermal_stability: 168, dielectric_constant: 11.0, bandgap: 4.2, solubility: 0.001, density: 1.02, source: "demo", is_verified: true },
                    { name: "Caffeine", smiles: "Cn1c(=O)c2c(ncn2C)n(C)c1=O", category: "organic", thermal_stability: 236, dielectric_constant: 4.5, bandgap: 4.8, solubility: 0.64, density: 1.23, source: "demo", is_verified: false },
                    { name: "Carbazole", smiles: "c1ccc2c(c1)[nH]c1ccccc12", category: "organic", thermal_stability: 340, dielectric_constant: 2.9, bandgap: 3.6, solubility: 0.15, density: 1.30, source: "demo", is_verified: false },
                    { name: "Naphthalene", smiles: "c1ccc2ccccc2c1", category: "organic", thermal_stability: 250, dielectric_constant: 2.5, bandgap: 3.6, solubility: 0.35, density: 1.14, source: "demo", is_verified: false },
                  ];

                  const { data, error } = await supabase.from("materials").insert(demoData).select("id");
                  if (error) {
                    setUploadResult({ error: error.message });
                  } else {
                    const { count } = await supabase.from("materials").select("*", { count: "exact", head: true });
                    setUploadResult({ filename: "데모 데이터", inserted: data?.length || 0, rows_total: demoData.length, total_in_db: count || 0 });
                    fetchSummary();
                  }
                } catch (err) {
                  setUploadResult({ error: String(err) });
                } finally {
                  setDemoLoading(false);
                }
              }}
              disabled={uploading || demoLoading}
              variant="outline"
              className="w-full border-[#00B4D8]/20 text-[#00B4D8] hover:bg-[#00B4D8]/10 text-xs"
            >
              {demoLoading ? (
                <div className="flex items-center gap-2">
                  <div className="w-3.5 h-3.5 border-2 border-[#00B4D8] border-t-transparent rounded-full animate-spin" />
                  저장 중...
                </div>
              ) : (
                <div className="flex items-center gap-2">
                  <svg className="w-4 h-4" fill="none" viewBox="0 0 24 24" stroke="currentColor" strokeWidth={1.5}>
                    <path strokeLinecap="round" strokeLinejoin="round" d="M20.25 6.375c0 2.278-3.694 4.125-8.25 4.125S3.75 8.653 3.75 6.375m16.5 0c0-2.278-3.694-4.125-8.25-4.125S3.75 4.097 3.75 6.375m16.5 0v11.25c0 2.278-3.694 4.125-8.25 4.125s-8.25-1.847-8.25-4.125V6.375" />
                  </svg>
                  데모 데이터 불러오기 (15건)
                </div>
              )}
            </Button>

            <div className="flex items-center gap-2">
              <div className="flex-1 h-px bg-white/[0.06]" />
              <span className="text-[10px] text-[#8892B0]">또는 직접 업로드</span>
              <div className="flex-1 h-px bg-white/[0.06]" />
            </div>

            <div className="space-y-2">
              <Label className="text-xs text-[#8892B0]">소재 카테고리</Label>
              <Select value={category} onValueChange={setCategory}>
                <SelectTrigger className="bg-[#112240] border-white/[0.08] text-white">
                  <SelectValue placeholder="자동 감지" />
                </SelectTrigger>
                <SelectContent className="bg-[#112240] border-white/[0.08]">
                  <SelectItem value="auto" className="text-white">자동 감지</SelectItem>
                  <SelectItem value="OLED" className="text-white">OLED</SelectItem>
                  <SelectItem value="battery" className="text-white">배터리</SelectItem>
                  <SelectItem value="semiconductor" className="text-white">반도체</SelectItem>
                  <SelectItem value="hard_coating" className="text-white">하드코팅</SelectItem>
                  <SelectItem value="display" className="text-white">디스플레이</SelectItem>
                </SelectContent>
              </Select>
            </div>

            <div className="space-y-2">
              <Label className="text-xs text-[#8892B0]">연구 데이터 파일</Label>
              <div className="relative">
                <input
                  type="file"
                  accept=".csv,.xlsx,.xls,.tsv"
                  onChange={handleUpload}
                  disabled={uploading}
                  className="absolute inset-0 opacity-0 cursor-pointer z-10"
                />
                <div className="p-6 border-2 border-dashed border-white/[0.08] rounded-xl text-center hover:border-[#00B4D8]/30 transition-colors">
                  {uploading ? (
                    <div className="flex flex-col items-center gap-2">
                      <div className="w-8 h-8 border-2 border-[#00B4D8] border-t-transparent rounded-full animate-spin" />
                      <span className="text-xs text-[#8892B0]">파싱 및 저장 중...</span>
                    </div>
                  ) : (
                    <>
                      <svg className="w-8 h-8 text-[#8892B0]/50 mx-auto mb-2" fill="none" viewBox="0 0 24 24" stroke="currentColor" strokeWidth={1.5}>
                        <path strokeLinecap="round" strokeLinejoin="round" d="M3 16.5v2.25A2.25 2.25 0 0 0 5.25 21h13.5A2.25 2.25 0 0 0 21 18.75V16.5m-13.5-9L12 3m0 0 4.5 4.5M12 3v13.5" />
                      </svg>
                      <p className="text-xs text-[#8892B0]">클릭하여 파일 선택</p>
                      <p className="text-[10px] text-[#8892B0]/50 mt-1">CSV, Excel 지원</p>
                    </>
                  )}
                </div>
              </div>
            </div>

            {uploadResult && (
              <motion.div initial={{ opacity: 0, y: 10 }} animate={{ opacity: 1, y: 0 }}
                className={`p-3 rounded-lg border text-xs ${
                  "error" in uploadResult
                    ? "bg-[#FF5252]/10 border-[#FF5252]/20 text-[#FF5252]"
                    : "bg-[#00E676]/10 border-[#00E676]/20 text-[#00E676]"
                }`}
              >
                {"error" in uploadResult ? (
                  <p>업로드 실패: {String(uploadResult.error)}</p>
                ) : (
                  <div className="space-y-1">
                    <p className="font-medium">업로드 완료!</p>
                    <p>파일: {String(uploadResult.filename)}</p>
                    <p>등록: {String(uploadResult.inserted)}건 / 전체 {String(uploadResult.rows_total)}행</p>
                    <p>DB 총 소재: {String(uploadResult.total_in_db)}건</p>
                  </div>
                )}
              </motion.div>
            )}
          </CardContent>
        </Card>

        {/* Data Summary */}
        <div className="lg:col-span-2 space-y-6">
          <Card className="bg-[#0D1B2A] border-white/[0.06]">
            <CardHeader>
              <CardTitle className="text-sm font-semibold text-white">데이터 현황</CardTitle>
            </CardHeader>
            <CardContent>
              {summary ? (
                <div className="space-y-6">
                  {/* Stats */}
                  <div className="grid grid-cols-3 gap-4">
                    <div className="p-4 rounded-lg bg-[#112240]/50 text-center">
                      <p className="text-2xl font-bold font-mono text-[#00B4D8]">{summary.total_materials}</p>
                      <p className="text-[10px] text-[#8892B0]">전체 소재</p>
                    </div>
                    <div className="p-4 rounded-lg bg-[#112240]/50 text-center">
                      <p className="text-2xl font-bold font-mono text-[#00E676]">{summary.verified}</p>
                      <p className="text-[10px] text-[#8892B0]">검증됨</p>
                    </div>
                    <div className="p-4 rounded-lg bg-[#112240]/50 text-center">
                      <p className="text-2xl font-bold font-mono text-[#7C4DFF]">{Object.keys(summary.by_category).length}</p>
                      <p className="text-[10px] text-[#8892B0]">카테고리</p>
                    </div>
                  </div>

                  {/* Source breakdown */}
                  <div>
                    <p className="text-xs text-[#8892B0] mb-2">데이터 출처</p>
                    <div className="space-y-2">
                      {Object.entries(summary.by_source).map(([src, cnt]) => (
                        <div key={src} className="flex items-center justify-between">
                          <div className="flex items-center gap-2">
                            <div className={`w-2 h-2 rounded-full ${src === "dongjin_internal" ? "bg-[#00B4D8]" : "bg-[#8892B0]/40"}`} />
                            <span className="text-xs text-[#8892B0]">{SOURCE_LABELS[src] || src}</span>
                          </div>
                          <span className="text-xs font-mono text-white">{cnt}</span>
                        </div>
                      ))}
                    </div>
                  </div>

                  {/* Property coverage */}
                  <div>
                    <p className="text-xs text-[#8892B0] mb-2">물성 데이터 커버리지</p>
                    <div className="space-y-2">
                      {Object.entries(summary.property_coverage).map(([prop, info]) => (
                        <div key={prop}>
                          <div className="flex justify-between mb-1">
                            <span className="text-[10px] text-[#8892B0]">{PROP_LABELS[prop] || prop}</span>
                            <span className="text-[10px] font-mono text-[#00B4D8]">{info.ratio}% ({info.count}건)</span>
                          </div>
                          <Progress value={info.ratio} className="h-1.5" />
                        </div>
                      ))}
                    </div>
                  </div>

                  {/* Category breakdown */}
                  <div>
                    <p className="text-xs text-[#8892B0] mb-2">카테고리별 분포</p>
                    <div className="flex flex-wrap gap-2">
                      {Object.entries(summary.by_category).map(([cat, cnt]) => (
                        <Badge key={cat} variant="outline" className="text-[10px] border-white/10 text-[#8892B0]">
                          {cat}: {cnt}
                        </Badge>
                      ))}
                    </div>
                  </div>
                </div>
              ) : (
                <div className="text-center py-10 text-[#8892B0] text-sm">
                  백엔드 API에 연결 중...
                </div>
              )}
            </CardContent>
          </Card>
        </div>
      </div>
    </div>
  );
}
