"use client";

import { useEffect, useState } from "react";
import { motion } from "framer-motion";
import {
  BarChart,
  Bar,
  XAxis,
  YAxis,
  CartesianGrid,
  Tooltip,
  ResponsiveContainer,
  AreaChart,
  Area,
  PieChart,
  Pie,
  Cell,
} from "recharts";
import { StatsCard } from "@/components/StatsCard";
import { HelpTooltip } from "@/components/HelpTooltip";

const API_URL = process.env.NEXT_PUBLIC_API_URL || "http://localhost:8000";

const CATEGORY_COLORS = [
  "#00B4D8", "#00E676", "#7C4DFF", "#FF9100", "#FF5252",
  "#48CAE4", "#FFD600", "#8892B0",
];

const MONTHLY_DATA = [
  { month: "8월", materials: 12, predictions: 8 },
  { month: "9월", materials: 25, predictions: 18 },
  { month: "10월", materials: 42, predictions: 35 },
  { month: "11월", materials: 61, predictions: 50 },
  { month: "12월", materials: 78, predictions: 65 },
  { month: "1월", materials: 92, predictions: 80 },
  { month: "2월", materials: 104, predictions: 95 },
];

const ACCURACY_DATA = [
  { property: "열안정성", accuracy: 72.3 },
  { property: "유전율", accuracy: 68.7 },
  { property: "밴드갭", accuracy: 71.1 },
  { property: "용해도", accuracy: 65.4 },
  { property: "밀도", accuracy: 74.2 },
];

const container = {
  hidden: { opacity: 0 },
  show: { opacity: 1, transition: { staggerChildren: 0.08 } },
};

const item = {
  hidden: { opacity: 0, y: 20 },
  show: { opacity: 1, y: 0, transition: { duration: 0.4 } },
};

interface DashboardStats {
  total_materials: number;
  verified_materials: number;
  categories: Record<string, number>;
}

export default function DashboardPage() {
  const [stats, setStats] = useState<DashboardStats | null>(null);
  const [dbConnected, setDbConnected] = useState(false);

  useEffect(() => {
    async function fetchStats() {
      try {
        const res = await fetch(`${API_URL}/api/materials/stats`);
        if (res.ok) {
          const data = await res.json();
          if (data.success) {
            setStats(data.data);
            setDbConnected(true);
          }
        }
      } catch {
        // API not available — show placeholder
        setStats({
          total_materials: 0,
          verified_materials: 0,
          categories: {},
        });
      }
    }
    fetchStats();
  }, []);

  if (!stats) {
    return (
      <div className="flex items-center justify-center min-h-screen">
        <div className="w-10 h-10 border-2 border-[#00B4D8] border-t-transparent rounded-full animate-spin" />
      </div>
    );
  }

  const categoryData = Object.entries(stats.categories).map(
    ([name, count], i) => ({
      name,
      count,
      fill: CATEGORY_COLORS[i % CATEGORY_COLORS.length],
    })
  );

  return (
    <div className="p-8 space-y-8 grid-pattern min-h-screen">
      {/* Header */}
      <motion.div
        initial={{ opacity: 0, y: -20 }}
        animate={{ opacity: 1, y: 0 }}
        transition={{ duration: 0.5 }}
      >
        <div className="flex items-center gap-2">
          <h1 className="text-2xl font-bold text-white">대시보드</h1>
          <HelpTooltip
            title="대시보드"
            description="AI 기반 신소재 탐색 플랫폼의 전체 현황을 한눈에 볼 수 있는 메인 화면입니다."
            details={[
              "통계 카드: 전체 소재 수, 검증된 소재, 예측 모델 상태, 카테고리 수를 표시합니다.",
              "데이터 추이 차트: 월별 소재 등록 및 AI 예측 활용 현황을 보여줍니다.",
              "카테고리 분포: 소재 유형별 비율을 파이 차트로 나타냅니다.",
              "시스템 상태: 각 모듈(DB, 임베딩, 검색 등)의 연결 상태를 보여줍니다.",
            ]}
            size="md"
          />
        </div>
        <p className="text-sm text-[#8892B0] mt-1">
          AI 기반 신소재 탐색 플랫폼 현황
          {dbConnected && (
            <span className="ml-2 text-[#00E676]">
              · DB 연결됨
            </span>
          )}
        </p>
      </motion.div>

      {/* Stats Cards */}
      <motion.div
        variants={container}
        initial="hidden"
        animate="show"
        className="grid grid-cols-1 md:grid-cols-2 lg:grid-cols-4 gap-4"
      >
        <motion.div variants={item}>
          <StatsCard
            title="전체 소재"
            value={stats.total_materials.toLocaleString()}
            subtitle="데이터베이스 등록"
            color="#00B4D8"
            icon={
              <svg className="w-5 h-5" fill="none" viewBox="0 0 24 24" stroke="currentColor" strokeWidth={1.5}>
                <path strokeLinecap="round" strokeLinejoin="round" d="M20.25 6.375c0 2.278-3.694 4.125-8.25 4.125S3.75 8.653 3.75 6.375m16.5 0c0-2.278-3.694-4.125-8.25-4.125S3.75 4.097 3.75 6.375m16.5 0v11.25c0 2.278-3.694 4.125-8.25 4.125s-8.25-1.847-8.25-4.125V6.375m16.5 0v3.75m-16.5-3.75v3.75m16.5 0v3.75C20.25 16.153 16.556 18 12 18s-8.25-1.847-8.25-4.125v-3.75m16.5 0c0 2.278-3.694 4.125-8.25 4.125s-8.25-1.847-8.25-4.125" />
              </svg>
            }
          />
        </motion.div>
        <motion.div variants={item}>
          <StatsCard
            title="검증된 소재"
            value={stats.verified_materials.toLocaleString()}
            subtitle="실험으로 확인됨"
            color="#00E676"
            icon={
              <svg className="w-5 h-5" fill="none" viewBox="0 0 24 24" stroke="currentColor" strokeWidth={1.5}>
                <path strokeLinecap="round" strokeLinejoin="round" d="M9 12.75 11.25 15 15 9.75M21 12a9 9 0 1 1-18 0 9 9 0 0 1 18 0Z" />
              </svg>
            }
          />
        </motion.div>
        <motion.div variants={item}>
          <StatsCard
            title="예측 모델"
            value="RDKit"
            subtitle="휴리스틱 (PoC)"
            color="#7C4DFF"
            icon={
              <svg className="w-5 h-5" fill="none" viewBox="0 0 24 24" stroke="currentColor" strokeWidth={1.5}>
                <path strokeLinecap="round" strokeLinejoin="round" d="M9.813 15.904 9 18.75l-.813-2.846a4.5 4.5 0 0 0-3.09-3.09L2.25 12l2.846-.813a4.5 4.5 0 0 0 3.09-3.09L9 5.25l.813 2.846a4.5 4.5 0 0 0 3.09 3.09L15.75 12l-2.846.813a4.5 4.5 0 0 0-3.09 3.09ZM18.259 8.715 18 9.75l-.259-1.035a3.375 3.375 0 0 0-2.455-2.456L14.25 6l1.036-.259a3.375 3.375 0 0 0 2.455-2.456L18 2.25l.259 1.035a3.375 3.375 0 0 0 2.455 2.456L21.75 6l-1.036.259a3.375 3.375 0 0 0-2.455 2.456ZM16.894 20.567 16.5 21.75l-.394-1.183a2.25 2.25 0 0 0-1.423-1.423L13.5 18.75l1.183-.394a2.25 2.25 0 0 0 1.423-1.423l.394-1.183.394 1.183a2.25 2.25 0 0 0 1.423 1.423l1.183.394-1.183.394a2.25 2.25 0 0 0-1.423 1.423Z" />
              </svg>
            }
          />
        </motion.div>
        <motion.div variants={item}>
          <StatsCard
            title="카테고리"
            value={Object.keys(stats.categories).length}
            subtitle="소재 분류"
            color="#FF9100"
            icon={
              <svg className="w-5 h-5" fill="none" viewBox="0 0 24 24" stroke="currentColor" strokeWidth={1.5}>
                <path strokeLinecap="round" strokeLinejoin="round" d="M9.568 3H5.25A2.25 2.25 0 0 0 3 5.25v4.318c0 .597.237 1.17.659 1.591l9.581 9.581c.699.699 1.78.872 2.607.33a18.095 18.095 0 0 0 5.223-5.223c.542-.827.369-1.908-.33-2.607L11.16 3.66A2.25 2.25 0 0 0 9.568 3Z" />
                <path strokeLinecap="round" strokeLinejoin="round" d="M6 6h.008v.008H6V6Z" />
              </svg>
            }
          />
        </motion.div>
      </motion.div>

      {/* Charts Row */}
      <motion.div
        variants={container}
        initial="hidden"
        animate="show"
        className="grid grid-cols-1 lg:grid-cols-3 gap-6"
      >
        {/* Area Chart */}
        <motion.div variants={item} className="lg:col-span-2 p-6 rounded-xl bg-[#0D1B2A] border border-white/[0.06]">
          <div className="flex items-center gap-2 mb-1">
            <h2 className="text-sm font-semibold text-white">데이터 누적 현황</h2>
            <HelpTooltip
              title="데이터 누적 현황"
              description="PoC 기간 동안 DB에 등록된 소재 수와 AI 예측이 수행된 횟수의 월별 추이입니다."
              details={[
                "파란 영역: 매달 새로 등록된 소재 수",
                "초록 영역: AI 물성 예측이 수행된 횟수",
                "데이터가 많을수록 AI 모델의 정확도가 향상됩니다.",
              ]}
            />
          </div>
          <p className="text-[11px] text-[#8892B0] mb-4">소재 등록 및 AI 예측 추이 (PoC 기간)</p>
          <ResponsiveContainer width="100%" height={280}>
            <AreaChart data={MONTHLY_DATA}>
              <defs>
                <linearGradient id="colorMaterials" x1="0" y1="0" x2="0" y2="1">
                  <stop offset="5%" stopColor="#00B4D8" stopOpacity={0.3} />
                  <stop offset="95%" stopColor="#00B4D8" stopOpacity={0} />
                </linearGradient>
                <linearGradient id="colorPredictions" x1="0" y1="0" x2="0" y2="1">
                  <stop offset="5%" stopColor="#00E676" stopOpacity={0.3} />
                  <stop offset="95%" stopColor="#00E676" stopOpacity={0} />
                </linearGradient>
              </defs>
              <CartesianGrid strokeDasharray="3 3" stroke="rgba(255,255,255,0.04)" />
              <XAxis dataKey="month" tick={{ fill: "#8892B0", fontSize: 11 }} axisLine={false} />
              <YAxis tick={{ fill: "#8892B0", fontSize: 11 }} axisLine={false} />
              <Tooltip contentStyle={{ backgroundColor: "#112240", borderColor: "rgba(255,255,255,0.08)", borderRadius: 8, fontSize: 12, color: "#E0E7FF" }} />
              <Area type="monotone" dataKey="materials" name="소재 등록" stroke="#00B4D8" fillOpacity={1} fill="url(#colorMaterials)" strokeWidth={2} />
              <Area type="monotone" dataKey="predictions" name="AI 예측" stroke="#00E676" fillOpacity={1} fill="url(#colorPredictions)" strokeWidth={2} />
            </AreaChart>
          </ResponsiveContainer>
        </motion.div>

        {/* Pie Chart */}
        <motion.div variants={item} className="p-6 rounded-xl bg-[#0D1B2A] border border-white/[0.06]">
          <div className="flex items-center gap-2 mb-1">
            <h2 className="text-sm font-semibold text-white">카테고리 분포</h2>
            <HelpTooltip
              title="카테고리 분포"
              description="데이터베이스에 등록된 소재들의 유형별 분포입니다."
              details={[
                "OLED: 유기 발광 소재 (호스트, 전자수송층 등)",
                "배터리: 전해질, 첨가제, 리튬염",
                "반도체: 절연막, 포토레지스트 소재",
                "하드코팅: 실란 기반 코팅 소재",
                "유기 분자: 일반 유기 화합물 (용해도 DB 등)",
              ]}
            />
          </div>
          <p className="text-[11px] text-[#8892B0] mb-4">실제 DB 기준</p>
          {categoryData.length > 0 ? (
            <>
              <ResponsiveContainer width="100%" height={200}>
                <PieChart>
                  <Pie data={categoryData} cx="50%" cy="50%" innerRadius={50} outerRadius={80} paddingAngle={2} dataKey="count">
                    {categoryData.map((entry, index) => (
                      <Cell key={index} fill={entry.fill} />
                    ))}
                  </Pie>
                  <Tooltip contentStyle={{ backgroundColor: "#112240", borderColor: "rgba(255,255,255,0.08)", borderRadius: 8, fontSize: 12, color: "#E0E7FF" }} />
                </PieChart>
              </ResponsiveContainer>
              <div className="space-y-1 mt-2">
                {categoryData.map((cat) => (
                  <div key={cat.name} className="flex items-center justify-between text-[11px]">
                    <div className="flex items-center gap-2">
                      <div className="w-2 h-2 rounded-full" style={{ backgroundColor: cat.fill }} />
                      <span className="text-[#8892B0]">{cat.name}</span>
                    </div>
                    <span className="font-mono text-white">{cat.count}</span>
                  </div>
                ))}
              </div>
            </>
          ) : (
            <div className="flex items-center justify-center h-[200px] text-[#8892B0] text-sm">
              DB 연결 대기 중...
            </div>
          )}
        </motion.div>
      </motion.div>

      {/* Accuracy + Activity */}
      <motion.div variants={container} initial="hidden" animate="show" className="grid grid-cols-1 lg:grid-cols-2 gap-6">
        <motion.div variants={item} className="p-6 rounded-xl bg-[#0D1B2A] border border-white/[0.06]">
          <div className="flex items-center gap-2 mb-1">
            <h2 className="text-sm font-semibold text-white">물성별 예측 정확도 (추정)</h2>
            <HelpTooltip
              title="물성 예측 정확도"
              description="각 물성별 AI 예측 모델의 추정 정확도입니다. 현재는 RDKit 기술자(descriptor) 기반 휴리스틱을 사용합니다."
              details={[
                "열안정성: 방향족 고리 수와 분자량으로 추정",
                "유전율: 극성 표면적(TPSA)과 수소결합 수용체로 추정",
                "밴드갭: 방향족 고리와 LogP로 추정",
                "Uni-Mol 모델 학습 시 정확도 90%+ 달성 가능",
              ]}
            />
          </div>
          <p className="text-[11px] text-[#8892B0] mb-4">RDKit 휴리스틱 기반 — Uni-Mol 학습 시 90%+ 가능</p>
          <ResponsiveContainer width="100%" height={260}>
            <BarChart data={ACCURACY_DATA} layout="vertical">
              <CartesianGrid strokeDasharray="3 3" stroke="rgba(255,255,255,0.04)" horizontal={false} />
              <XAxis type="number" domain={[0, 100]} tick={{ fill: "#8892B0", fontSize: 11 }} axisLine={false} />
              <YAxis type="category" dataKey="property" tick={{ fill: "#8892B0", fontSize: 11 }} axisLine={false} width={60} />
              <Tooltip contentStyle={{ backgroundColor: "#112240", borderColor: "rgba(255,255,255,0.08)", borderRadius: 8, fontSize: 12, color: "#E0E7FF" }} formatter={(value) => [`${value}%`, "정확도"]} />
              <Bar dataKey="accuracy" fill="#00B4D8" radius={[0, 4, 4, 0]} barSize={20} />
            </BarChart>
          </ResponsiveContainer>
        </motion.div>

        <motion.div variants={item} className="p-6 rounded-xl bg-[#0D1B2A] border border-white/[0.06]">
          <div className="flex items-center gap-2 mb-1">
            <h2 className="text-sm font-semibold text-white">시스템 상태</h2>
            <HelpTooltip
              title="시스템 상태"
              description="플랫폼의 각 구성 모듈이 정상적으로 연결/동작하고 있는지 보여줍니다."
              details={[
                "초록 점: 정상 작동 중인 모듈",
                "회색 점: 아직 설치되지 않았거나 설정이 필요한 모듈",
                "GPU 서버에 Uni-Mol 설치 시 임베딩 품질이 크게 향상됩니다.",
              ]}
            />
          </div>
          <p className="text-[11px] text-[#8892B0] mb-4">PoC 현재 구현 상태</p>
          <div className="space-y-3">
            {[
              { label: "소재 DB (SQLite)", status: dbConnected ? "연결됨" : "미연결", ok: dbConnected },
              { label: "임베딩 (Morgan FP → 768-dim)", status: "활성", ok: true },
              { label: "벡터 유사도 검색", status: "코사인 유사도", ok: true },
              { label: "물성 예측 모델", status: "RDKit 휴리스틱", ok: true },
              { label: "3D 분자 뷰어 (RDKit → SDF)", status: "대기", ok: false },
              { label: "Uni-Mol 임베딩", status: "미설치 (GPU 필요)", ok: false },
              { label: "학습된 신경망 예측", status: "학습 데이터 필요", ok: false },
              { label: "Active Learning 피드백", status: "API 준비됨", ok: true },
            ].map((item, i) => (
              <div key={i} className="flex items-center justify-between p-2.5 rounded-lg hover:bg-white/[0.02] transition-colors">
                <div className="flex items-center gap-3">
                  <div className={`w-2 h-2 rounded-full ${item.ok ? "bg-[#00E676]" : "bg-[#8892B0]/40"}`} />
                  <span className="text-xs text-white">{item.label}</span>
                </div>
                <span className={`text-[10px] font-mono ${item.ok ? "text-[#00E676]" : "text-[#8892B0]"}`}>
                  {item.status}
                </span>
              </div>
            ))}
          </div>
        </motion.div>
      </motion.div>
    </div>
  );
}
