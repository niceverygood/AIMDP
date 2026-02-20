/**
 * POST /api/setup-feedback
 * Supabase에 researcher_feedback 테이블이 없으면 생성하는 안내를 반환합니다.
 * 테이블은 Supabase Dashboard → SQL Editor에서 직접 생성해야 합니다.
 */
import { NextResponse } from "next/server";

export async function GET() {
  return NextResponse.json({
    message: "Supabase SQL Editor에서 아래 SQL을 실행하세요",
    sql: `
CREATE TABLE IF NOT EXISTS researcher_feedback (
    id BIGSERIAL PRIMARY KEY,
    material_name TEXT,
    material_smiles TEXT,
    category TEXT,
    rating INTEGER CHECK (rating >= 1 AND rating <= 5),
    experiment_success BOOLEAN,
    actual_thermal_stability REAL,
    actual_dielectric_constant REAL,
    actual_bandgap REAL,
    actual_solubility REAL,
    actual_density REAL,
    pros TEXT,
    cons TEXT,
    notes TEXT,
    researcher TEXT,
    feedback_type TEXT DEFAULT 'evaluation',
    created_at TIMESTAMPTZ DEFAULT NOW()
);

ALTER TABLE researcher_feedback ENABLE ROW LEVEL SECURITY;
CREATE POLICY "Allow all access to researcher_feedback" ON researcher_feedback FOR ALL USING (true) WITH CHECK (true);
CREATE INDEX idx_feedback_category ON researcher_feedback(category);
CREATE INDEX idx_feedback_rating ON researcher_feedback(rating);
    `,
  });
}
