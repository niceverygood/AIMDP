-- =====================================================
-- Materials AI Platform — Supabase Database Setup
-- Run this in Supabase Dashboard → SQL Editor
-- =====================================================

-- 1. Materials table
CREATE TABLE IF NOT EXISTS materials (
    id BIGSERIAL PRIMARY KEY,
    name TEXT,
    smiles TEXT NOT NULL,
    category TEXT,
    molecular_weight REAL,
    logp REAL,
    hbd INTEGER,
    hba INTEGER,
    tpsa REAL,
    rotatable_bonds INTEGER,
    aromatic_rings INTEGER,
    thermal_stability REAL,
    dielectric_constant REAL,
    bandgap REAL,
    solubility REAL,
    density REAL,
    source TEXT DEFAULT 'upload',
    is_verified BOOLEAN DEFAULT false,
    created_at TIMESTAMPTZ DEFAULT NOW()
);

-- 2. Experiments table (feedback loop)
CREATE TABLE IF NOT EXISTS experiments (
    id BIGSERIAL PRIMARY KEY,
    material_id BIGINT REFERENCES materials(id),
    actual_thermal_stability REAL,
    actual_dielectric_constant REAL,
    actual_bandgap REAL,
    actual_solubility REAL,
    actual_density REAL,
    success BOOLEAN DEFAULT false,
    notes TEXT,
    researcher TEXT,
    created_at TIMESTAMPTZ DEFAULT NOW()
);

-- 3. Enable Row Level Security (allow all for now)
ALTER TABLE materials ENABLE ROW LEVEL SECURITY;
ALTER TABLE experiments ENABLE ROW LEVEL SECURITY;

CREATE POLICY "Allow all access to materials" ON materials FOR ALL USING (true) WITH CHECK (true);
CREATE POLICY "Allow all access to experiments" ON experiments FOR ALL USING (true) WITH CHECK (true);

-- 4. Indexes
CREATE INDEX IF NOT EXISTS idx_materials_category ON materials(category);
CREATE INDEX IF NOT EXISTS idx_materials_thermal ON materials(thermal_stability);
CREATE INDEX IF NOT EXISTS idx_materials_verified ON materials(is_verified);
