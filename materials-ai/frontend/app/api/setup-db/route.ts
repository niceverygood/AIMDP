import { NextResponse } from "next/server";
import { supabase } from "@/lib/supabase";

const SEED_DATA = [
  // OLED (14)
  { name: "BPhen", smiles: "c1ccc(-c2ccnc3c2ccc2c(-c4ccccc4)ccnc23)cc1", category: "OLED", thermal_stability: 380, dielectric_constant: 3.2, bandgap: 3.5, solubility: 0.002, density: 1.25, molecular_weight: 332.41, source: "curated", is_verified: true },
  { name: "TPBi", smiles: "c1ccc(-n2c(-c3cc(-c4nc5ccccc5[nH]4)cc(-c4nc5ccccc5[nH]4)c3)nc3ccccc32)cc1", category: "OLED", thermal_stability: 420, dielectric_constant: 3.0, bandgap: 3.2, solubility: 0.001, density: 1.32, molecular_weight: 502.58, source: "curated", is_verified: true },
  { name: "Ir(ppy)3", smiles: "c1ccc(-c2ccccn2)cc1", category: "OLED", thermal_stability: 450, dielectric_constant: 3.3, bandgap: 2.4, solubility: 0.0003, density: 1.68, molecular_weight: 155.2, source: "curated", is_verified: true },
  { name: "TCTA", smiles: "c1ccc(-n2c3ccccc3c3ccccc32)cc1", category: "OLED", thermal_stability: 395, dielectric_constant: 3.1, bandgap: 3.4, solubility: 0.001, density: 1.30, molecular_weight: 243.31, source: "curated", is_verified: true },
  { name: "CBP", smiles: "c1ccc2c(c1)[nH]c1ccccc12", category: "OLED", thermal_stability: 365, dielectric_constant: 2.9, bandgap: 3.6, solubility: 0.001, density: 1.18, molecular_weight: 168.19, source: "curated", is_verified: true },
  { name: "NPB", smiles: "c1ccc(N(c2ccccc2)c2ccc(-c3ccc(N(c4ccccc4)c4ccc5ccccc5c4)cc3)cc2)cc1", category: "OLED", thermal_stability: 380, dielectric_constant: 3.0, bandgap: 3.1, solubility: 0.001, density: 1.23, molecular_weight: 588.74, source: "curated", is_verified: true },
  { name: "DPEPO", smiles: "O=P(c1ccccc1)(c1ccccc1)c1ccccc1", category: "OLED", thermal_stability: 350, dielectric_constant: 3.8, bandgap: 4.1, solubility: 0.003, density: 1.28, molecular_weight: 278.29, source: "curated", is_verified: true },
  { name: "TAPC", smiles: "CC(c1ccc(N(c2ccccc2)c2ccccc2)cc1)c1ccc(N(c2ccccc2)c2ccccc2)cc1", category: "OLED", thermal_stability: 360, dielectric_constant: 2.8, bandgap: 3.5, solubility: 0.001, density: 1.15, molecular_weight: 466.62, source: "curated", is_verified: true },
  { name: "PhCz", smiles: "c1ccc(-c2ccc(-n3c4ccccc4c4ccccc43)cc2)cc1", category: "OLED", thermal_stability: 355, dielectric_constant: 2.85, bandgap: 3.55, solubility: 0.0015, density: 1.20, molecular_weight: 319.40, source: "curated", is_verified: true },
  { name: "Fluoranthene", smiles: "c1ccc2c(c1)c1ccccc1c1ccccc12", category: "OLED", thermal_stability: 375, dielectric_constant: 2.7, bandgap: 3.3, solubility: 0.002, density: 1.27, molecular_weight: 202.25, source: "curated", is_verified: true },
  { name: "2-PhQ", smiles: "c1ccc(-c2ccc3ccccc3n2)cc1", category: "OLED", thermal_stability: 340, dielectric_constant: 3.15, bandgap: 3.0, solubility: 0.003, density: 1.22, molecular_weight: 205.26, source: "curated", is_verified: true },
  { name: "BBI", smiles: "c1ccc2[nH]c(-c3ccc(-c4nc5ccccc5[nH]4)cc3)nc2c1", category: "OLED", thermal_stability: 410, dielectric_constant: 3.05, bandgap: 3.25, solubility: 0.0012, density: 1.35, molecular_weight: 310.36, source: "curated", is_verified: true },
  { name: "Fluorene", smiles: "c1ccc2c(c1)-c1ccccc1-2", category: "OLED", thermal_stability: 310, dielectric_constant: 2.75, bandgap: 3.8, solubility: 0.003, density: 1.10, molecular_weight: 166.22, source: "curated", is_verified: false },
  { name: "Dibenzofuran", smiles: "c1ccc2oc3ccccc3c2c1", category: "OLED", thermal_stability: 330, dielectric_constant: 2.8, bandgap: 3.9, solubility: 0.004, density: 1.17, molecular_weight: 168.19, source: "curated", is_verified: false },
  // Battery (12)
  { name: "EC", smiles: "O=C1OCCO1", category: "battery", thermal_stability: 248, dielectric_constant: 89.8, bandgap: 6.7, solubility: 0.9, density: 1.32, molecular_weight: 88.06, source: "curated", is_verified: true },
  { name: "DMC", smiles: "COC(=O)OC", category: "battery", thermal_stability: 90, dielectric_constant: 3.1, bandgap: 6.5, solubility: 0.8, density: 1.07, molecular_weight: 90.08, source: "curated", is_verified: true },
  { name: "EMC", smiles: "CCOC(=O)OC", category: "battery", thermal_stability: 107, dielectric_constant: 2.9, bandgap: 6.4, solubility: 0.7, density: 1.01, molecular_weight: 104.1, source: "curated", is_verified: true },
  { name: "FEC", smiles: "O=C1OCC(F)O1", category: "battery", thermal_stability: 212, dielectric_constant: 78.4, bandgap: 5.9, solubility: 0.85, density: 1.45, molecular_weight: 106.05, source: "curated", is_verified: true },
  { name: "VC", smiles: "O=c1occo1", category: "battery", thermal_stability: 162, dielectric_constant: 126.0, bandgap: 5.5, solubility: 0.9, density: 1.36, molecular_weight: 86.05, source: "curated", is_verified: true },
  { name: "LiTFSI", smiles: "O=S(=O)(C(F)(F)F)N([Li])S(=O)(=O)C(F)(F)F", category: "battery", thermal_stability: 360, dielectric_constant: null, bandgap: null, solubility: 0.95, density: 1.33, molecular_weight: 287.09, source: "curated", is_verified: true },
  { name: "PC", smiles: "CC1COC(=O)O1", category: "battery", thermal_stability: 240, dielectric_constant: 64.9, bandgap: 6.6, solubility: 0.88, density: 1.20, molecular_weight: 102.09, source: "curated", is_verified: true },
  { name: "DEC", smiles: "CCOC(=O)OCC", category: "battery", thermal_stability: 126, dielectric_constant: 2.8, bandgap: 6.3, solubility: 0.6, density: 0.98, molecular_weight: 118.13, source: "curated", is_verified: true },
  { name: "DMS", smiles: "CS(=O)(=O)C", category: "battery", thermal_stability: 285, dielectric_constant: 47.2, bandgap: 5.8, solubility: 0.75, density: 1.13, molecular_weight: 94.13, source: "curated", is_verified: true },
  { name: "Sulfolane", smiles: "O=S1(=O)CCCC1", category: "battery", thermal_stability: 287, dielectric_constant: 43.4, bandgap: 5.6, solubility: 0.78, density: 1.26, molecular_weight: 120.17, source: "curated", is_verified: true },
  { name: "TEP", smiles: "CCOP(=O)(OCC)OCC", category: "battery", thermal_stability: 215, dielectric_constant: 13.2, bandgap: 6.8, solubility: 0.50, density: 1.07, molecular_weight: 182.15, source: "curated", is_verified: true },
  { name: "GBL", smiles: "O=C1OCCC1", category: "battery", thermal_stability: 204, dielectric_constant: 39.1, bandgap: 6.2, solubility: 0.85, density: 1.13, molecular_weight: 86.09, source: "curated", is_verified: true },
  // Semiconductor & Coating (10)
  { name: "PMDA", smiles: "O=c1oc(=O)c2cc3c(=O)oc(=O)c3cc12", category: "semiconductor", thermal_stability: 480, dielectric_constant: 3.4, bandgap: 3.0, solubility: 0.0001, density: 1.42, molecular_weight: 218.12, source: "curated", is_verified: true },
  { name: "BCB", smiles: "C1=CC2=CC=CC2=C1", category: "semiconductor", thermal_stability: 350, dielectric_constant: 2.65, bandgap: 4.0, solubility: 0.005, density: 1.05, molecular_weight: 104.15, source: "curated", is_verified: true },
  { name: "TEOS", smiles: "CCO[Si](OCC)(OCC)OCC", category: "hard_coating", thermal_stability: 165, dielectric_constant: 3.9, bandgap: 8.9, solubility: 0.01, density: 0.93, molecular_weight: 208.33, source: "curated", is_verified: true },
  { name: "GPTMS", smiles: "CO[Si](OC)(OC)CCCOCC1CO1", category: "hard_coating", thermal_stability: 200, dielectric_constant: 3.5, bandgap: 7.0, solubility: 0.02, density: 1.07, molecular_weight: 236.34, source: "curated", is_verified: true },
  { name: "TMSPM", smiles: "CO[Si](OC)(OC)CCCOC(=O)C(=C)C", category: "hard_coating", thermal_stability: 220, dielectric_constant: 3.2, bandgap: 6.5, solubility: 0.01, density: 1.04, molecular_weight: 248.35, source: "curated", is_verified: true },
  { name: "5CB", smiles: "CCCCCC1=CC=C(C#N)C=C1", category: "display", thermal_stability: 168, dielectric_constant: 11.0, bandgap: 4.2, solubility: 0.001, density: 1.02, molecular_weight: 249.36, source: "curated", is_verified: true },
  { name: "7CB", smiles: "CCCCCCCC1=CC=C(C#N)C=C1", category: "display", thermal_stability: 185, dielectric_constant: 12.5, bandgap: 4.1, solubility: 0.001, density: 1.01, molecular_weight: 277.41, source: "curated", is_verified: true },
  // Common organics (10)
  { name: "Benzene", smiles: "c1ccccc1", category: "organic", thermal_stability: 172, dielectric_constant: 2.3, bandgap: 4.0, solubility: 0.6, density: 0.88, molecular_weight: 78.11, source: "esol", is_verified: false },
  { name: "Naphthalene", smiles: "c1ccc2ccccc2c1", category: "organic", thermal_stability: 250, dielectric_constant: 2.5, bandgap: 3.6, solubility: 0.35, density: 1.14, molecular_weight: 128.17, source: "esol", is_verified: false },
  { name: "Caffeine", smiles: "Cn1c(=O)c2c(ncn2C)n(C)c1=O", category: "organic", thermal_stability: 236, dielectric_constant: 4.5, bandgap: 4.8, solubility: 0.64, density: 1.23, molecular_weight: 194.19, source: "esol", is_verified: false },
  { name: "Aspirin", smiles: "CC(=O)Oc1ccccc1C(O)=O", category: "organic", thermal_stability: 210, dielectric_constant: 3.5, bandgap: 4.2, solubility: 0.4, density: 1.40, molecular_weight: 180.16, source: "esol", is_verified: false },
  { name: "Carbazole", smiles: "c1ccc2c(c1)[nH]c1ccccc12", category: "organic", thermal_stability: 340, dielectric_constant: 2.9, bandgap: 3.6, solubility: 0.15, density: 1.30, molecular_weight: 167.21, source: "esol", is_verified: false },
  { name: "Phenol", smiles: "Oc1ccccc1", category: "organic", thermal_stability: 182, dielectric_constant: 4.5, bandgap: 4.3, solubility: 0.75, density: 1.07, molecular_weight: 94.11, source: "esol", is_verified: false },
  { name: "Aniline", smiles: "Nc1ccccc1", category: "organic", thermal_stability: 184, dielectric_constant: 3.8, bandgap: 4.1, solubility: 0.68, density: 1.02, molecular_weight: 93.13, source: "esol", is_verified: false },
  { name: "Pyridine", smiles: "c1ccncc1", category: "organic", thermal_stability: 170, dielectric_constant: 4.2, bandgap: 4.5, solubility: 0.77, density: 0.98, molecular_weight: 79.10, source: "esol", is_verified: false },
  { name: "Indole", smiles: "c1ccc2[nH]ccc2c1", category: "organic", thermal_stability: 253, dielectric_constant: 3.2, bandgap: 3.7, solubility: 0.35, density: 1.17, molecular_weight: 117.15, source: "esol", is_verified: false },
  { name: "Quinoline", smiles: "c1ccc2ncccc2c1", category: "organic", thermal_stability: 238, dielectric_constant: 3.5, bandgap: 3.8, solubility: 0.52, density: 1.09, molecular_weight: 129.16, source: "esol", is_verified: false },
];

export async function POST() {
  // Check if table exists by trying a select
  const { error: checkError } = await supabase.from("materials").select("id").limit(1);

  if (checkError?.code === "PGRST205") {
    return NextResponse.json({
      success: false,
      error: "테이블이 없습니다. Supabase Dashboard → SQL Editor에서 supabase_setup.sql을 실행하세요.",
      sql_url: "https://supabase.com/dashboard/project/hkzefauqitxychpnjkra/sql",
    });
  }

  // Check current count
  const { count } = await supabase.from("materials").select("*", { count: "exact", head: true });

  if (count && count > 0) {
    return NextResponse.json({
      success: true,
      message: `이미 ${count}건의 데이터가 있습니다. 시드를 건너뜁니다.`,
      count,
    });
  }

  // Insert seed data
  const { data, error } = await supabase.from("materials").insert(SEED_DATA).select();

  if (error) {
    return NextResponse.json({ success: false, error: error.message });
  }

  return NextResponse.json({
    success: true,
    message: `${data.length}건의 시드 데이터가 등록되었습니다.`,
    count: data.length,
  });
}
