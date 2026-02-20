/**
 * Database access layer — Supabase를 통한 소재 데이터 CRUD.
 * Vercel에서 별도 백엔드 없이 직접 DB 접근.
 */

import { supabase, type MaterialRow } from "./supabase";

export async function getMaterialStats() {
  const { count: total } = await supabase
    .from("materials")
    .select("*", { count: "exact", head: true });

  const { count: verified } = await supabase
    .from("materials")
    .select("*", { count: "exact", head: true })
    .eq("is_verified", true);

  const { data: catData } = await supabase
    .from("materials")
    .select("category");

  const categories: Record<string, number> = {};
  catData?.forEach((row) => {
    const cat = row.category || "미분류";
    categories[cat] = (categories[cat] || 0) + 1;
  });

  return {
    total_materials: total || 0,
    verified_materials: verified || 0,
    categories,
  };
}

export async function searchMaterials(params: {
  category?: string;
  min_thermal?: number;
  max_thermal?: number;
  min_bandgap?: number;
  max_bandgap?: number;
  min_dielectric?: number;
  max_dielectric?: number;
  limit?: number;
}) {
  let query = supabase.from("materials").select("*");

  if (params.category && params.category !== "all") {
    query = query.ilike("category", params.category);
  }
  if (params.min_thermal != null) {
    query = query.gte("thermal_stability", params.min_thermal);
  }
  if (params.max_thermal != null) {
    query = query.lte("thermal_stability", params.max_thermal);
  }
  if (params.min_bandgap != null) {
    query = query.gte("bandgap", params.min_bandgap);
  }
  if (params.max_bandgap != null) {
    query = query.lte("bandgap", params.max_bandgap);
  }
  if (params.min_dielectric != null) {
    query = query.gte("dielectric_constant", params.min_dielectric);
  }
  if (params.max_dielectric != null) {
    query = query.lte("dielectric_constant", params.max_dielectric);
  }

  query = query.order("thermal_stability", { ascending: false, nullsFirst: false });
  query = query.limit(params.limit || 30);

  const { data, error, count } = await query;

  if (error) {
    console.error("Search error:", error);
    return { results: [], total: 0 };
  }

  // Calculate match scores
  const results = (data || []).map((row) => {
    const thermal_norm = row.thermal_stability ? row.thermal_stability / 500 : 0;
    const dielectric_norm = row.dielectric_constant ? Math.min(row.dielectric_constant / 25, 1) : 0;
    const bandgap_norm = row.bandgap ? row.bandgap / 8 : 0;
    const match_score = (thermal_norm * 0.33 + dielectric_norm * 0.33 + bandgap_norm * 0.34);

    return {
      ...row,
      similarity: 0,
      match_score: Math.round(match_score * 10000) / 10000,
    };
  });

  results.sort((a, b) => b.match_score - a.match_score);

  return { results, total: data?.length || 0 };
}

export async function getDataSummary() {
  const { count: total } = await supabase
    .from("materials")
    .select("*", { count: "exact", head: true });

  const { count: verified } = await supabase
    .from("materials")
    .select("*", { count: "exact", head: true })
    .eq("is_verified", true);

  const { data: allData } = await supabase
    .from("materials")
    .select("source, category, thermal_stability, dielectric_constant, bandgap, solubility, density");

  const by_source: Record<string, number> = {};
  const by_category: Record<string, number> = {};
  const property_counts: Record<string, number> = {
    thermal_stability: 0,
    dielectric_constant: 0,
    bandgap: 0,
    solubility: 0,
    density: 0,
  };

  allData?.forEach((row) => {
    const src = row.source || "unknown";
    by_source[src] = (by_source[src] || 0) + 1;
    const cat = row.category || "미분류";
    by_category[cat] = (by_category[cat] || 0) + 1;
    if (row.thermal_stability != null) property_counts.thermal_stability++;
    if (row.dielectric_constant != null) property_counts.dielectric_constant++;
    if (row.bandgap != null) property_counts.bandgap++;
    if (row.solubility != null) property_counts.solubility++;
    if (row.density != null) property_counts.density++;
  });

  const t = total || 1;
  const property_coverage: Record<string, { count: number; ratio: number }> = {};
  for (const [k, v] of Object.entries(property_counts)) {
    property_coverage[k] = { count: v, ratio: Math.round((v / t) * 1000) / 10 };
  }

  return {
    total_materials: total || 0,
    verified: verified || 0,
    by_source,
    by_category,
    property_coverage,
  };
}

export async function insertMaterials(rows: Partial<MaterialRow>[]) {
  const { data, error } = await supabase.from("materials").insert(rows).select();
  if (error) throw error;
  return data;
}
