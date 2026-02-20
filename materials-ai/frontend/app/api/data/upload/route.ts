import { NextRequest, NextResponse } from "next/server";
import { supabase } from "@/lib/supabase";
import { parse } from "csv-parse/sync";

const COL_MAP: Record<string, string> = {
  smiles: "smiles", smi: "smiles", SMILES: "smiles",
  name: "name", "이름": "name", "소재명": "name", material: "name",
  category: "category", "카테고리": "category", "분류": "category",
  thermal_stability: "thermal_stability", td: "thermal_stability", "열안정성": "thermal_stability",
  dielectric_constant: "dielectric_constant", dk: "dielectric_constant", "유전율": "dielectric_constant",
  bandgap: "bandgap", eg: "bandgap", "밴드갭": "bandgap", band_gap: "bandgap",
  solubility: "solubility", "용해도": "solubility",
  density: "density", "밀도": "density",
  molecular_weight: "molecular_weight", "분자량": "molecular_weight",
};

export async function POST(request: NextRequest) {
  try {
    const formData = await request.formData();
    const file = formData.get("file") as File;
    const category = formData.get("category") as string || "";
    const source = formData.get("source") as string || "dongjin_internal";

    if (!file) {
      return NextResponse.json({ success: false, error: "파일이 없습니다" }, { status: 400 });
    }

    const text = await file.text();
    let records: Record<string, string>[];

    try {
      records = parse(text, {
        columns: true,
        skip_empty_lines: true,
        trim: true,
        bom: true,
      });
    } catch {
      // Try tab-separated
      records = parse(text, {
        columns: true,
        skip_empty_lines: true,
        trim: true,
        delimiter: "\t",
        bom: true,
      });
    }

    if (!records.length) {
      return NextResponse.json({ success: false, error: "파일에 데이터가 없습니다" }, { status: 400 });
    }

    // Map columns
    const mappedRecords = records.map((row) => {
      const mapped: Record<string, unknown> = {};
      for (const [key, value] of Object.entries(row)) {
        const canonical = COL_MAP[key.trim().toLowerCase()] || COL_MAP[key.trim()];
        if (canonical) {
          mapped[canonical] = value;
        }
      }
      return mapped;
    });

    // Validate SMILES exists
    const withSmiles = mappedRecords.filter((r) => r.smiles && String(r.smiles).trim());

    if (!withSmiles.length) {
      return NextResponse.json({
        success: false,
        error: "SMILES 컬럼을 찾을 수 없습니다. 'SMILES' 또는 'smiles' 컬럼이 필요합니다.",
      }, { status: 400 });
    }

    // Build insert rows
    const insertRows = withSmiles.map((r) => ({
      name: r.name ? String(r.name) : null,
      smiles: String(r.smiles).trim(),
      category: category && category !== "auto" ? category : (r.category ? String(r.category) : null),
      thermal_stability: r.thermal_stability ? parseFloat(String(r.thermal_stability)) || null : null,
      dielectric_constant: r.dielectric_constant ? parseFloat(String(r.dielectric_constant)) || null : null,
      bandgap: r.bandgap ? parseFloat(String(r.bandgap)) || null : null,
      solubility: r.solubility ? parseFloat(String(r.solubility)) || null : null,
      density: r.density ? parseFloat(String(r.density)) || null : null,
      molecular_weight: r.molecular_weight ? parseFloat(String(r.molecular_weight)) || null : null,
      source,
      is_verified: false,
    }));

    // Insert into Supabase
    const { data, error } = await supabase.from("materials").insert(insertRows).select();

    if (error) {
      return NextResponse.json({ success: false, error: error.message }, { status: 500 });
    }

    // Get total count
    const { count } = await supabase.from("materials").select("*", { count: "exact", head: true });

    return NextResponse.json({
      success: true,
      data: {
        filename: file.name,
        rows_total: records.length,
        inserted: data?.length || 0,
        skipped: records.length - withSmiles.length,
        total_in_db: count || 0,
      },
    });
  } catch (e) {
    return NextResponse.json({ success: false, error: String(e) }, { status: 500 });
  }
}
