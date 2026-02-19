/**
 * TypeScript type definitions for the Materials AI Platform.
 */

// ---------- API Response ----------

export interface APIResponse<T> {
  success: boolean;
  data: T | null;
  error?: string;
}

// ---------- Material ----------

export interface Material {
  id: number;
  name: string | null;
  smiles: string;
  category: string | null;
  molecular_weight: number | null;
  logp: number | null;
  hbd: number | null;
  hba: number | null;
  tpsa: number | null;
  rotatable_bonds: number | null;
  aromatic_rings: number | null;
  thermal_stability: number | null;
  dielectric_constant: number | null;
  bandgap: number | null;
  solubility: number | null;
  density: number | null;
  source: string | null;
  is_verified: boolean;
  created_at: string;
  updated_at: string;
}

export interface MaterialSearchResult extends Material {
  similarity: number;
  match_score: number;
}

// ---------- Search ----------

export interface SearchQuery {
  smiles?: string;
  category?: string;
  min_thermal_stability?: number;
  max_thermal_stability?: number;
  min_dielectric?: number;
  max_dielectric?: number;
  min_bandgap?: number;
  max_bandgap?: number;
  min_solubility?: number;
  max_solubility?: number;
  weight_thermal?: number;
  weight_dielectric?: number;
  weight_bandgap?: number;
  weight_similarity?: number;
  limit?: number;
  offset?: number;
}

export interface SearchResponse {
  results: MaterialSearchResult[];
  total: number;
  query_time_ms: number;
}

// ---------- Pipeline ----------

export type PipelineStepStatusType = "pending" | "running" | "completed" | "failed";

export interface PipelineStepStatus {
  name: string;
  label: string;
  status: PipelineStepStatusType;
  progress: number;
  rows_processed: number;
  total_rows: number;
  quality_score: number | null;
  time_elapsed_sec: number;
  error: string | null;
}

export interface PipelineStatus {
  pipeline_id: string;
  status: string;
  steps: PipelineStepStatus[];
  started_at: string | null;
  completed_at: string | null;
}

// ---------- Property Prediction ----------

export interface PropertyPrediction {
  smiles: string;
  thermal_stability: number | null;
  dielectric_constant: number | null;
  bandgap: number | null;
  solubility: number | null;
  density: number | null;
  confidence: number | null;
}

// ---------- Experiment Feedback ----------

export interface ExperimentResult {
  material_id: number;
  actual_thermal_stability?: number;
  actual_dielectric_constant?: number;
  actual_bandgap?: number;
  actual_solubility?: number;
  actual_density?: number;
  success: boolean;
  notes?: string;
  researcher?: string;
}

export interface FeedbackAnalysis {
  material_id: number;
  prediction_errors: Record<string, number>;
  model_improvement: string;
  retrain_triggered: boolean;
}

// ---------- Dashboard Stats ----------

export interface DashboardStats {
  total_materials: number;
  verified_materials: number;
  categories: Record<string, number>;
}

// ---------- Search Filter Config ----------

export interface PropertyFilter {
  key: string;
  label: string;
  unit: string;
  min: number;
  max: number;
  step: number;
}

export const PROPERTY_FILTERS: PropertyFilter[] = [
  { key: "thermal_stability", label: "열안정성", unit: "°C", min: 0, max: 500, step: 10 },
  { key: "dielectric", label: "유전율", unit: "", min: 0, max: 25, step: 0.5 },
  { key: "bandgap", label: "밴드갭", unit: "eV", min: 0, max: 8, step: 0.1 },
  { key: "solubility", label: "용해도", unit: "", min: 0, max: 1, step: 0.01 },
];

// ---------- Material Categories ----------

export const MATERIAL_CATEGORIES = [
  { value: "OLED", label: "OLED" },
  { value: "organic", label: "유기 분자" },
  { value: "battery", label: "배터리" },
  { value: "semiconductor", label: "반도체" },
  { value: "hard_coating", label: "하드코팅" },
  { value: "display", label: "디스플레이" },
] as const;
