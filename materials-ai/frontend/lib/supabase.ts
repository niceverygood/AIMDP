import { createClient } from "@supabase/supabase-js";

const SUPABASE_URL = process.env.NEXT_PUBLIC_SUPABASE_URL || "https://hkzefauqitxychpnjkra.supabase.co";
const SUPABASE_KEY = process.env.NEXT_PUBLIC_SUPABASE_KEY || "sb_publishable_ryae9bsqfx_hko78SbbcTg_LrNEaL3U";

export const supabase = createClient(SUPABASE_URL, SUPABASE_KEY);

export type MaterialRow = {
  id?: number;
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
  created_at?: string;
};
