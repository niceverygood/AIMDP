/**
 * API client for the Materials AI Platform backend.
 */

import type {
  APIResponse,
  SearchQuery,
  SearchResponse,
  Material,
  PipelineStatus,
  PropertyPrediction,
  ExperimentResult,
  FeedbackAnalysis,
  DashboardStats,
} from "./types";

const API_BASE = process.env.NEXT_PUBLIC_API_URL || "http://localhost:8000";

// ---------- Generic Fetch Helper ----------

async function apiFetch<T>(
  endpoint: string,
  options?: RequestInit
): Promise<APIResponse<T>> {
  try {
    const res = await fetch(`${API_BASE}${endpoint}`, {
      headers: {
        "Content-Type": "application/json",
        ...options?.headers,
      },
      ...options,
    });

    if (!res.ok) {
      const errorData = await res.json().catch(() => ({}));
      return {
        success: false,
        data: null,
        error: errorData.detail || `HTTP ${res.status}`,
      };
    }

    return await res.json();
  } catch (error) {
    return {
      success: false,
      data: null,
      error: error instanceof Error ? error.message : "Network error",
    };
  }
}

// ---------- Search API ----------

export async function searchMaterials(
  query: SearchQuery
): Promise<APIResponse<SearchResponse>> {
  return apiFetch<SearchResponse>("/api/search", {
    method: "POST",
    body: JSON.stringify(query),
  });
}

export async function predictProperties(
  smiles: string
): Promise<APIResponse<PropertyPrediction>> {
  return apiFetch<PropertyPrediction>(
    `/api/search/predict?smiles=${encodeURIComponent(smiles)}`,
    { method: "POST" }
  );
}

// ---------- Materials CRUD ----------

export async function getMaterials(params?: {
  category?: string;
  limit?: number;
  offset?: number;
}): Promise<APIResponse<Material[]>> {
  const searchParams = new URLSearchParams();
  if (params?.category) searchParams.set("category", params.category);
  if (params?.limit) searchParams.set("limit", params.limit.toString());
  if (params?.offset) searchParams.set("offset", params.offset.toString());

  const query = searchParams.toString();
  return apiFetch<Material[]>(`/api/materials${query ? `?${query}` : ""}`);
}

export async function getMaterial(
  id: number
): Promise<APIResponse<Material>> {
  return apiFetch<Material>(`/api/materials/${id}`);
}

export async function getMaterialStats(): Promise<APIResponse<DashboardStats>> {
  return apiFetch<DashboardStats>("/api/materials/stats");
}

export async function createMaterial(
  data: Partial<Material>
): Promise<APIResponse<Material>> {
  return apiFetch<Material>("/api/materials", {
    method: "POST",
    body: JSON.stringify(data),
  });
}

// ---------- Pipeline API ----------

export async function startPipeline(params: {
  source_type: string;
  file_path: string;
  options?: Record<string, unknown>;
}): Promise<APIResponse<PipelineStatus>> {
  return apiFetch<PipelineStatus>("/api/pipeline/start", {
    method: "POST",
    body: JSON.stringify(params),
  });
}

export async function getPipelineStatus(
  pipelineId: string
): Promise<APIResponse<PipelineStatus>> {
  return apiFetch<PipelineStatus>(`/api/pipeline/${pipelineId}`);
}

export async function listPipelines(): Promise<APIResponse<PipelineStatus[]>> {
  return apiFetch<PipelineStatus[]>("/api/pipeline");
}

// ---------- Feedback API ----------

export async function submitFeedback(
  result: ExperimentResult
): Promise<APIResponse<FeedbackAnalysis>> {
  return apiFetch<FeedbackAnalysis>("/api/materials/feedback", {
    method: "POST",
    body: JSON.stringify(result),
  });
}

// ---------- Pipeline WebSocket ----------

export function connectPipelineWS(
  pipelineId: string,
  onMessage: (status: PipelineStatus) => void,
  onError?: (error: Event) => void
): WebSocket {
  const wsUrl = API_BASE.replace("http", "ws");
  const ws = new WebSocket(`${wsUrl}/api/pipeline/ws/${pipelineId}`);

  ws.onmessage = (event) => {
    try {
      const data = JSON.parse(event.data);
      onMessage(data);
    } catch {
      console.error("Failed to parse WebSocket message");
    }
  };

  ws.onerror = (event) => {
    console.error("WebSocket error:", event);
    onError?.(event);
  };

  return ws;
}
