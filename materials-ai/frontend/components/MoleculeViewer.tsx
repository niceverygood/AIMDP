"use client";

import { useEffect, useRef, useState } from "react";

interface MoleculeViewerProps {
  smiles: string;
  width?: number;
  height?: number;
}

export function MoleculeViewer({
  smiles,
  width = 300,
  height = 300,
}: MoleculeViewerProps) {
  const containerRef = useRef<HTMLDivElement>(null);
  const [loading, setLoading] = useState(true);
  const [error, setError] = useState<string | null>(null);

  useEffect(() => {
    if (!containerRef.current) return;

    const loadViewer = async () => {
      setLoading(true);
      setError(null);

      try {
        // Load 3Dmol.js dynamically
        if (!(window as unknown as Record<string, unknown>).$3Dmol) {
          await new Promise<void>((resolve, reject) => {
            const script = document.createElement("script");
            script.src =
              "https://cdnjs.cloudflare.com/ajax/libs/3Dmol/2.0.4/3Dmol-min.js";
            script.onload = () => resolve();
            script.onerror = () => reject(new Error("Failed to load 3Dmol.js"));
            document.head.appendChild(script);
          });
        }

        const $3Dmol = (window as unknown as Record<string, unknown>).$3Dmol as {
          createViewer: (
            el: HTMLDivElement,
            opts: Record<string, unknown>
          ) => {
            addModel: (data: string, format: string) => void;
            setStyle: (sel: Record<string, unknown>, style: Record<string, unknown>) => void;
            zoomTo: () => void;
            render: () => void;
            spin: (axis: string, speed?: number) => void;
            clear: () => void;
          };
        };

        const viewer = $3Dmol.createViewer(containerRef.current!, {
          backgroundColor: "#0A1628",
        });

        // Fetch SDF from backend
        const apiUrl =
          process.env.NEXT_PUBLIC_API_URL || "http://localhost:8000";
        const res = await fetch(
          `${apiUrl}/api/materials/0/sdf?smiles=${encodeURIComponent(smiles)}`
        );

        if (res.ok) {
          const sdf = await res.text();
          viewer.addModel(sdf, "sdf");
        } else {
          // Fallback: show a simple representation using SMILES
          // Generate a basic 3D model from SMILES via 3Dmol
          viewer.addModel(smiles, "smi");
        }

        viewer.setStyle(
          {},
          { stick: { colorscheme: "Jmol", radius: 0.15 } }
        );
        viewer.zoomTo();
        viewer.render();
        viewer.spin("y", 0.5);
        setLoading(false);
      } catch (err) {
        setError(err instanceof Error ? err.message : "Viewer error");
        setLoading(false);
      }
    };

    loadViewer();
  }, [smiles]);

  return (
    <div
      className="relative rounded-xl overflow-hidden border border-white/[0.06]"
      style={{ width, height }}
    >
      {loading && (
        <div className="absolute inset-0 flex items-center justify-center bg-[#0A1628] z-10">
          <div className="flex flex-col items-center gap-2">
            <div className="w-8 h-8 border-2 border-[#00B4D8] border-t-transparent rounded-full animate-spin" />
            <span className="text-xs text-[#8892B0]">분자 로딩중...</span>
          </div>
        </div>
      )}
      {error && (
        <div className="absolute inset-0 flex items-center justify-center bg-[#0A1628] z-10">
          <div className="text-center px-4">
            <svg
              className="w-8 h-8 text-[#FF9100] mx-auto mb-2"
              fill="none"
              viewBox="0 0 24 24"
              stroke="currentColor"
              strokeWidth={1.5}
            >
              <path
                strokeLinecap="round"
                strokeLinejoin="round"
                d="M12 9v3.75m9-.75a9 9 0 1 1-18 0 9 9 0 0 1 18 0Zm-9 3.75h.008v.008H12v-.008Z"
              />
            </svg>
            <p className="text-xs text-[#8892B0]">3D 뷰어를 로드할 수 없습니다</p>
          </div>
        </div>
      )}
      <div ref={containerRef} style={{ width, height }} />
      <div className="absolute bottom-2 left-2 bg-black/50 backdrop-blur-sm rounded px-2 py-1">
        <code className="text-[10px] text-[#00B4D8] font-mono">
          {smiles.length > 30 ? smiles.slice(0, 30) + "..." : smiles}
        </code>
      </div>
    </div>
  );
}
