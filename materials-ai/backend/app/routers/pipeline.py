"""
Data Pipeline API — monitor and control the data processing pipeline.
Supports WebSocket for real-time progress tracking.
"""

import uuid
from fastapi import APIRouter, Depends, HTTPException, WebSocket, WebSocketDisconnect
from sqlalchemy.ext.asyncio import AsyncSession

from app.core.database import get_db
from app.models.schemas import (
    APIResponse,
    PipelineStatus,
    PipelineStepStatus,
    PipelineStartRequest,
)
from app.services.pipeline import PipelineOrchestrator

router = APIRouter()

# In-memory pipeline state (use Redis for production)
pipeline_states: dict[str, PipelineStatus] = {}
active_connections: dict[str, list[WebSocket]] = {}

orchestrator = PipelineOrchestrator()


def _default_steps() -> list[PipelineStepStatus]:
    return [
        PipelineStepStatus(name="collect", label="데이터 수집", status="pending"),
        PipelineStepStatus(name="normalize", label="정규화", status="pending"),
        PipelineStepStatus(name="clean", label="데이터 정제", status="pending"),
        PipelineStepStatus(name="engineer", label="특성 추출", status="pending"),
        PipelineStepStatus(name="embed", label="임베딩 생성", status="pending"),
        PipelineStepStatus(name="dataset", label="데이터셋 생성", status="pending"),
    ]


@router.post("/start", response_model=APIResponse[PipelineStatus])
async def start_pipeline(
    request: PipelineStartRequest,
    db: AsyncSession = Depends(get_db),
):
    """
    새 데이터 파이프라인 실행 시작.
    WebSocket으로 실시간 진행률 모니터링 가능.
    """
    pipeline_id = str(uuid.uuid4())[:8]

    status = PipelineStatus(
        pipeline_id=pipeline_id,
        status="running",
        steps=_default_steps(),
    )
    pipeline_states[pipeline_id] = status

    # Start pipeline in background
    import asyncio

    asyncio.create_task(
        orchestrator.run_pipeline(
            pipeline_id=pipeline_id,
            source_type=request.source_type,
            file_path=request.file_path,
            options=request.options,
            progress_callback=_on_progress,
        )
    )

    return APIResponse(success=True, data=status)


@router.get("/{pipeline_id}", response_model=APIResponse[PipelineStatus])
async def get_pipeline_status(pipeline_id: str):
    """파이프라인 진행 상태 조회."""
    status = pipeline_states.get(pipeline_id)
    if not status:
        raise HTTPException(status_code=404, detail="Pipeline not found")
    return APIResponse(success=True, data=status)


@router.get("", response_model=APIResponse[list[PipelineStatus]])
async def list_pipelines():
    """전체 파이프라인 목록 조회."""
    return APIResponse(success=True, data=list(pipeline_states.values()))


@router.websocket("/ws/{pipeline_id}")
async def pipeline_websocket(websocket: WebSocket, pipeline_id: str):
    """실시간 파이프라인 진행 상태 WebSocket 엔드포인트."""
    await websocket.accept()

    if pipeline_id not in active_connections:
        active_connections[pipeline_id] = []
    active_connections[pipeline_id].append(websocket)

    try:
        # Send current state immediately
        if pipeline_id in pipeline_states:
            await websocket.send_json(pipeline_states[pipeline_id].model_dump(mode="json"))

        # Keep connection alive, listening for close
        while True:
            try:
                await websocket.receive_text()
            except WebSocketDisconnect:
                break
    finally:
        if pipeline_id in active_connections:
            active_connections[pipeline_id].remove(websocket)


async def _on_progress(pipeline_id: str, step_name: str, step_status: dict):
    """Callback invoked by the orchestrator to report progress."""
    if pipeline_id not in pipeline_states:
        return

    status = pipeline_states[pipeline_id]
    for step in status.steps:
        if step.name == step_name:
            step.status = step_status.get("status", step.status)
            step.progress = step_status.get("progress", step.progress)
            step.rows_processed = step_status.get("rows_processed", step.rows_processed)
            step.total_rows = step_status.get("total_rows", step.total_rows)
            step.quality_score = step_status.get("quality_score", step.quality_score)
            step.time_elapsed_sec = step_status.get("time_elapsed_sec", step.time_elapsed_sec)
            step.error = step_status.get("error")
            break

    # Check if all steps completed
    if all(s.status in ("completed", "failed") for s in status.steps):
        status.status = "completed" if all(s.status == "completed" for s in status.steps) else "failed"

    # Broadcast to WebSocket clients
    if pipeline_id in active_connections:
        msg = status.model_dump(mode="json")
        for ws in active_connections[pipeline_id]:
            try:
                await ws.send_json(msg)
            except Exception:
                pass
