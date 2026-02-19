"""
FastAPI application entry point.
AI-powered Materials Discovery Platform — Backend API
"""

from contextlib import asynccontextmanager
from fastapi import FastAPI
from fastapi.middleware.cors import CORSMiddleware

from app.core.config import get_settings
from app.core.database import init_db, close_db
from app.routers import search, pipeline, materials

settings = get_settings()


@asynccontextmanager
async def lifespan(app: FastAPI):
    """Startup / shutdown lifecycle."""
    # Startup
    await init_db()
    yield
    # Shutdown
    await close_db()


app = FastAPI(
    title=settings.APP_NAME,
    version=settings.APP_VERSION,
    description="AI-powered materials discovery platform for R&D",
    lifespan=lifespan,
)

# CORS
app.add_middleware(
    CORSMiddleware,
    allow_origins=settings.CORS_ORIGINS,
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)


# ---------- Routers ----------
app.include_router(search.router, prefix="/api/search", tags=["Search"])
app.include_router(pipeline.router, prefix="/api/pipeline", tags=["Pipeline"])
app.include_router(materials.router, prefix="/api/materials", tags=["Materials"])


# ---------- Health Check ----------
@app.get("/health")
async def health():
    return {"status": "healthy", "version": settings.APP_VERSION}


@app.get("/")
async def root():
    return {
        "message": "Materials AI Platform API",
        "docs": "/docs",
        "version": settings.APP_VERSION,
    }
