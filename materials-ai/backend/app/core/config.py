"""
Application configuration — loads from environment variables.
"""

from pydantic_settings import BaseSettings
from functools import lru_cache
from typing import Optional


class Settings(BaseSettings):
    """Global configuration loaded from .env or environment variables."""

    # Application
    APP_NAME: str = "Materials AI Platform"
    APP_VERSION: str = "0.1.0"
    DEBUG: bool = True

    # Database
    DATABASE_URL: str = "postgresql+asyncpg://user:pass@localhost:5432/materials"
    DATABASE_ECHO: bool = False

    # Redis
    REDIS_URL: str = "redis://localhost:6379"

    # CORS
    CORS_ORIGINS: list[str] = ["http://localhost:3000"]

    # ML Models
    UNIMOL_MODEL_PATH: str = "./ml/unimol"
    PROPERTY_MODEL_PATH: str = "./ml/property_models"

    # Embeddings
    EMBEDDING_DIM: int = 768

    # Search
    DEFAULT_SEARCH_LIMIT: int = 20
    MAX_SEARCH_LIMIT: int = 100

    # Pipeline
    UPLOAD_DIR: str = "./data/uploads"
    PROCESSED_DIR: str = "./data/processed"

    # Security
    API_KEY: Optional[str] = None

    model_config = {
        "env_file": ".env",
        "env_file_encoding": "utf-8",
        "case_sensitive": True,
    }


@lru_cache()
def get_settings() -> Settings:
    return Settings()
