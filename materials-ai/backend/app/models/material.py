"""
SQLAlchemy ORM models for the materials database.
"""

from datetime import datetime
from sqlalchemy import (
    Column,
    Integer,
    String,
    Float,
    Text,
    Boolean,
    DateTime,
    ForeignKey,
    Index,
)
from sqlalchemy.orm import relationship
from pgvector.sqlalchemy import Vector

from app.core.database import Base


class Material(Base):
    """A chemical material with molecular properties and embeddings."""

    __tablename__ = "materials"

    id = Column(Integer, primary_key=True, autoincrement=True)
    name = Column(String(200), nullable=True, index=True)
    smiles = Column(Text, nullable=False)
    category = Column(String(100), nullable=True, index=True)

    # Molecular descriptors
    molecular_weight = Column(Float, nullable=True)
    logp = Column(Float, nullable=True)
    hbd = Column(Integer, nullable=True)  # hydrogen bond donors
    hba = Column(Integer, nullable=True)  # hydrogen bond acceptors
    tpsa = Column(Float, nullable=True)  # topological polar surface area
    rotatable_bonds = Column(Integer, nullable=True)
    aromatic_rings = Column(Integer, nullable=True)

    # Target properties (predicted or measured)
    thermal_stability = Column(Float, nullable=True)  # °C
    dielectric_constant = Column(Float, nullable=True)
    bandgap = Column(Float, nullable=True)  # eV
    solubility = Column(Float, nullable=True)
    density = Column(Float, nullable=True)  # g/cm³

    # Vector embeddings
    embedding = Column(Vector(768), nullable=True)  # Uni-Mol embedding

    # Metadata
    source = Column(String(100), nullable=True)  # patent, experiment, database, etc.
    is_verified = Column(Boolean, default=False)
    created_at = Column(DateTime, default=datetime.utcnow)
    updated_at = Column(DateTime, default=datetime.utcnow, onupdate=datetime.utcnow)

    # Relationships
    experiments = relationship("Experiment", back_populates="material")

    __table_args__ = (
        Index(
            "ix_materials_embedding_cosine",
            "embedding",
            postgresql_using="ivfflat",
            postgresql_with={"lists": 100},
            postgresql_ops={"embedding": "vector_cosine_ops"},
        ),
    )

    def __repr__(self):
        return f"<Material(id={self.id}, name={self.name}, smiles={self.smiles[:30]})>"


class Experiment(Base):
    """
    Experiment result submitted by a researcher.
    Used for the feedback loop / active learning.
    """

    __tablename__ = "experiments"

    id = Column(Integer, primary_key=True, autoincrement=True)
    material_id = Column(Integer, ForeignKey("materials.id"), nullable=False)

    # Actual measured properties
    actual_thermal_stability = Column(Float, nullable=True)
    actual_dielectric_constant = Column(Float, nullable=True)
    actual_bandgap = Column(Float, nullable=True)
    actual_solubility = Column(Float, nullable=True)
    actual_density = Column(Float, nullable=True)

    success = Column(Boolean, default=False)
    notes = Column(Text, nullable=True)
    researcher = Column(String(100), nullable=True)

    created_at = Column(DateTime, default=datetime.utcnow)

    # Relationships
    material = relationship("Material", back_populates="experiments")

    def __repr__(self):
        return f"<Experiment(id={self.id}, material_id={self.material_id}, success={self.success})>"
