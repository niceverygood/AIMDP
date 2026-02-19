"""
Data normalizer — standardizes units, formats, and naming conventions.
"""

import logging
import re
from typing import Optional

logger = logging.getLogger(__name__)


class DataNormalizer:
    """
    Normalize raw material data:
    - Standardize SMILES notation
    - Unify temperature units (°F → °C, K → °C)
    - Normalize property names
    - Handle missing values
    """

    # Property name aliases → canonical names
    PROPERTY_ALIASES = {
        # Thermal stability
        "td": "thermal_stability",
        "decomposition_temp": "thermal_stability",
        "thermal_decomposition": "thermal_stability",
        "열안정성": "thermal_stability",
        "분해온도": "thermal_stability",
        # Dielectric constant
        "dk": "dielectric_constant",
        "permittivity": "dielectric_constant",
        "유전율": "dielectric_constant",
        # Bandgap
        "band_gap": "bandgap",
        "eg": "bandgap",
        "energy_gap": "bandgap",
        "밴드갭": "bandgap",
        # Solubility
        "water_solubility": "solubility",
        "용해도": "solubility",
        # Density
        "밀도": "density",
        "specific_gravity": "density",
    }

    def normalize(self, record: dict) -> dict:
        """Normalize a single data record."""
        normalized = {}

        for key, value in record.items():
            # Normalize property names
            canonical_key = self._normalize_key(key)
            normalized[canonical_key] = value

        # Standardize SMILES
        if "smiles" in normalized and normalized["smiles"]:
            normalized["smiles"] = self._normalize_smiles(normalized["smiles"])

        # Convert temperature units
        if "thermal_stability" in normalized:
            normalized["thermal_stability"] = self._normalize_temperature(
                normalized["thermal_stability"],
                record.get("thermal_unit", "C"),
            )

        return normalized

    def _normalize_key(self, key: str) -> str:
        """Convert property name to canonical form."""
        k = key.strip().lower().replace(" ", "_").replace("-", "_")
        return self.PROPERTY_ALIASES.get(k, k)

    def _normalize_smiles(self, smiles: str) -> str:
        """Standardize SMILES notation using RDKit canonical SMILES."""
        try:
            from rdkit import Chem

            mol = Chem.MolFromSmiles(smiles.strip())
            if mol:
                return Chem.MolToSmiles(mol)
        except Exception:
            pass
        return smiles.strip()

    def _normalize_temperature(
        self, value: Optional[float], unit: str = "C"
    ) -> Optional[float]:
        """Convert temperature to Celsius."""
        if value is None:
            return None

        unit = unit.upper().strip()
        if unit in ("F", "°F", "FAHRENHEIT"):
            return round((value - 32) * 5 / 9, 1)
        elif unit in ("K", "KELVIN"):
            return round(value - 273.15, 1)
        return value  # Already Celsius
