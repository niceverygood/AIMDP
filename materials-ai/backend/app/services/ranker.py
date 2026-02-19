"""
Multi-objective ranking service.
Combines vector similarity with property-matching to produce a single match score.
"""

from typing import Optional


class MultiObjectiveRanker:
    """
    Computes a composite score for material search results.

    Score = Σ (weight_i × normalized_score_i) for each objective.

    Objectives:
      - similarity: cosine similarity from vector search
      - thermal_stability: how close to target
      - dielectric_constant: how close to target
      - bandgap: how close to target
    """

    # Normalization ranges for each property
    PROPERTY_RANGES = {
        "thermal_stability": (0.0, 500.0),
        "dielectric_constant": (0.0, 25.0),
        "bandgap": (0.0, 8.0),
        "solubility": (0.0, 1.0),
        "density": (0.5, 3.0),
    }

    def compute_score(
        self,
        similarity: float = 0.0,
        thermal_stability: Optional[float] = None,
        dielectric_constant: Optional[float] = None,
        bandgap: Optional[float] = None,
        solubility: Optional[float] = None,
        density: Optional[float] = None,
        weights: Optional[dict[str, float]] = None,
        targets: Optional[dict[str, float]] = None,
    ) -> float:
        """
        Compute multi-objective match score in [0, 1].

        Args:
            similarity: Vector cosine similarity (0-1)
            thermal_stability: Material's thermal stability value
            dielectric_constant: Material's dielectric constant
            bandgap: Material's bandgap
            solubility: Material's solubility
            density: Material's density
            weights: Custom weights per objective
            targets: Target values for each property (if proximity scoring needed)

        Returns:
            Combined score in [0, 1]
        """
        if weights is None:
            weights = {
                "similarity": 0.4,
                "thermal": 0.2,
                "dielectric": 0.2,
                "bandgap": 0.2,
            }

        total_score = 0.0
        total_weight = 0.0

        # Similarity score (already 0-1)
        sim_weight = weights.get("similarity", 0.0)
        if sim_weight > 0:
            total_score += sim_weight * max(0, min(1, similarity))
            total_weight += sim_weight

        # Property scores — normalize to [0, 1]
        property_values = {
            "thermal": thermal_stability,
            "dielectric": dielectric_constant,
            "bandgap": bandgap,
            "solubility": solubility,
            "density": density,
        }
        property_range_keys = {
            "thermal": "thermal_stability",
            "dielectric": "dielectric_constant",
            "bandgap": "bandgap",
            "solubility": "solubility",
            "density": "density",
        }

        for key, value in property_values.items():
            w = weights.get(key, 0.0)
            if w > 0 and value is not None:
                range_key = property_range_keys[key]
                low, high = self.PROPERTY_RANGES[range_key]
                normalized = (value - low) / (high - low) if high != low else 0.0
                normalized = max(0.0, min(1.0, normalized))

                if targets and key in targets:
                    # Proximity scoring: penalize distance from target
                    target_norm = (targets[key] - low) / (high - low)
                    proximity = 1.0 - abs(normalized - target_norm)
                    total_score += w * max(0, proximity)
                else:
                    total_score += w * normalized

                total_weight += w

        # Normalize by total weight
        if total_weight > 0:
            return round(total_score / total_weight, 4)

        return 0.0

    def rank_results(
        self,
        results: list[dict],
        weights: Optional[dict[str, float]] = None,
    ) -> list[dict]:
        """
        Rank a list of result dicts by multi-objective score.
        """
        for r in results:
            r["match_score"] = self.compute_score(
                similarity=r.get("similarity", 0.0),
                thermal_stability=r.get("thermal_stability"),
                dielectric_constant=r.get("dielectric_constant"),
                bandgap=r.get("bandgap"),
                weights=weights,
            )
        return sorted(results, key=lambda r: r["match_score"], reverse=True)
