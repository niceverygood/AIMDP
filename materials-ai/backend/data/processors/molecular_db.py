"""
Molecular database connector — fetches data from external public databases.
Supports PubChem, Materials Project, QM9, etc.
"""

import logging
from typing import Optional

logger = logging.getLogger(__name__)


class MolecularDBConnector:
    """
    Fetch molecular data from external databases for dataset enrichment.
    """

    async def fetch_pubchem(
        self,
        query: str,
        max_results: int = 100,
    ) -> list[dict]:
        """
        Search PubChem for molecules matching a query.
        Uses PubChem PUG REST API.
        """
        import httpx

        base_url = "https://pubchem.ncbi.nlm.nih.gov/rest/pug"

        try:
            async with httpx.AsyncClient(timeout=30) as client:
                # Search by name or SMILES
                search_url = f"{base_url}/compound/name/{query}/property/MolecularFormula,MolecularWeight,CanonicalSMILES,IsomericSMILES,XLogP,TPSA,HBondDonorCount,HBondAcceptorCount/JSON"
                resp = await client.get(search_url)

                if resp.status_code == 200:
                    data = resp.json()
                    properties = data.get("PropertyTable", {}).get("Properties", [])

                    results = []
                    for prop in properties[:max_results]:
                        results.append({
                            "smiles": prop.get("CanonicalSMILES", ""),
                            "name": query,
                            "molecular_weight": prop.get("MolecularWeight"),
                            "logp": prop.get("XLogP"),
                            "tpsa": prop.get("TPSA"),
                            "hbd": prop.get("HBondDonorCount"),
                            "hba": prop.get("HBondAcceptorCount"),
                            "source": "pubchem",
                        })

                    return results
                else:
                    logger.warning(f"PubChem search failed: {resp.status_code}")
                    return []
        except Exception as e:
            logger.error(f"PubChem fetch error: {e}")
            return []

    async def fetch_materials_project(
        self,
        api_key: str,
        elements: Optional[list[str]] = None,
        max_results: int = 100,
    ) -> list[dict]:
        """
        Fetch materials from Materials Project API.
        Requires an API key from materialsproject.org.
        """
        try:
            from mp_api.client import MPRester

            with MPRester(api_key) as mpr:
                if elements:
                    docs = mpr.materials.summary.search(
                        elements=elements,
                        num_chunks=1,
                        chunk_size=max_results,
                    )
                else:
                    docs = mpr.materials.summary.search(
                        num_chunks=1,
                        chunk_size=max_results,
                    )

                results = []
                for doc in docs:
                    results.append({
                        "name": str(doc.material_id),
                        "formula": str(doc.composition),
                        "bandgap": doc.band_gap if hasattr(doc, "band_gap") else None,
                        "density": doc.density if hasattr(doc, "density") else None,
                        "source": "materials_project",
                    })

                return results
        except ImportError:
            logger.error("mp-api not installed. pip install mp-api")
            return []
        except Exception as e:
            logger.error(f"Materials Project fetch error: {e}")
            return []

    async def load_qm9_dataset(self, data_dir: str = "./data/qm9") -> list[dict]:
        """
        Load QM9 dataset (134k molecules with 19 quantum properties).
        Great for PoC property prediction models.
        """
        try:
            from torch_geometric.datasets import QM9

            dataset = QM9(root=data_dir)
            results = []

            for data in dataset:
                # QM9 properties:
                # 0: dipole_moment, 1: isotropic_polarizability,
                # 2: HOMO, 3: LUMO, 4: gap, ...
                results.append({
                    "smiles": data.smiles if hasattr(data, "smiles") else "",
                    "dipole_moment": float(data.y[0][0]) if data.y is not None else None,
                    "homo": float(data.y[0][2]) if data.y is not None else None,
                    "lumo": float(data.y[0][3]) if data.y is not None else None,
                    "bandgap": float(data.y[0][4]) if data.y is not None else None,
                    "source": "qm9",
                })

            return results
        except ImportError:
            logger.error("torch-geometric not installed")
            return []
