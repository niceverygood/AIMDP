"""
Patent data extractor — parses patent documents (PDF, XML)
to extract molecular formulas, SMILES, and property data.
"""

import logging
import re
from typing import Optional

logger = logging.getLogger(__name__)


class PatentParser:
    """
    Extract molecular information from patent documents.
    Uses regex patterns and NER for chemical entity recognition.
    """

    # Common patterns for chemical formulas
    FORMULA_PATTERN = re.compile(
        r"\b[A-Z][a-z]?\d*(?:[A-Z][a-z]?\d*)*\b"
    )
    # SMILES-like patterns
    SMILES_PATTERN = re.compile(
        r"[A-Za-z0-9@+\-\[\]\(\)\\\/=#%\.:]+(?:\.[A-Za-z0-9@+\-\[\]\(\)\\\/=#%\.:]+)*"
    )

    def parse_text(self, text: str) -> list[dict]:
        """
        Extract chemical data from patent text.

        Returns list of dicts with:
          - formula: molecular formula
          - smiles: SMILES string (if found)
          - properties: extracted property values
          - context: surrounding text
        """
        results = []

        # Split into paragraphs
        paragraphs = text.split("\n\n")

        for para in paragraphs:
            formulas = self._extract_formulas(para)
            properties = self._extract_properties(para)

            if formulas or properties:
                results.append({
                    "formulas": formulas,
                    "properties": properties,
                    "context": para[:500],
                })

        return results

    def _extract_formulas(self, text: str) -> list[str]:
        """Extract chemical formulas from text."""
        candidates = self.FORMULA_PATTERN.findall(text)
        # Filter to likely molecular formulas (must contain C)
        return [f for f in candidates if "C" in f and len(f) > 2]

    def _extract_properties(self, text: str) -> dict:
        """Extract property values from text using pattern matching."""
        properties = {}

        # Temperature patterns (e.g., "thermal stability of 350°C")
        temp_match = re.search(
            r"(?:thermal|열안정|분해온도|Td)\s*[:\=]?\s*(\d+(?:\.\d+)?)\s*°?C",
            text,
            re.IGNORECASE,
        )
        if temp_match:
            properties["thermal_stability"] = float(temp_match.group(1))

        # Dielectric constant patterns
        dk_match = re.search(
            r"(?:dielectric|유전율|Dk|ε)\s*[:\=]?\s*(\d+(?:\.\d+)?)",
            text,
            re.IGNORECASE,
        )
        if dk_match:
            properties["dielectric_constant"] = float(dk_match.group(1))

        # Bandgap patterns
        bg_match = re.search(
            r"(?:bandgap|band\s*gap|밴드갭|Eg)\s*[:\=]?\s*(\d+(?:\.\d+)?)\s*eV",
            text,
            re.IGNORECASE,
        )
        if bg_match:
            properties["bandgap"] = float(bg_match.group(1))

        return properties

    async def parse_pdf(self, file_path: str) -> list[dict]:
        """Parse a patent PDF file."""
        try:
            import fitz  # PyMuPDF

            doc = fitz.open(file_path)
            full_text = ""
            for page in doc:
                full_text += page.get_text() + "\n\n"
            doc.close()

            return self.parse_text(full_text)
        except ImportError:
            logger.error("PyMuPDF not installed. pip install pymupdf")
            return []

    async def parse_xml(self, file_path: str) -> list[dict]:
        """Parse a patent XML file."""
        try:
            import xml.etree.ElementTree as ET

            tree = ET.parse(file_path)
            root = tree.getroot()

            # Extract text from all text elements
            texts = [elem.text for elem in root.iter() if elem.text]
            full_text = "\n\n".join(texts)

            return self.parse_text(full_text)
        except Exception as e:
            logger.error(f"XML parsing failed: {e}")
            return []
