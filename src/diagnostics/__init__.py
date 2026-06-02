"""Adequacy diagnostics: finesse (lineshape & cavity) and the round-trip truncation factor."""
from .finesse import (
    lineshape_finesse, finesse_to_reflectance, reflectance_to_finesse,
    cavity_finesse_from_R, cavity_finesse_from_roundtrip_amplitude,
    band_q_statistics, two_beam_relative_error,
)

__all__ = [
    "lineshape_finesse", "finesse_to_reflectance", "reflectance_to_finesse",
    "cavity_finesse_from_R", "cavity_finesse_from_roundtrip_amplitude",
    "band_q_statistics", "two_beam_relative_error",
]
