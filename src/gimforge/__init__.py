"""GIMForge: construction of Genetically Influenced Metabotypes."""

from .parameters import GIMParameters, parameters_from_args
from .pipeline import run_gim
from .components import components_from_matrix
from .regions import clump_sentinels

__all__ = [
    "GIMParameters",
    "parameters_from_args",
    "run_gim",
    "clump_sentinels",
    "components_from_matrix",
]
__version__ = "0.3.0"
