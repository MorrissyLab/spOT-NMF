from ._version import __version__

# Expose important submodules for easier access in notebooks
from . import io, models, pl, hvg, annotate, enrichment, niche_networks, utils, eval, gscore, cli

__all__ = [
    "__version__", "io", "models", "pl", "hvg", "annotate",
    "enrichment", "niche_networks", "utils", "eval", "gscore", "cli"
]