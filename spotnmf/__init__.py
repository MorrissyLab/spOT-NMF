from ._version import __version__

# Import lightweight internal modules
from . import io, models, hvg, annotate, enrichment, niche_networks, utils, eval, gscore, cli

# Lazy-loading for heavy modules
def __getattr__(name):
    if name == "pl":
        from . import pl
        return pl
    if name == "cnmf":
        from .external import cnmf
        return cnmf
    raise AttributeError(name)

__all__ = [
    "__version__",
    "io", "models", "hvg", "annotate", "enrichment",
    "niche_networks", "utils", "eval", "gscore", "cli",
    "pl", "cnmf"
]
