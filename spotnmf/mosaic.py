"""Adapter that exposes spOT-NMF as a factorization backend for mosaicMPI.

`mosaicMPI <https://github.com/MorrissyLab/mosaicMPI>`_ drives a consensus
factorization pipeline whose only algorithm-specific step is factorizing a
single (cells x genes) matrix into a spectra matrix (programs x genes) and a
usage matrix (cells x programs). mosaicMPI 3.x lets that step be swapped for an
alternative backend registered through ``mosaicmpi.factorization``.

This module provides such a backend, :func:`factorize_kernel`, that runs the
optimal-transport NMF from :func:`spotnmf.models.run_spotnmf`. Importing this
module registers the backend under the name ``"spotnmf"`` (if mosaicMPI is
installed), so the following just works::

    mosaicmpi initialize-cnmf -n mydata -k 10 --algorithm spotnmf
    mosaicmpi factorize -n mydata
    mosaicmpi postprocess -n mydata

The consensus, k-selection and on-disk output format are entirely mosaicMPI's;
spOT-NMF only supplies the per-iteration factorization. Each mosaicMPI iteration
(a distinct random seed) produces one OT-NMF solution, and mosaicMPI clusters
them into a consensus solution exactly as it does for scikit-learn NMF.

Note on inputs: mosaicMPI passes a unit-variance-scaled highly-variable-gene
matrix (non-negative). spOT-NMF column-normalizes it internally onto the
probability simplex before computing the entropic-OT objective.
"""

import numpy as np
import pandas as pd
import anndata as ad

from .models import run_spotnmf


# Default OT-NMF parameters. run_spotnmf() requires these keys; callers may
# override any of them via mosaicMPI's --algorithm_param KEY VALUE (or the
# ``factorizer_params`` dict passed to ``Dataset.initialize_cnmf``).
DEFAULT_PARAMS = {
    "h": 1e-2,             # entropic regularization on the spectra (gene loadings)
    "w": 1e-2,             # entropic regularization on the usage (spot loadings)
    "eps": 2e-2,           # entropic-OT smoothing
    "lr": 1e-2,            # optimizer learning rate (tuned for the adam default)
    "normalize_rows": True,
    "cost": "cosine",      # ground-cost metric between genes
}


def factorize_kernel(X, n_components, random_state, var_names=None, params=None):
    """Factorize a single matrix with spOT-NMF, matching mosaicMPI's contract.

    :param X: Non-negative expression matrix, cells x genes (dense or sparse).
    :param n_components: Rank k (number of programs) to learn.
    :type n_components: int
    :param random_state: Random seed for reproducibility of this iteration.
    :type random_state: int
    :param var_names: Optional gene labels for the columns of ``X``. Used to keep
        the returned spectra aligned to mosaicMPI's gene order.
    :param params: Optional backend-specific overrides for :data:`DEFAULT_PARAMS`
        and any other keyword accepted by :func:`spotnmf.models.run_spotnmf`
        (e.g. ``cost``, ``optim_name``, ``max_iter``, ``device``).
    :type params: dict, optional
    :return: ``(spectra, usages)`` where ``spectra`` is a (k x genes) ndarray and
        ``usages`` is a (cells x k) ndarray -- the same orientation as
        scikit-learn's ``non_negative_factorization`` return value.
    :rtype: tuple(numpy.ndarray, numpy.ndarray)
    """
    # Densify so AnnData/torch see a plain non-negative float array.
    if hasattr(X, "toarray"):
        X = X.toarray()
    X = np.asarray(X, dtype=np.float64)

    n_obs, n_var = X.shape
    if var_names is None:
        var_names = [str(i) for i in range(n_var)]
    else:
        var_names = list(map(str, var_names))

    # Wrap into an AnnData that spOT-NMF understands: every gene is treated as
    # highly variable (mosaicMPI has already performed HVG selection upstream).
    adata = ad.AnnData(
        X=X,
        obs=pd.DataFrame(index=[str(i) for i in range(n_obs)]),
        var=pd.DataFrame(index=var_names),
    )
    adata.var["highly_variable"] = True

    # Merge defaults with caller-supplied overrides.
    model_params = dict(DEFAULT_PARAMS)
    if params:
        model_params.update(params)

    results, _losses = run_spotnmf(
        adata,
        components=int(n_components),
        seed=int(random_state),
        **model_params,
    )

    # run_spotnmf returns:
    #   topics_per_spot: cells x programs (usage W)
    #   genes_per_topic: genes x programs (spectra H, transposed vs. mosaicMPI)
    usages = np.asarray(results["topics_per_spot"].to_numpy(), dtype=np.float64)

    genes_per_topic = results["genes_per_topic"]
    # Realign spectra rows to mosaicMPI's gene order, then transpose to k x genes.
    genes_per_topic = genes_per_topic.reindex(var_names)
    spectra = np.asarray(genes_per_topic.to_numpy().T, dtype=np.float64)

    return spectra, usages


def register(register_fn=None, name="spotnmf"):
    """Register :func:`factorize_kernel` as a mosaicMPI factorization backend.

    Called automatically on import when mosaicMPI is available. Also exposed for
    explicit use::

        from mosaicmpi.factorization import register_factorizer
        import spotnmf.mosaic
        spotnmf.mosaic.register(register_factorizer)

    :param register_fn: mosaicMPI's ``register_factorizer``. If ``None``, it is
        imported from :mod:`mosaicmpi.factorization`.
    :param name: Backend name to register under. Defaults to ``"spotnmf"``.
    """
    if register_fn is None:
        from mosaicmpi.factorization import register_factorizer as register_fn
    register_fn(name, factorize_kernel)


# Eagerly self-register when mosaicMPI is installed, so that simply importing
# this module (which mosaicMPI does when --algorithm spotnmf is requested) is
# enough to make the backend available. Stay silent if mosaicMPI is absent.
try:
    register()
except Exception:
    pass
