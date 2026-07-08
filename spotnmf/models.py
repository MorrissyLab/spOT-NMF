"""
This file implements a modified version of the Mowgli model, transforming it into a 
1 mode of Optimal Transport Non-negative Matrix Factorization (OT-NMF) approach.
The original Mowgli model was developed for paired single-cell multi-omics data integration, 
as described in the publication: Huizing, G.-J., Deutschmann, I. M., Peyré, G., & Cantini, L. (2023). 
Paired single-cell multi-omics data integration with Mowgli. Nature Communications, 14(1), 7711.
"""

from typing import Callable, List
import numpy as np
import torch
import torch.nn.functional as F
from torch import optim
from tqdm import tqdm
import pandas as pd 
from spotnmf import utils

def seed_all(seed):
    """Seed all relevant random number generators for reproducibility.

    Seeds NumPy (which also covers SciPy) and PyTorch on both CPU and,
    when available, GPU. On GPU, cuDNN is also configured for deterministic
    behavior.

    Args:
        seed (int): The random seed to apply across NumPy and PyTorch.
    """
    np.random.seed(seed)  # Seed numpy (covers scipy)
    torch.manual_seed(seed)  # Seed pytorch (CPU)
    if torch.cuda.is_available():
        torch.cuda.manual_seed(seed)  # Seed pytorch (GPU)
        torch.cuda.manual_seed_all(seed)  # Seed all GPUs
        torch.backends.cudnn.deterministic = True  # Ensure deterministic behavior
        torch.backends.cudnn.benchmark = False

class spotnmf:
    """The spotnmf model, which performs OT-NMF.

    Args:
        factors (int, optional):
            The factors of the model. Defaults to 15.
        highly_variable (bool, optional):
            Whether to use highly variable features. Defaults to True.
            For now, only True is supported.
        use_mod_weight (bool, optional):
            Whether to use a different weight for each modality and each
            cell. If `True`, the weights are expected in the `mod_weight`
            obs field of each modality. Defaults to False.
        h_regularization (float, optional):
            The entropy parameter for the spectra. We advise setting values
            between 0.001 (biological signal driven by very few features) and 1.0
            (very diffuse biological signals).
        w_regularization (float, optional):
            The entropy parameter for the usage. As with `h_regularization`,
            small values mean sparse vectors. Defaults to 1e-2.
        eps (float, optional):
            The entropy parameter for epsilon transport. Large values
            decrease importance of individual genes. Defaults to 2e-2.
        cost (str, optional):
            The function used to compute an emprical ground cost. All
            metrics from Scipy's `cdist` are allowed. Defaults to 'cosine'.
        pca_cost (bool, optional):
            If True, the emprical ground cost will be computed on PCA
            embeddings rather than raw data. Defaults to False.
        cost_path (dict, optional):
            Will look for an existing cost as a `.npy` file at this
            path. If not found, the cost will be computed then saved
            there. Defaults to None.
    """

    def __init__(
        self,
        factors: int = 50,
        highly_variable: bool = True,
        h_regularization: float = 1e-2,
        w_regularization: float = 1e-2,
        eps: float = 2e-2,
        cost: str = "cosine",
        pca_cost: bool = False,
        cost_path: dict = None,
    ):

        # Check that the user-defined parameters are valid.
        assert factors > 0
        assert w_regularization > 0
        assert h_regularization > 0
        assert eps > 0
        assert highly_variable is True

        # Save arguments as attributes.
        self.factors = factors
        self.h_regularization = h_regularization
        self.w_regularization = w_regularization
        self.eps = eps
        self.cost = cost
        self.cost_path = cost_path
        self.pca_cost = pca_cost

        # Initialize the loss and statistics histories.
        self.losses_w, self.losses_h, self.losses = [], [], []

        self.A, self.H, self.G, self.K = None, None, None, None

    def init_parameters(
        self,
        adata_spatial,
        dtype: torch.dtype,
        device: torch.device,
        force_recompute: bool = False,
        normalize_rows: bool = False,
        impute: bool = False,
    ) -> None:
        """Initialize the model parameters based on the input data.

        Builds the reference dataset ``A`` and ground cost ``K`` from the
        spatial data, then initializes the spectra ``H``, the usage ``W``,
        and the dual variable ``G``.

        Args:
            adata_spatial (anndata.AnnData): The input spatial AnnData object.
            dtype (torch.dtype): The dtype to work with.
            device (torch.device): The device to work on.
            force_recompute (bool, optional): Whether to recompute the ground
                cost even if a cached cost is available. Defaults to False.
            normalize_rows (bool, optional): Whether to normalize the rows of
                the reference dataset before column normalization.
                Defaults to False.
            impute (bool, optional): Whether to impute the data using FastICA
                (with ``factors`` components) before building the reference
                dataset. Defaults to False.
        """

        # Set some attributes.
        self.n_obs = adata_spatial.n_obs

        # Select the highly variable features.
        if "highly_variable" not in adata_spatial.var.columns:
            keep_idx = np.ones(len(adata_spatial.var), dtype=bool)
        else:
            keep_idx = adata_spatial.var["highly_variable"].to_numpy()
        
        if impute:
            from sklearn.decomposition import FastICA
            X_data = adata_spatial.X
            X_data = X_data[:, keep_idx]

            ica = FastICA(n_components=self.factors, random_state=0)
            S = ica.fit_transform(X_data)  # Source matrix (independent components)
            X_data = ica.inverse_transform(S)
            X_data = np.clip(X_data, 0, None)  # Clip to remove negative values
        else:
            X_data = adata_spatial.X
            X_data = X_data[:, keep_idx]

        self.A = utils.reference_dataset(X_data, dtype, device)
        self.n_var = self.A.shape[0]

        # Normalize the reference dataset, and add a small value
        # for numerical stability.
        self.A += 1e-6
        if normalize_rows:
            mean_row_sum = self.A.sum(1).mean()
            self.A /= self.A.sum(1).reshape(-1, 1) * mean_row_sum
        self.A /= self.A.sum(0)

        # Determine which cost function to use.
        cost = self.cost if isinstance(self.cost, str) else self.cost
        try:
            cost_path = self.cost_path
        except Exception:
            cost_path = None

        # Define the features that the ground cost will be computed on.
        features = 1e-6 + self.A.cpu().numpy()
        if self.pca_cost:
            from sklearn.decomposition import PCA
            pca = PCA(n_components=self.factors)
            features = pca.fit_transform(features)

        # Compute ground cost, using the specified cost function.
        self.K = utils.compute_ground_cost(
            features, cost, self.eps, force_recompute, cost_path, dtype, device
        )

        # Initialize the matrices `H`, which should be normalized.
        self.H = torch.rand(
            self.n_var, self.factors, device=device, dtype=dtype
        )
        self.H = utils.normalize_tensor(self.H)

        # Initialize the dual variable `G`
        self.G = torch.zeros_like(self.A, requires_grad=True)

        # Initialize the shared factor `W`, which should be normalized.
        self.W = torch.rand(self.factors, self.n_obs, device=device, dtype=dtype)
        self.W = utils.normalize_tensor(self.W)

        # Clean up.
        del keep_idx, features

    def train(
        self,
        adata_spatial,
        max_iter_inner: int = 1000,
        max_iter: int = 50,
        device: torch.device = "cpu",
        dtype: torch.dtype = torch.double,
        lr: float = 1e-2,
        optim_name: str = "adam",
        tol_inner: float = 1e-12,
        tol_outer: float = 1e-4,
        normalize_rows: bool = True,
        impute: bool = False,
        batch_size = 512,
    ) -> None:
        """Train the spotnmf model on an input AnnData object.

        Initializes the parameters, then alternately optimizes the usage
        ``W`` and the spectra ``H`` until convergence or ``max_iter`` is
        reached. On completion, the learned ``H`` is stored in
        ``adata_spatial.uns["H_OT"]`` and ``W`` in
        ``adata_spatial.obsm["W_OT"]``.

        Args:
            adata_spatial (anndata.AnnData): The input spatial AnnData object.
            max_iter_inner (int, optional): The maximum number of iterations
                for the inner optimization loop (optimizing H or W).
                Defaults to 1000.
            max_iter (int, optional): The maximum number of iterations for the
                outer optimization loop (successive optimizations of H and W).
                Defaults to 50.
            device (torch.device, optional): The device to work on.
                Defaults to "cpu".
            dtype (torch.dtype, optional): The dtype to work with.
                Defaults to torch.double.
            lr (float, optional): The learning rate for the optimizer. The
                default (1e-2) is tuned for adam; use ~1 for LBFGS.
                Defaults to 1e-2.
            optim_name (str, optional): The optimizer to use (``"adam"``,
                ``"lbfgs"`` or ``"sgd"``). Adam is the default and recommended:
                it optimizes the OT dual variable stably. LBFGS is *not* advised
                -- over many outer iterations it lets the dual variable drift
                toward zero, collapsing the usage matrix to a uniform per-spot
                mixture (no deconvolution signal). Defaults to "adam".
            tol_inner (float, optional): The tolerance for the inner iterations
                before early stopping. Defaults to 1e-12.
            tol_outer (float, optional): The tolerance for the outer iterations
                before early stopping. Defaults to 1e-4.
            normalize_rows (bool, optional): Whether to normalize the rows of
                the reference dataset during initialization. Defaults to True.
            impute (bool, optional): Whether to impute the data using FastICA
                during initialization. Defaults to False.
            batch_size (int, optional): The batch size. Defaults to 512.
        """

        # First, initialize the different parameters.
        self.init_parameters(
            adata_spatial,
            dtype=dtype,
            device=device,
            normalize_rows=normalize_rows,
            impute=impute
        )

        # This is needed to save things in uns if it doesn't exist.
        if adata_spatial.uns is None:
            adata_spatial.uns = {}

        self.lr = lr
        self.optim_name = optim_name

        # Initialize the loss histories.
        self.losses_w, self.losses_h, self.losses = [], [], []

        # Set up the progress bar.
        pbar = tqdm(total=2 * max_iter, position=0, leave=True)


        # This is the main loop, with at most `max_iter` iterations.
        try:
            for _ in range(max_iter):

                # Perform the `W` optimization step.
                self.optimize(
                    loss_fn=self.loss_fn_w,
                    max_iter=max_iter_inner,
                    tol=tol_inner,
                    history=self.losses_h,
                    pbar=pbar,
                    device=device,
                )

                # Update the shared factor `W`.
                htgw = self.H.T @ self.G
                coef = np.log(self.factors) / (self.w_regularization)
                self.W = F.softmin(coef * htgw.detach(), dim=0)
                # Clean up.
                del htgw

                # Update the progress bar.
                pbar.update(1)

                # Save the total dual loss and statistics.
                self.losses.append(self.total_dual_loss().cpu().detach())

                # Perform the `H` optimization step.
                self.optimize(
                    loss_fn=self.loss_fn_h,
                    device=device,
                    max_iter=max_iter_inner,
                    tol=tol_inner,
                    history=self.losses_h,
                    pbar=pbar,
                )

                # Update the omic specific factors `H`.
                coef = self.factors * np.log(self.n_var)
                coef /= self.n_obs * self.h_regularization

                self.H = self.G.detach()
                self.H = self.H @ self.W.T
                self.H = F.softmin(coef * self.H, dim=0)

                # Update the progress bar.
                pbar.update(1)

                # Save the total dual loss and statistics.
                self.losses.append(self.total_dual_loss().cpu().detach())

                # Early stopping
                if utils.early_stop(self.losses, tol_outer, nonincreasing=True):
                    break

        except KeyboardInterrupt:
            print("Training interrupted.")

        # Add H and W to the adata object.
        adata_spatial.uns["H_OT"] = self.H.cpu().numpy()
        adata_spatial.obsm["W_OT"] = self.W.T.cpu().numpy()

    def transform(
        self,
        adata_spatial,
        spectra,
        max_iter_inner: int = 1000,
        max_iter: int = 50,
        device: torch.device = "cpu",
        dtype: torch.dtype = torch.double,
        lr: float = 1e-2,
        optim_name: str = "adam",
        tol_inner: float = 1e-12,
        tol_outer: float = 1e-4,
        normalize_rows: bool = True,
        impute: bool = False,
    ) -> None:
        """Infer OT usages for fixed spectra (usage-only refit).

        This is the optimal-transport analogue of cNMF's NNLS usage refit
        (scikit-learn's ``non_negative_factorization(..., update_H=False)``):
        the spectra ``H`` are *frozen* to ``spectra`` -- e.g. the consensus
        programs obtained by clustering replicate factorizations -- and only
        the usage ``W`` is solved for. It runs the same entropic-OT dual ascent
        as :meth:`train`, but performs the ``W`` step only and never updates
        ``H``. On completion ``W`` is stored in ``adata_spatial.obsm["W_OT"]``
        and the (unchanged) ``H`` in ``adata_spatial.uns["H_OT"]``.

        Args:
            adata_spatial (anndata.AnnData): The input spatial AnnData object.
                Its highly-variable columns must align with the rows/columns of
                ``spectra`` (see below).
            spectra (numpy.ndarray or torch.Tensor): The fixed spectra, either
                as ``(genes x factors)`` or ``(factors x genes)``; it is
                oriented to ``(genes x factors)`` and column-normalized like a
                trained ``H``. ``genes`` must match the number of
                highly-variable features in ``adata_spatial``.
            max_iter_inner (int, optional): Max inner iterations for the ``W``
                (dual-variable) optimization. Defaults to 1000.
            max_iter (int, optional): Max outer iterations. Defaults to 50.
            device (torch.device, optional): Device to work on. Defaults "cpu".
            dtype (torch.dtype, optional): Working dtype. Defaults torch.double.
            lr (float, optional): Learning rate (tuned for adam). Defaults 1e-2.
            optim_name (str, optional): Optimizer; adam recommended (see
                :meth:`train`). Defaults "adam".
            tol_inner (float, optional): Inner early-stopping tolerance.
            tol_outer (float, optional): Outer early-stopping tolerance.
            normalize_rows (bool, optional): Row-normalize the reference dataset
                during initialization. Defaults True.
            impute (bool, optional): FastICA imputation during init. Defaults
                False.
        """

        # Initialize A, K, G (and throwaway random H/W) from the data.
        self.init_parameters(
            adata_spatial,
            dtype=dtype,
            device=device,
            normalize_rows=normalize_rows,
            impute=impute,
        )

        if adata_spatial.uns is None:
            adata_spatial.uns = {}

        # Freeze the spectra to the supplied (consensus) programs. Accept either
        # orientation and normalize columns the way a trained H is normalized.
        spectra = torch.as_tensor(np.asarray(spectra), dtype=dtype, device=device)
        if spectra.shape[0] != self.n_var:
            if spectra.shape[1] == self.n_var:
                spectra = spectra.T
            else:
                raise ValueError(
                    f"spectra shape {tuple(spectra.shape)} is incompatible with "
                    f"{self.n_var} highly-variable genes"
                )
        if spectra.shape[1] != self.factors:
            raise ValueError(
                f"spectra has {spectra.shape[1]} factors, expected {self.factors}"
            )
        # Guard against zero columns before column-normalization.
        spectra = spectra.clamp_min(0) + 1e-12
        self.H = utils.normalize_tensor(spectra)

        self.lr = lr
        self.optim_name = optim_name

        # Reset loss histories (a fresh solve).
        self.losses_w, self.losses_h, self.losses = [], [], []

        pbar = tqdm(total=max_iter, position=0, leave=True)

        try:
            for _ in range(max_iter):

                # Optimize the dual variable for the (frozen-H) W subproblem.
                self.optimize(
                    loss_fn=self.loss_fn_w,
                    max_iter=max_iter_inner,
                    tol=tol_inner,
                    history=self.losses_w,
                    pbar=pbar,
                    device=device,
                )

                # Update the usage W from the dual variable. H stays fixed.
                htgw = self.H.T @ self.G
                coef = np.log(self.factors) / (self.w_regularization)
                self.W = F.softmin(coef * htgw.detach(), dim=0)
                del htgw

                pbar.update(1)
                self.losses.append(self.total_dual_loss().cpu().detach())

                # Early stopping on the (W-only) outer objective.
                if utils.early_stop(self.losses, tol_outer, nonincreasing=True):
                    break

        except KeyboardInterrupt:
            print("Transform interrupted.")

        adata_spatial.uns["H_OT"] = self.H.cpu().numpy()
        adata_spatial.obsm["W_OT"] = self.W.T.cpu().numpy()

    def build_optimizer(
        self, params, lr: float, optim_name: str
    ) -> torch.optim.Optimizer:
        """Build the optimizer used to update the dual variable.

        The PyTorch LBFGS implementation is parametrized following the
        discussion in https://discuss.pytorch.org/t/unclear-purpose-of-max-iter-kwarg-in-the-lbfgs-optimizer/65695.

        Args:
            params (Iterable[torch.Tensor]): The parameters to be optimized.
            lr (float): The learning rate of the optimizer.
            optim_name (str): The name of the optimizer, among ``"lbfgs"``,
                ``"sgd"`` or ``"adam"``.

        Returns:
            torch.optim.Optimizer: The configured optimizer.
        """
        if optim_name == "lbfgs":
            return optim.LBFGS(
                params,
                lr=lr,
                history_size=5,
                max_iter=1,
                line_search_fn="strong_wolfe",
            )
        elif optim_name == "sgd":
            return optim.SGD(params, lr=lr)
        elif optim_name == "adam":
            return optim.Adam(params, lr=lr)

    def optimize(
        self,
        loss_fn: Callable,
        max_iter: int,
        history: List,
        tol: float,
        pbar,
        device: str,
    ) -> None:
        """Optimize the dual variable with respect to a given loss function.

        Runs the inner optimization loop over the dual variable ``G``,
        updating the progress bar and checking for early stopping every
        few steps.

        Args:
            loss_fn (Callable): The loss function to minimize.
            max_iter (int): The maximum number of iterations.
            history (list): A list to which the losses are appended.
            tol (float): The tolerance before early stopping.
            pbar (tqdm.tqdm): The progress bar to update.
            device (str): The device to work on.
        """

        # Build the optimizer.
        optimizer = self.build_optimizer(
            [self.G], lr=self.lr, optim_name=self.optim_name
        )

        # This value will be initially be displayed in the progress bar
        if len(self.losses) > 0:
            total_loss = self.losses[-1].cpu().numpy()
        else:
            total_loss = "?"

        # This is the main optimization loop.
        for i in range(max_iter):

            # Define the closure function required by the optimizer.
            def closure():
                optimizer.zero_grad()
                loss = loss_fn()
                loss.backward()
                return loss.detach()

            # Perform an optimization step.
            optimizer.step(closure)

            # Every x steps, update the progress bar.
            if i % 10 == 0:

                # Add a value to the loss history.
                history.append(loss_fn().cpu().detach())
                gpu_mem_alloc = torch.cuda.memory_allocated(device=device) if torch.cuda.is_available() else 0

                # Populate the progress bar.
                pbar.set_postfix(
                    {
                        "loss": total_loss,
                        "loss_inner": history[-1].cpu().numpy(),
                        "inner_steps": i,
                        "gpu_memory_allocated": gpu_mem_alloc,
                    }
                )

                # Attempt early stopping.
                if utils.early_stop(history, tol):
                    break

    @torch.no_grad()
    def total_dual_loss(self) -> torch.Tensor:
        """Compute the total dual loss.

        This is used for monitoring and early stopping only, not by the
        optimization algorithm itself. It combines the OT dual loss, the
        Lagrange multiplier term, and the entropy terms for ``H`` and ``W``.

        Returns:
            torch.Tensor: The total dual loss.
        """

        # Initialize the loss to zero.
        loss = 0

        # Add the OT dual loss.
        loss -= (
            utils.ot_dual_loss(
                self.A,
                self.G,
                self.K,
                self.eps
            )
            / self.n_obs
        )

        # Add the Lagrange multiplier term.
        lagrange = self.H @ self.W
        lagrange *= self.G
        lagrange = lagrange.sum()
        loss += lagrange / self.n_obs

        # Add the `H` entropy term.
        coef = self.h_regularization / (
            self.factors * np.log(self.n_var)
        )
        loss -= coef * utils.entropy(self.H, min_one=True)

        # Add the `W` entropy term.
        coef = (
            self.w_regularization / (self.n_obs * np.log(self.factors))
        )
        loss -= coef * utils.entropy(self.W, min_one=True)

        # Return the full loss.
        return loss

    def loss_fn_h(self) -> torch.Tensor:
        """Compute the loss for the update of the spectra ``H``.

        Returns:
            torch.Tensor: The loss for the ``H`` optimization step.
        """
        loss_h = 0

        # OT dual loss term
        loss_h += (
            utils.ot_dual_loss(
                self.A,
                self.G,
                self.K,
                self.eps
            )
            / self.n_obs
        )

        # Entropy dual loss term
        coef = self.h_regularization / (
            self.factors * np.log(self.n_var)
        )
        gwt = self.G @ self.W.T
        gwt /= self.n_obs * coef
        loss_h -= coef * utils.entropy_dual_loss(-gwt)

        # Clean up.
        del gwt

        # Return the loss.
        return loss_h

    def loss_fn_w(self) -> torch.Tensor:
        """Compute the loss for the update of the usage ``W``.

        Returns:
            torch.Tensor: The loss for the ``W`` optimization step.
        """
        loss_w, htgw = 0, 0

        # For the entropy dual loss term.
        htgw += self.H.T @ (self.G)

        # OT dual loss term.
        loss_w += (
            utils.ot_dual_loss(
                self.A,
                self.G,
                self.K,
                self.eps,
            )
            / self.n_obs
        )

        # Entropy dual loss term.
        coef = self.w_regularization
        coef /= self.n_obs * np.log(self.factors)
        htgw /= coef * self.n_obs
        loss_w -= coef * utils.entropy_dual_loss(-htgw)

        # Clean up.
        del htgw

        # Return the loss.
        return loss_w
    

def run_spotnmf(adata_spatial, components, seed=42, **kwargs):
    """Run the spotnmf model with flexible parameters provided via kwargs.

    Seeds the random number generators, builds and trains a
    :class:`spotnmf` model, then returns the learned factors as DataFrames.

    Args:
        adata_spatial (anndata.AnnData): Input spatial data as an AnnData
            object.
        components (int): The number of factors to learn.
        seed (int, optional): The random seed for reproducibility.
            Defaults to 42.
        **kwargs: Additional parameters for the model and training, all
            optional (sensible tuned defaults are provided): ``h`` (entropy
            regularization for the spectra, default 1e-2), ``w`` (entropy
            regularization for the usage, default 1e-2), ``eps`` (entropic
            regularization, default 2e-2), ``lr`` (learning rate, default
            1e-2), ``normalize_rows`` (default True), ``cost`` (default
            ``"cosine"``), ``optim_name`` (default ``"adam"``), plus
            ``tol_inner``, ``tol_outer``, ``max_iter``, ``max_iter_inner``
            and ``device``.

    Returns:
        tuple: A tuple ``(results, losses)`` where ``results`` is a dict with
        keys ``"topics_per_spot"`` and ``"genes_per_topic"`` mapping to
        pandas.DataFrame objects, and ``losses`` is the list of recorded
        total dual losses during training.
    """
    seed_all(seed)
    
    # Define default values. The regularization/optimizer defaults below were
    # tuned on the STARmap deconvolution benchmark (see scripts/benchmark/); they
    # can still be overridden via kwargs.
    defaults = {
        "h": 1e-2,             # entropic regularization on the spectra
        "w": 1e-2,             # entropic regularization on the usage
        "eps": 2e-2,           # entropic-OT smoothing
        "lr": 1e-2,            # optimizer learning rate (tuned for adam)
        "normalize_rows": True,
        "cost": "cosine",
        "optim_name": "adam",
        "tol_inner": 1e-12,
        "tol_outer": 0.00001,
        "max_iter": 100,
        "max_iter_inner": 1000,
        # Use the GPU when available, otherwise fall back to CPU.
        "device": "cuda" if torch.cuda.is_available() else "cpu",
    }
    defaults.update(kwargs)
    print("Model Params:", defaults)

    # Define the model with the selected parameters
    model = spotnmf(
        factors=int(components),
        h_regularization=defaults["h"],
        w_regularization=defaults["w"],
        eps=defaults["eps"],
        cost=defaults["cost"],
        pca_cost=False,
    )

    # Train the model with selected parameters
    model.train(
        adata_spatial,
        lr=defaults["lr"],
        optim_name= defaults["optim_name"],
        tol_inner=defaults["tol_inner"],
        tol_outer=defaults["tol_outer"],
        normalize_rows=defaults["normalize_rows"],
        max_iter=defaults["max_iter"],
        max_iter_inner=defaults["max_iter_inner"],
        impute=False,
        device=defaults["device"],
    )

    # Extract results
    gene_list = adata_spatial.var[adata_spatial.var.highly_variable].index
    df_genes_per_topic = pd.DataFrame(adata_spatial.uns["H_OT"], index=gene_list)
    df_topics_per_spot = pd.DataFrame(adata_spatial.obsm["W_OT"], index=adata_spatial.obs.index)

    results = {
        "topics_per_spot": df_topics_per_spot,
        "genes_per_topic": df_genes_per_topic
    }

    # Format column names
    for key_matrix in results:
        results[key_matrix].columns = [f"ot_{x+1}" for x in results[key_matrix].columns]

    return results, model.losses


def _cluster_consensus_spectra(spectra_pool, k, density_threshold, local_neighborhood_size):
    """Cluster pooled replicate spectra into ``k`` consensus programs.

    Mirrors cNMF's consensus clustering (:meth:`spotnmf.external.cnmf.cNMF.consensus`):
    L2-normalize every replicate program, optionally drop local-density
    outliers, KMeans into ``k`` clusters, and take the per-cluster median as the
    consensus spectrum (renormalized to a probability distribution).

    Args:
        spectra_pool (pandas.DataFrame): ``(n_iter*k) x genes`` pooled replicate
            spectra (each row is one program from one replicate).
        k (int): Number of consensus programs.
        density_threshold (float): Local-density outlier threshold (>=2 disables
            filtering), matching cNMF's ``local_density_threshold``.
        local_neighborhood_size (float): Fraction of replicate programs used as
            neighbors for the density estimate, matching cNMF.

    Returns:
        tuple: ``(median_spectra, stability)`` where ``median_spectra`` is a
        ``k x genes`` DataFrame and ``stability`` is the KMeans silhouette score
        (``nan`` if it cannot be computed).
    """
    from sklearn.cluster import KMeans
    from sklearn.metrics import silhouette_score
    from sklearn.metrics.pairwise import euclidean_distances

    # Rescale every replicate program to unit L2 norm (cNMF's l2_spectra).
    l2 = spectra_pool.div(np.sqrt((spectra_pool ** 2).sum(axis=1)), axis=0)

    # Local-density outlier filtering (skipped when threshold >= 2).
    n_iter_est = max(1, l2.shape[0] // k)
    n_neighbors = int(local_neighborhood_size * l2.shape[0] / k)
    if density_threshold < 2 and n_neighbors >= 1 and l2.shape[0] > n_neighbors + 1:
        dists = euclidean_distances(l2.values)
        part = np.argpartition(dists, n_neighbors + 1)[:, : n_neighbors + 1]
        nn = dists[np.arange(dists.shape[0])[:, None], part]
        local_density = nn.sum(1) / n_neighbors
        keep = local_density < density_threshold
        if keep.sum() >= k:
            l2 = l2.loc[keep]

    km = KMeans(n_clusters=k, n_init=10, random_state=1)
    labels = km.fit_predict(l2.values)
    labels = pd.Series(labels, index=l2.index)

    stability = float("nan")
    try:
        stability = float(silhouette_score(l2.values, labels, metric="euclidean"))
    except Exception:
        pass

    median_spectra = l2.groupby(labels).median()
    median_spectra = median_spectra.div(median_spectra.sum(axis=1), axis=0)
    return median_spectra, stability


def run_spotnmf_consensus(
    adata_spatial,
    components,
    n_iter: int = 10,
    seed: int = 42,
    density_threshold: float = 0.5,
    local_neighborhood_size: float = 0.30,
    return_stability: bool = False,
    **kwargs,
):
    """Consensus spOT-NMF: cluster replicate OT spectra, then OT-refit usages.

    This is spOT-NMF's proper consensus, structured exactly like cNMF
    (:func:`spotnmf.external.cnmf.run_cnmf`) with one deliberate difference: the
    final usage refit uses spOT-NMF's *fixed-spectra optimal-transport* solve
    (:meth:`spotnmf.transform`) instead of cNMF's NNLS. Concretely:

    1. Run ``n_iter`` independent spOT-NMF factorizations (seeds ``seed`` ..
       ``seed + n_iter - 1``), collecting each replicate's spectra.
    2. Pool and cluster the replicate spectra into ``components`` consensus
       programs (L2-normalize, density-filter, KMeans, per-cluster median) --
       identical to cNMF.
    3. Freeze those consensus programs and solve only the usages with the
       entropic-OT dual ascent, so the consensus usages keep the OT geometry
       rather than reverting to least squares.

    Args:
        adata_spatial (anndata.AnnData): Input spatial data.
        components (int): Number of programs/factors ``k``.
        n_iter (int, optional): Number of replicate factorizations. Defaults 10.
        seed (int, optional): Base random seed. Defaults 42.
        density_threshold (float, optional): cNMF local-density outlier
            threshold (>=2 disables filtering). Defaults 0.5.
        local_neighborhood_size (float, optional): cNMF neighbor fraction.
            Defaults 0.30.
        return_stability (bool, optional): If True, also return the KMeans
            silhouette stability. Defaults False.
        **kwargs: Forwarded to :func:`run_spotnmf` / :meth:`spotnmf.transform`
            (``h``, ``w``, ``eps``, ``lr``, ``normalize_rows``, ``cost``,
            ``optim_name``, ``max_iter``, ``max_iter_inner``, ``device`` ...);
            the same tuned defaults apply.

    Returns:
        tuple: ``(results, losses)`` -- or ``(results, losses, stability)`` when
        ``return_stability`` is True -- where ``results`` has the same
        ``"topics_per_spot"`` / ``"genes_per_topic"`` DataFrames as
        :func:`run_spotnmf`.
    """
    k = int(components)

    # Resolve the same tuned defaults used by run_spotnmf / transform.
    defaults = {
        "h": 1e-2, "w": 1e-2, "eps": 2e-2, "lr": 1e-2,
        "normalize_rows": True, "cost": "cosine", "optim_name": "adam",
        "tol_inner": 1e-12, "tol_outer": 0.00001,
        "max_iter": 100, "max_iter_inner": 1000,
        "device": "cuda" if torch.cuda.is_available() else "cpu",
    }
    defaults.update(kwargs)

    gene_list = adata_spatial.var[adata_spatial.var.highly_variable].index

    # 1. Replicate factorizations (each on a fresh copy so uns/obsm aren't shared).
    spectra_frames = []
    for i in range(n_iter):
        rep = adata_spatial.copy()
        res, _ = run_spotnmf(rep, components=k, seed=seed + i, **kwargs)
        # genes_per_topic is (genes x k); transpose to (k x genes) programs.
        prog = res["genes_per_topic"].T
        prog.index = [f"iter{i}_topic{t + 1}" for t in range(k)]
        prog.columns = gene_list
        spectra_frames.append(prog)
    spectra_pool = pd.concat(spectra_frames, axis=0)

    # 2. Cluster into consensus programs (cNMF-style).
    median_spectra, stability = _cluster_consensus_spectra(
        spectra_pool, k, density_threshold, local_neighborhood_size
    )

    # 3. OT usage refit with the consensus spectra held fixed.
    seed_all(seed)
    model = spotnmf(
        factors=k,
        h_regularization=defaults["h"],
        w_regularization=defaults["w"],
        eps=defaults["eps"],
        cost=defaults["cost"],
        pca_cost=False,
    )
    refit = adata_spatial.copy()
    model.transform(
        refit,
        spectra=median_spectra.reindex(columns=gene_list).to_numpy(),  # (k x genes)
        lr=defaults["lr"],
        optim_name=defaults["optim_name"],
        tol_inner=defaults["tol_inner"],
        tol_outer=defaults["tol_outer"],
        normalize_rows=defaults["normalize_rows"],
        max_iter=defaults["max_iter"],
        max_iter_inner=defaults["max_iter_inner"],
        impute=False,
        device=defaults["device"],
    )

    df_topics_per_spot = pd.DataFrame(refit.obsm["W_OT"], index=adata_spatial.obs.index)
    df_genes_per_topic = pd.DataFrame(refit.uns["H_OT"], index=gene_list)
    results = {
        "topics_per_spot": df_topics_per_spot,
        "genes_per_topic": df_genes_per_topic,
    }
    for key_matrix in results:
        results[key_matrix].columns = [f"ot_{x + 1}" for x in results[key_matrix].columns]

    if return_stability:
        return results, model.losses, stability
    return results, model.losses
