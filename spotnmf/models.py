"""
This file implements a modified version of the Mowgli model, transforming it into a 
1 mode of Optimal Transport Non-negative Matrix Factorization (OT-NMF) approach.
The original Mowgli model was developed for paired single-cell multi-omics data integration, 
as described in the publication: Huizing, G.-J., Deutschmann, I. M., Peyré, G., & Cantini, L. (2023). 
Paired single-cell multi-omics data integration with Mowgli. Nature Communications, 14(1), 7711.
"""

from typing import Callable
import numpy as np
import torch
import torch.nn.functional as F
from torch import optim
from tqdm import tqdm
import pandas as pd 
from spotnmf import utils

def _resolve_dtype(dtype):
    """Map a dtype spec to a torch.dtype.

    Accepts a torch.dtype directly, or a string alias so that callers going
    through ``run_spotnmf(..., dtype="float32")`` don't need to import torch.
    The whole OT solve is log-sum-exp stabilized, so ``"float32"`` is safe and,
    on consumer GPUs where FP64 runs at a fraction of FP32 throughput, much
    faster and half the memory. Defaults elsewhere remain ``"float64"`` so
    behaviour is unchanged unless explicitly requested.
    """
    if isinstance(dtype, torch.dtype):
        return dtype
    aliases = {
        "float64": torch.float64, "double": torch.float64, "fp64": torch.float64,
        "float32": torch.float32, "float": torch.float32, "fp32": torch.float32,
    }
    key = str(dtype).lower()
    if key not in aliases:
        raise ValueError(f"unknown dtype {dtype!r}; use float32 or float64")
    return aliases[key]


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

        # Initialize the (outer) loss history.
        self.losses = []

        self.A, self.H, self.G, self.K = None, None, None, None

    def init_parameters(
        self,
        adata_spatial,
        dtype: torch.dtype,
        device: torch.device,
        force_recompute: bool = False,
        normalize_rows: bool = False,
        impute: bool = False,
        init: str = "random",
        row_power: float = 1.0,
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
            # Give every gene the same total mass across spots before the column
            # normalization below (so no gene dominates the OT reference). This
            # is the single-modality case: the original Mowgli code also scaled by
            # ``mean_row_sum`` to balance *across modalities*, but with one
            # modality that uniform factor cancels exactly under the subsequent
            # column normalization, so it is dropped (and avoids dividing by a
            # possibly-large scalar, which is slightly safer in float32).
            # ``row_power`` in (0, 1] interpolates between no row scaling (0) and
            # full gene-mass equalization (1, the default).
            self.A /= self.A.sum(1).reshape(-1, 1) ** row_power
        self.A /= self.A.sum(0)

        # The cost metric and any cached-cost path are set in __init__.
        cost = self.cost
        cost_path = self.cost_path

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

        # Initialize the spectra `H` (genes x k) and usage `W` (k x spots).
        # ``init="nndsvd"`` seeds them from NNDSVD (the deterministic SVD-based
        # NMF initialization); random NMF starts have high run-to-run variance and
        # land in poorer local optima, which is much of why the consensus helps.
        # Both are column-normalized onto the probability simplex like a trained
        # H/W. ``init="random"`` (default) preserves the original behaviour.
        if init == "nndsvd":
            from sklearn.decomposition._nmf import _initialize_nmf

            Xnn = X_data.toarray() if hasattr(X_data, "toarray") else np.asarray(X_data)
            Xnn = np.clip(np.asarray(Xnn, dtype=np.float64), 0, None)
            W0, H0 = _initialize_nmf(Xnn, self.factors, init="nndsvda", random_state=0)
            self.H = utils.normalize_tensor(
                torch.as_tensor(H0.T, device=device, dtype=dtype).clamp_min(1e-12))
            self.W = utils.normalize_tensor(
                torch.as_tensor(W0.T, device=device, dtype=dtype).clamp_min(1e-12))
        else:
            self.H = torch.rand(self.n_var, self.factors, device=device, dtype=dtype)
            self.H = utils.normalize_tensor(self.H)
            self.W = torch.rand(self.factors, self.n_obs, device=device, dtype=dtype)
            self.W = utils.normalize_tensor(self.W)

        # Initialize the dual variable `G`.
        self.G = torch.zeros_like(self.A, requires_grad=True)

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
        init: str = "random",
        damping: float = 1.0,
        persist_optimizer: bool = False,
        w_anneal: float = 1.0,
        analytic_grad: bool = None,
        row_power: float = 1.0,
        batch_size: int = None,
        g_tail: int = 0,
        readout_temp: float = 1.0,
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
        """

        # First, initialize the different parameters.
        self.init_parameters(
            adata_spatial,
            dtype=dtype,
            device=device,
            normalize_rows=normalize_rows,
            impute=impute,
            init=init,
            row_power=row_power,
        )

        # This is needed to save things in uns if it doesn't exist.
        if adata_spatial.uns is None:
            adata_spatial.uns = {}

        self.lr = lr
        self.optim_name = optim_name

        # Analytic (closed-form) dual gradient vs autograd. It is the same
        # gradient (to ~1e-14 in float64) but skips the autograd backward, whose
        # per-step overhead dominates in float32 (1.8-2.2x faster). In float64 the
        # step is matmul-bound, so it is no faster (occasionally slower on large
        # gene panels). Default (``None``) therefore auto-selects: on for float32,
        # off for float64.
        if analytic_grad is None:
            analytic_grad = dtype == torch.float32
        grad_w = self._grad_w if analytic_grad else None
        grad_h = self._grad_h if analytic_grad else None

        # Optional persistent dual optimizer: reuse one optimizer (and its Adam
        # moment estimates) across every inner solve instead of rebuilding it each
        # block. Because G is warm-started, rebuilding throws away momentum and
        # re-warms every block; persisting can converge the dual with fewer steps.
        self.persist_optimizer = persist_optimizer
        self._persist_opt = None

        # Entropy/temperature annealing on the usage regularizer w. The W update
        # temperature is coef = log(k)/w, which is very aggressive for small w and
        # can make an over-converged dual collapse spots onto a shared program.
        # With w_anneal > 1 we start at w*w_anneal (gentler, explore broad
        # structure) and linearly decay to the target w (sharpen, commit late).
        w_target = self.w_regularization

        # Initialize the (outer) loss history.
        self.losses = []

        # Set up the progress bar.
        pbar = tqdm(total=2 * max_iter, position=0, leave=True)


        # This is the main loop, with at most `max_iter` iterations.
        try:
            for t in range(max_iter):

                # Anneal the usage temperature toward the target w.
                if w_anneal != 1.0 and max_iter > 1:
                    frac = t / (max_iter - 1)
                    self.w_regularization = w_target * (1 + (w_anneal - 1) * (1 - frac))

                # Perform the `W` optimization step.
                self.optimize(
                    loss_fn=self.loss_fn_w,
                    max_iter=max_iter_inner,
                    tol=tol_inner,
                    pbar=pbar,
                    grad_fn=grad_w,
                    batch_size=batch_size,
                    g_tail=g_tail,
                )

                # Update the shared factor `W` (optionally damped). Damping
                # relaxes the closed-form softmin update toward the previous value,
                # H_new -> (1-d)*H_old + d*H_new; since both are column-simplex
                # distributions the convex combination stays valid. This stabilizes
                # the outer alternation: the peaked softmin (coef ~ log(k)/w) can
                # otherwise let an over-converged dual collapse spots onto a shared
                # program. damping=1.0 recovers the original exact update.
                W_new = self._usage_from_dual()
                self.W = W_new if damping >= 1.0 else (1 - damping) * self.W + damping * W_new

                # Update the progress bar.
                pbar.update(1)

                # Save the total dual loss and statistics.
                self.losses.append(self.total_dual_loss().cpu().detach())

                # Perform the `H` optimization step.
                self.optimize(
                    loss_fn=self.loss_fn_h,
                    max_iter=max_iter_inner,
                    tol=tol_inner,
                    pbar=pbar,
                    grad_fn=grad_h,
                    batch_size=batch_size,
                    g_tail=g_tail,
                )

                # Update the omic specific factors `H`.
                coef = self.factors * np.log(self.n_var)
                coef /= self.n_obs * self.h_regularization

                H_new = F.softmin(coef * (self.G.detach() @ self.W.T), dim=0)
                self.H = H_new if damping >= 1.0 else (1 - damping) * self.H + damping * H_new

                # Update the progress bar.
                pbar.update(1)

                # Save the total dual loss and statistics.
                self.losses.append(self.total_dual_loss().cpu().detach())

                # Early stopping
                if utils.early_stop(self.losses, tol_outer, nonincreasing=True):
                    break

        except KeyboardInterrupt:
            print("Training interrupted.")

        # Read out the FINAL usages from the converged dual and spectra, rescaling
        # the softmin inverse temperature by ``readout_temp`` (tau). This recompute
        # (rather than keeping the loop's last W, which was formed before the final
        # H update) keeps W consistent with the stored H and makes tau a continuous
        # knob through 1.0. That single temperature is best for the *optimization*
        # but not necessarily the deconvolution metric: tau>1 sharpens (crisp
        # panels, also lowers seed variance), tau<1 softens toward graded
        # proportions (mixed spots, e.g. MOB); tau=1.0 is the plain softmin readout.
        self.W = self._usage_from_dual(readout_temp)

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
        analytic_grad: bool = None,
        readout_temp: float = 1.0,
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

        # Reset the (outer) loss history (a fresh solve).
        self.losses = []

        # Auto-select the analytic gradient for float32 (see :meth:`train`).
        if analytic_grad is None:
            analytic_grad = dtype == torch.float32
        grad_w = self._grad_w if analytic_grad else None

        pbar = tqdm(total=max_iter, position=0, leave=True)

        try:
            for _ in range(max_iter):

                # Optimize the dual variable for the (frozen-H) W subproblem.
                self.optimize(
                    loss_fn=self.loss_fn_w,
                    max_iter=max_iter_inner,
                    tol=tol_inner,
                    pbar=pbar,
                    grad_fn=grad_w,
                )

                # Update the usage W from the dual variable. H stays fixed.
                self.W = self._usage_from_dual()

                pbar.update(1)
                self.losses.append(self.total_dual_loss().cpu().detach())

                # Early stopping on the (W-only) outer objective.
                if utils.early_stop(self.losses, tol_outer, nonincreasing=True):
                    break

        except KeyboardInterrupt:
            print("Transform interrupted.")

        # Read out the final (frozen-H) usages, rescaled by ``readout_temp``; see
        # :meth:`train`. Here the W step is the last loop op, so tau=1.0 reproduces
        # the loop's usages exactly.
        self.W = self._usage_from_dual(readout_temp)

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
        elif optim_name == "adamw":
            return optim.AdamW(params, lr=lr)
        elif optim_name == "nadam":
            return optim.NAdam(params, lr=lr)
        elif optim_name == "rmsprop":
            return optim.RMSprop(params, lr=lr)
        raise ValueError(f"unknown optim_name {optim_name!r}")

    def optimize(
        self,
        loss_fn: Callable,
        max_iter: int,
        tol: float,
        pbar,
        grad_fn: Callable = None,
        batch_size: int = None,
        g_tail: int = 0,
    ) -> None:
        """Optimize the dual variable with respect to a given loss function.

        Runs the inner optimization loop over the dual variable ``G``,
        updating the progress bar and checking for early stopping every
        few steps.

        Args:
            loss_fn (Callable): The loss function to minimize (also used to
                report the loss for monitoring/early-stopping).
            max_iter (int): The maximum number of iterations.
            tol (float): The tolerance before early stopping.
            pbar (tqdm.tqdm): The progress bar to update.
            grad_fn (Callable, optional): If given, a function returning the
                closed-form gradient of ``loss_fn`` w.r.t. ``G``. It is assigned
                to ``G.grad`` directly, skipping autograd's forward+backward
                (the dominant per-step cost). ``loss_fn`` is then evaluated only
                at the monitoring cadence. ``None`` uses the autograd path.
            g_tail (int, optional): If >0, replace ``G`` on exit with the mean of
                the last ``g_tail`` inner iterates (Polyak/tail averaging), to
                cancel the dual's near-optimum jitter before the peaked-softmin
                block update. Defaults 0 (off); benchmarked as a no-op (see the
                ``g_tail`` sweep), kept only as an option.
        """

        # Build the optimizer (or reuse a persistent one across inner solves).
        if getattr(self, "persist_optimizer", False):
            if getattr(self, "_persist_opt", None) is None:
                self._persist_opt = self.build_optimizer(
                    [self.G], lr=self.lr, optim_name=self.optim_name)
            optimizer = self._persist_opt
        else:
            optimizer = self.build_optimizer(
                [self.G], lr=self.lr, optim_name=self.optim_name
            )

        # This value will be initially be displayed in the progress bar
        if len(self.losses) > 0:
            total_loss = self.losses[-1].cpu().numpy()
        else:
            total_loss = "?"

        # The closure is defined once (nothing it closes over changes across
        # steps) rather than rebuilt every iteration. ``optimizer.step`` returns
        # this closure's (detached) loss, so we reuse that value for monitoring
        # and early stopping instead of a second forward pass through ``loss_fn``.
        def closure():
            optimizer.zero_grad()
            loss = loss_fn()
            loss.backward()
            return loss.detach()

        # Fresh convergence history for THIS inner solve. The early-stopping
        # decision must depend only on the current subproblem's progress -- a
        # persistent list would compare against stale losses from a previous
        # (different) subproblem, making the check meaningless.
        local_hist = []

        # Tail-averaging accumulator (Polyak): running sum of the last g_tail
        # iterates of G, kept as a single buffer rather than g_tail copies.
        g_sum, g_n = None, 0
        tail_start = max(0, max_iter - g_tail) if g_tail > 0 else max_iter

        # This is the main optimization loop.
        for i in range(max_iter):

            if grad_fn is None:
                # Autograd path: step returns the loss the closure just evaluated.
                loss = optimizer.step(closure)
            else:
                # Analytic path: set G.grad directly and step; the loss value is
                # only needed at the monitoring cadence below. With batch_size a
                # random spot-column subset is updated each step (block-coordinate
                # / stochastic dual ascent).
                if batch_size is not None and batch_size < self.n_obs:
                    idx = torch.randperm(self.n_obs, device=self.G.device)[:batch_size]
                    self.G.grad = grad_fn(idx)
                else:
                    self.G.grad = grad_fn()
                optimizer.step()

            # Accumulate the tail iterates of G for Polyak averaging.
            if i >= tail_start:
                g_det = self.G.detach()
                if g_sum is None:
                    g_sum = g_det.clone()
                else:
                    g_sum += g_det
                g_n += 1

            # Every x steps, check convergence and update the progress bar.
            if i % 10 == 0:
                loss = loss_fn().detach() if grad_fn is not None else loss
                local_hist.append(loss)

                # Attempt early stopping (based on this solve's history only).
                if utils.early_stop(local_hist, tol):
                    break

                # Populate the progress bar (one host sync per 10 steps).
                pbar.set_postfix(
                    {"loss": total_loss, "loss_inner": loss.item(), "inner_steps": i}
                )

        # Replace G with the tail average (Polyak) before the block update reads
        # it through the peaked softmin. Done in-place so warm-starting and the
        # W/H softmin updates see the denoised dual.
        if g_n:
            self.G.data.copy_(g_sum.div_(g_n))
            # The dual just jumped to the tail average; a persistent optimizer's
            # Adam moments still describe the pre-average trajectory, so clear them
            # to avoid mis-scaled steps on the next inner solve.
            if getattr(self, "persist_optimizer", False):
                optimizer.state.clear()

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

    def _usage_from_dual(self, temp: float = 1.0) -> torch.Tensor:
        """Usage ``W`` read out from the current dual ``G`` and spectra ``H``.

        ``W = softmin(temp * (log k / w) * HᵀG)`` over the ``k`` factors -- the
        closed-form entropic-usage update shared by :meth:`train` and
        :meth:`transform`. ``temp`` (the readout temperature) rescales the softmin
        inverse temperature ``log(k)/w``: ``temp=1.0`` is the in-loop update;
        ``temp>1`` sharpens the final usages (crisp panels, also lowers seed
        variance), ``temp<1`` softens toward graded proportions (mixed spots).
        """
        coef = temp * np.log(self.factors) / self.w_regularization
        return F.softmin(coef * (self.H.T @ self.G).detach(), dim=0)

    @torch.no_grad()
    def _grad_w(self, idx=None) -> torch.Tensor:
        """Closed-form gradient of :meth:`loss_fn_w` w.r.t. the dual ``G``.

        Analytic alternative to autograd's backward pass. The entropic-OT dual
        term has a Sinkhorn-marginal gradient ``v ⊙ (K @ (A / (K @ v)))`` with
        ``v = exp(G/eps)`` (``K`` is symmetric), and the usage-entropy term's
        gradient is exactly the ``softmin`` used for the ``W`` update. Verified
        to match autograd to ~1e-14 (float64); it just avoids building/traversing
        the autograd graph, which is the dominant per-step cost in float32.

        The W subproblem is separable over spots, so ``idx`` (a spot-column
        subset) yields the *exact* gradient restricted to those columns (a
        random block-coordinate step); other columns are left zero.
        """
        coef = np.log(self.factors) / self.w_regularization
        if idx is None:
            grad = self._ot_grad() - self.H @ F.softmin(coef * (self.H.T @ self.G), dim=0)
            return grad / self.n_obs
        Gb = self.G[:, idx]
        grad_b = self._ot_grad(idx) - self.H @ F.softmin(coef * (self.H.T @ Gb), dim=0)
        full = torch.zeros_like(self.G)
        full[:, idx] = grad_b / self.n_obs
        return full

    @torch.no_grad()
    def _grad_h(self, idx=None) -> torch.Tensor:
        """Closed-form gradient of :meth:`loss_fn_h` w.r.t. the dual ``G``.

        Same OT term as :meth:`_grad_w`; the spectra-entropy term's gradient is
        the ``softmin`` used for the ``H`` update. Verified against autograd to
        ~1e-14 (float64). See :meth:`_grad_w`. With ``idx`` the OT term is
        computed on the spot-column subset while the (spot-coupled) entropy term
        is evaluated on the full ``G`` and sliced -- still the exact partial
        gradient for those columns.
        """
        coef = self.factors * np.log(self.n_var) / (self.n_obs * self.h_regularization)
        grad_ent = F.softmin(coef * (self.G @ self.W.T), dim=0) @ self.W
        if idx is None:
            return (self._ot_grad() - grad_ent) / self.n_obs
        full = torch.zeros_like(self.G)
        full[:, idx] = (self._ot_grad(idx) - grad_ent[:, idx]) / self.n_obs
        return full

    @torch.no_grad()
    def _ot_grad(self, idx=None) -> torch.Tensor:
        """Gradient of the entropic-OT dual term w.r.t. ``G`` (shared by W/H).

        ``grad = v ⊙ (K @ (A / (K @ v)))`` with ``v = exp(G/eps)``. Computed with
        the same per-column log-sum-exp stabilization as :func:`utils.ot_dual_loss`
        (subtract the per-column max of ``G/eps`` before exp); that constant
        cancels exactly in the ratio, so the result is unchanged while ``exp`` no
        longer overflows in float32. The OT term is separable over spots, so
        ``idx`` restricts the computation to a subset of spot columns (returning
        ``genes x |idx|``).
        """
        G = self.G if idx is None else self.G[:, idx]
        A = self.A if idx is None else self.A[:, idx]
        log_v = G / self.eps
        log_v = log_v - log_v.max(0).values
        v = torch.exp(log_v)
        return v * (self.K @ (A / (self.K @ v)))


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
        "init": "random",
        "damping": 1.0,
        "persist_optimizer": False,
        "w_anneal": 1.0,
        "row_power": 1.0,
        # Stochastic spot mini-batching in the (analytic) inner dual ascent.
        # None = full batch (all spots each step). Only active with analytic_grad.
        "batch_size": None,
        # Polyak tail-averaging of the dual G: average the last g_tail inner
        # iterates before the peaked-softmin block update, to cancel near-optimum
        # dual jitter (within-run analogue of consensus). 0 = off.
        "g_tail": 0,
        # Readout-temperature (tau) for the FINAL stored usages only: rescales the
        # softmin inverse temperature log(k)/w. >1 sharpens (crisp panels, also
        # lowers seed variance), <1 softens toward graded proportions (mixed spots
        # e.g. MOB). 1.0 = usages exactly as trained. Per-dataset, like w.
        "readout_temp": 1.0,
        # Closed-form dual gradient instead of autograd (identical gradient to
        # ~1e-14; skips the backward pass). None auto-selects: on for float32
        # (1.8-2.2x faster), off for float64 (matmul-bound, no gain). True/False
        # force it.
        "analytic_grad": None,
        # Working precision. "float64" preserves the original behaviour; "float32"
        # halves memory (G, A, K and Adam's two G-moment buffers) and is much
        # faster on GPUs with limited FP64 throughput. The OT solve is
        # log-sum-exp stabilized, so float32 is numerically safe here.
        "dtype": "float64",
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
        init=defaults["init"],
        damping=defaults["damping"],
        persist_optimizer=defaults["persist_optimizer"],
        w_anneal=defaults["w_anneal"],
        analytic_grad=defaults["analytic_grad"],
        row_power=defaults["row_power"],
        batch_size=defaults["batch_size"],
        g_tail=defaults["g_tail"],
        readout_temp=defaults["readout_temp"],
        dtype=_resolve_dtype(defaults["dtype"]),
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


def _l1_normalize_rows(M):
    """Rescale each row of a non-negative array to sum to 1 (a distribution)."""
    return M / np.clip(M.sum(axis=1, keepdims=True), 1e-12, None)


def _gene_ground_cost(adata_spatial, gene_list, cost, normalize_rows):
    """Genes x genes ground cost on the model's normalized reference.

    Reproduces :meth:`spotnmf.init_parameters`'s reference construction (row/column
    normalization of the highly-variable expression matrix) and computes the
    pairwise gene distance with the model's ``cost`` metric, so the Wasserstein
    barycenter transports mass under the same geometry the factorization uses.
    """
    from scipy.spatial.distance import cdist

    X = adata_spatial[:, adata_spatial.var["highly_variable"]].X
    X = X.toarray() if hasattr(X, "toarray") else np.asarray(X)
    A = np.asarray(X, dtype=np.float64).T.copy()      # genes x spots
    A += 1e-6
    if normalize_rows:
        mean_row_sum = A.sum(1).mean()
        A /= A.sum(1, keepdims=True) * mean_row_sum
    A /= A.sum(0, keepdims=True)
    features = 1e-6 + A
    metric = cost if isinstance(cost, str) else "cosine"
    if metric == "ones":
        return 1.0 - np.eye(features.shape[0])
    return cdist(features, features, metric=metric)


def _wasserstein_barycenter(dists, kernel, weights=None,
                            max_iter=1000, tol=1e-5):
    """Entropic Wasserstein barycenter of gene-distributions (Benamou et al. IBP).

    Computes the fixed-support entropic Wasserstein barycenter of the columns of
    ``dists`` (each a distribution over the same gene support) using Iterative
    Bregman Projections with the Gibbs ``kernel = exp(-M/reg)`` derived from the
    gene ground cost ``M``. Because the barycenter transports mass between similar
    genes (per the OT geometry the model itself uses), it denoises the consensus
    program in an OT-native way that a coordinate-wise mean/geomean cannot.

    Args:
        dists (numpy.ndarray): ``(genes, S)`` non-negative columns (each summed to
            one) -- the ``S`` matched replicate programs for one consensus cluster.
        kernel (numpy.ndarray): ``(genes, genes)`` Gibbs kernel ``exp(-M/reg)``.
        weights (numpy.ndarray, optional): length-``S`` barycentric weights
            (uniform by default).
        max_iter (int, optional): Max IBP iterations. Defaults 1000.
        tol (float, optional): Convergence tolerance on the barycenter. Defaults 1e-5.

    Returns:
        numpy.ndarray: ``(genes,)`` barycenter distribution (sums to one).
    """
    n, S = dists.shape
    if weights is None:
        weights = np.ones(S) / S
    u = np.ones((n, S)) / n
    UKv = u * (kernel @ (dists / np.clip(kernel @ u, 1e-300, None)))
    bar = np.exp(np.log(np.clip(UKv, 1e-300, None)) @ weights)
    for cpt in range(max_iter):
        u = (bar[:, None] / np.clip(UKv, 1e-300, None)) * u
        UKv = u * (kernel @ (dists / np.clip(kernel @ u, 1e-300, None)))
        bar_new = np.exp(np.log(np.clip(UKv, 1e-300, None)) @ weights)
        if cpt % 10 == 0 and np.abs(bar_new - bar).sum() < tol:
            bar = bar_new
            break
        bar = bar_new
    s = bar.sum()
    return bar / s if s > 0 else np.ones(n) / n


def _bary_kernel(cost_matrix, bary_reg):
    """Gibbs kernel ``exp(-M/reg)`` from a (max-normalized) gene ground cost."""
    M = np.asarray(cost_matrix, dtype=np.float64)
    M = M / (M.max() + 1e-12)
    return np.exp(-M / bary_reg)


def _aggregate_programs(stacked, aggregate, cost_matrix=None, bary_reg=0.1):
    """Aggregate matched replicate programs ``(n_iter, k, genes)`` over replicates.

    spOT-NMF spectra are probability distributions over genes, so the natural
    consensus is a *distribution barycenter*, not cNMF's coordinate-wise median:

    * ``"wasserstein"`` -- entropic **Wasserstein barycenter** under the gene
      ground cost (the fully OT-native aggregation; requires ``cost_matrix``).
      Best on crisp/mixed panels in the spOT-NMF benchmark and lets the safe
      ``match`` clustering reach the accuracy that otherwise needed ``kmeans``.
    * ``"geomean"`` -- normalized geometric mean, i.e. the **KL barycenter** of
      the distributions. Best on low-k mixed data (e.g. MOB).
    * ``"mean"`` -- arithmetic mean (the mixture barycenter); robust on PCC.
    * ``"median"`` -- cNMF's choice; kept for comparison. Over-sparsifies
      distribution-valued spectra and can wipe out the consensus benefit.
    """
    if aggregate == "median":
        return np.median(stacked, axis=0)
    if aggregate == "mean":
        return stacked.mean(axis=0)
    if aggregate == "geomean":
        return np.exp(np.log(np.clip(stacked, 1e-12, None)).mean(axis=0))
    if aggregate == "wasserstein":
        if cost_matrix is None:
            raise ValueError("aggregate='wasserstein' requires cost_matrix")
        kernel = _bary_kernel(cost_matrix, bary_reg)
        k = stacked.shape[1]
        # stacked[:, c, :] is (S, genes); barycenter wants (genes, S).
        return np.stack([_wasserstein_barycenter(stacked[:, c, :].T, kernel)
                         for c in range(k)])
    raise ValueError(f"unknown aggregate {aggregate!r}")


def _cross_replicate_agreement(pool):
    """Reproducibility of replicate factorizations, and a medoid reference.

    For every pair of replicates, Hungarian-match their programs by cosine and
    average the matched similarities; the mean over all pairs is a direct
    reproducibility measure (1.0 = every replicate finds identical programs) that
    predicts consensus quality far better than a KMeans silhouette. The medoid
    replicate (highest total agreement) is returned for use as a matching anchor.

    Args:
        pool (numpy.ndarray): ``(n_iter, k, genes)`` replicate spectra.

    Returns:
        tuple: ``(agreement, ref_index, unit_rows, dist_rows)`` where ``unit_rows``
        are L2-normalized programs (for cosine matching) and ``dist_rows`` are the
        L1-normalized programs (distributions, for aggregation).
    """
    from scipy.optimize import linear_sum_assignment

    n_iter = pool.shape[0]
    dists = [_l1_normalize_rows(pool[i]) for i in range(n_iter)]
    units = [d / np.clip(np.linalg.norm(d, axis=1, keepdims=True), 1e-12, None)
             for d in dists]
    pair = np.zeros((n_iter, n_iter))
    for i in range(n_iter):
        for j in range(i + 1, n_iter):
            C = units[i] @ units[j].T
            r, c = linear_sum_assignment(-C)
            pair[i, j] = pair[j, i] = C[r, c].mean()
    ref = int(pair.sum(axis=1).argmax())
    agreement = float(pair[np.triu_indices(n_iter, 1)].mean()) if n_iter > 1 else 1.0
    return agreement, ref, units, dists


def _consensus_by_matching(k, aggregate, ref, units, dists,
                           cost_matrix=None, bary_reg=0.1):
    """Balanced consensus by cross-replicate program matching.

    Hungarian-matches every replicate's ``k`` programs to a medoid reference
    replicate, so each consensus program is the barycenter of exactly one matched
    program per replicate. Unlike free KMeans on the pooled programs, this
    guarantees balance (one program per replicate per cluster). It suits
    high-``k`` cases where KMeans clusters become unbalanced, but -- lacking
    KMeans' local-density outlier filter -- it can be pulled by an outlier
    replicate on noisier (crisp, lower-agreement) data, so it is not the default.
    """
    from scipy.optimize import linear_sum_assignment

    aligned = [dists[ref]]
    for j in range(len(dists)):
        if j == ref:
            continue
        C = units[ref] @ units[j].T
        r, c = linear_sum_assignment(-C)
        order = np.empty(k, dtype=int)
        order[r] = c
        aligned.append(dists[j][order])
    return _l1_normalize_rows(_aggregate_programs(
        np.stack(aligned), aggregate, cost_matrix=cost_matrix, bary_reg=bary_reg))


def _cluster_consensus_spectra(spectra_pool, k, density_threshold,
                               local_neighborhood_size, n_iter=None,
                               method="match", aggregate="mean",
                               cost_matrix=None, bary_reg=0.1):
    """Aggregate pooled replicate spectra into ``k`` consensus programs.

    Adapts cNMF's consensus to spOT-NMF's optimal-transport nature. The key
    change from cNMF is the **aggregation**: spOT-NMF spectra are probability
    distributions over genes, for which cNMF's coordinate-wise *median*
    over-sparsifies and can erase the consensus benefit entirely; a distribution
    barycenter (``"mean"`` = mixture, default, or ``"geomean"`` = KL barycenter)
    is both principled and empirically better (e.g. +0.03 PCC on the mixed-spot
    MOB benchmark, neutral on crisp panels). Clustering strategy (``method``):

    * ``"match"`` (default): balanced cross-replicate Hungarian matching
      (:func:`_consensus_by_matching`) -- guarantees one program per replicate per
      cluster. This is the safer default across ``k``: free KMeans becomes badly
      unbalanced at high ``k`` and can collapse the consensus *below* a single run
      (e.g. k=29 stereo-seq: KMeans 0.16 vs single-run 0.22 vs match 0.23 PCC).
      Requires ``n_iter`` and a pool of exactly ``n_iter*k`` rows.
    * ``"kmeans"``: cNMF-style L2-normalize + local-density outlier filter +
      KMeans, then aggregate each cluster. Competitive only at low ``k`` (its
      density filter can slightly beat matching on crisp low-``k`` panels); avoid
      at high ``k``.

    Regardless of ``method``, ``stability`` is reported as the cross-replicate
    program agreement (:func:`_cross_replicate_agreement`) when the pool can be
    reshaped, falling back to the KMeans silhouette otherwise.

    Args:
        spectra_pool (pandas.DataFrame): ``(n_iter*k) x genes`` pooled replicate
            spectra, ordered replicate-by-replicate (each block of ``k`` rows is
            one replicate's programs).
        k (int): Number of consensus programs.
        density_threshold (float): KMeans-path local-density outlier threshold
            (>=2 disables filtering), matching cNMF's ``local_density_threshold``.
        local_neighborhood_size (float): KMeans-path neighbor fraction (cNMF).
        n_iter (int, optional): Number of replicates (enables agreement +
            ``"match"``).
        method (str, optional): ``"kmeans"`` (default) or ``"match"``.
        aggregate (str, optional): ``"mean"`` (default), ``"geomean"`` (KL
            barycenter) or ``"median"`` (cNMF).

    Returns:
        tuple: ``(consensus_spectra (k x genes DataFrame), stability)``.
    """
    genes = spectra_pool.columns

    reshapeable = n_iter is not None and spectra_pool.shape[0] == n_iter * k
    agreement = float("nan")
    if reshapeable:
        pool = spectra_pool.to_numpy().reshape(n_iter, k, -1)
        agreement, ref, units, dists = _cross_replicate_agreement(pool)

    if method == "match":
        if not reshapeable:
            raise ValueError("method='match' needs n_iter and an n_iter*k pool")
        consensus = _consensus_by_matching(k, aggregate, ref, units, dists,
                                           cost_matrix=cost_matrix, bary_reg=bary_reg)
        return pd.DataFrame(consensus, columns=genes), agreement

    # -- cNMF-style KMeans (default) -----------------------------------------
    from sklearn.cluster import KMeans
    from sklearn.metrics import silhouette_score
    from sklearn.metrics.pairwise import euclidean_distances

    # Rescale every replicate program to unit L2 norm (cNMF's l2_spectra).
    l2 = spectra_pool.div(np.sqrt((spectra_pool ** 2).sum(axis=1)), axis=0)

    # Local-density outlier filtering (skipped when threshold >= 2).
    n_neighbors = int(local_neighborhood_size * l2.shape[0] / k)
    if density_threshold < 2 and n_neighbors >= 1 and l2.shape[0] > n_neighbors + 1:
        dmat = euclidean_distances(l2.values)
        part = np.argpartition(dmat, n_neighbors + 1)[:, : n_neighbors + 1]
        nn = dmat[np.arange(dmat.shape[0])[:, None], part]
        local_density = nn.sum(1) / n_neighbors
        keep = local_density < density_threshold
        if keep.sum() >= k:
            l2 = l2.loc[keep]

    km = KMeans(n_clusters=k, n_init=10, random_state=1)
    labels = pd.Series(km.fit_predict(l2.values), index=l2.index)

    if not reshapeable:
        try:
            agreement = float(silhouette_score(l2.values, labels, metric="euclidean"))
        except Exception:
            pass

    agg = l2.groupby(labels).apply(
        lambda g: pd.Series(_aggregate_rows(
            g.to_numpy(), aggregate, cost_matrix=cost_matrix, bary_reg=bary_reg),
            index=l2.columns))
    agg = agg.div(agg.sum(axis=1), axis=0)
    return agg, agreement


def _aggregate_rows(rows, aggregate, cost_matrix=None, bary_reg=0.1):
    """Aggregate a set of program rows ``(m, genes)`` into one ``(genes,)``."""
    if aggregate == "median":
        return np.median(rows, axis=0)
    if aggregate == "mean":
        return rows.mean(axis=0)
    if aggregate == "geomean":
        return np.exp(np.log(np.clip(rows, 1e-12, None)).mean(axis=0))
    if aggregate == "wasserstein":
        if cost_matrix is None:
            raise ValueError("aggregate='wasserstein' requires cost_matrix")
        dists = _l1_normalize_rows(np.clip(rows, 0, None) + 1e-12)
        kernel = _bary_kernel(cost_matrix, bary_reg)
        return _wasserstein_barycenter(dists.T, kernel)
    raise ValueError(f"unknown aggregate {aggregate!r}")


def run_spotnmf_consensus(
    adata_spatial,
    components,
    n_iter: int = 10,
    seed: int = 42,
    density_threshold: float = 0.5,
    local_neighborhood_size: float = 0.30,
    consensus_method: str = "match",
    aggregate: str = "mean",
    bary_reg: float = 0.1,
    refit: str = "nnls",
    refit_w: float = None,
    return_stability: bool = False,
    **kwargs,
):
    """Consensus spOT-NMF: aggregate replicate OT spectra, then OT-refit usages.

    Structured like cNMF's consensus but adapted to spOT-NMF's optimal-transport
    nature. The key departure from cNMF is the **aggregation**: spOT-NMF spectra
    are probability distributions over genes, so replicate programs are combined
    with a distribution barycenter -- ``aggregate="mean"`` (mixture barycenter,
    default) or ``"geomean"`` (KL barycenter) -- rather than cNMF's
    coordinate-wise ``"median"``, which over-sparsifies distributions and can
    erase the consensus benefit (e.g. +0.03 PCC on the mixed-spot MOB benchmark,
    neutral on crisp panels; see ``scripts/benchmark/``). ``consensus_method``
    selects clustering: ``"match"`` (default, balanced cross-replicate Hungarian
    matching -- the safer choice across ``k``, since free KMeans unbalances at
    high ``k`` and can collapse the consensus below a single run) or ``"kmeans"``
    (cNMF-style with density filtering, competitive only at low ``k``).
    ``stability`` is reported as the cross-replicate program agreement, a truer
    reproducibility measure than a KMeans silhouette.

    The final usage refit uses spOT-NMF's *fixed-spectra optimal-transport* solve
    (:meth:`spotnmf.transform`); note that on genuinely mixed-spot data a plain
    NNLS refit of these consensus spectra can deconvolve better than the OT refit
    (whose usages are more peaked). Concretely:

    1. Run ``n_iter`` independent spOT-NMF factorizations (seeds ``seed`` ..
       ``seed + n_iter - 1``), collecting each replicate's spectra.
    2. Pool and aggregate the replicate spectra into ``components`` consensus
       programs (balanced matching + mixture-barycenter mean by default).
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
        consensus_method (str, optional): Clustering for the consensus spectra --
            ``"match"`` (default, balanced cross-replicate matching) or
            ``"kmeans"`` (cNMF-style with density filtering). Defaults ``"match"``.
        aggregate (str, optional): Distribution barycenter for combining matched
            programs -- ``"mean"`` (default, mixture barycenter), ``"geomean"``
            (KL barycenter; best on low-k mixed data e.g. MOB), ``"wasserstein"``
            (entropic **Wasserstein barycenter** under the gene ground cost -- the
            fully OT-native aggregation; best on crisp/higher-k panels and lets the
            safe ``match`` clustering reach ``kmeans``-level accuracy) or
            ``"median"`` (cNMF). Defaults ``"mean"``.
        bary_reg (float, optional): Entropic regularization for the ``"wasserstein"``
            barycenter (larger = more OT smoothing across similar genes). ~0.1
            worked best in the benchmark. Defaults 0.1.
        refit (str, optional): How to solve the final usages against the frozen
            consensus spectra -- ``"nnls"`` (default; scikit-learn least-squares
            refit as in cNMF) or ``"ot"`` (fixed-spectra entropic-OT solve via
            :meth:`transform`). NNLS gives graded usages that deconvolve mixed
            spots better and the best PCC across the benchmark datasets; OT keeps
            the peaked OT usage geometry and wins RMSE/JSD on crisp panels.
            Defaults ``"nnls"``.
        refit_w (float, optional): Usage-entropy for the ``"ot"`` refit only,
            decoupled from the replicate training ``w`` (smaller = more peaked).
            ``None`` reuses the training ``w``. Defaults None.
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
        "dtype": "float64", "analytic_grad": None,
        # Applies to the refit="ot" usage readout only; NNLS usages have no
        # softmin temperature, so readout_temp is a no-op when refit="nnls".
        "readout_temp": 1.0,
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

    # 2. Aggregate into consensus programs. For the OT-native Wasserstein
    #    barycenter, build the gene ground cost (same cosine-on-normalized-
    #    reference geometry the model uses) so mass transports between similar
    #    genes; other aggregates ignore it.
    cost_matrix = None
    if aggregate == "wasserstein":
        cost_matrix = _gene_ground_cost(
            adata_spatial, gene_list, defaults["cost"], defaults["normalize_rows"])
    median_spectra, stability = _cluster_consensus_spectra(
        spectra_pool, k, density_threshold, local_neighborhood_size,
        n_iter=n_iter, method=consensus_method, aggregate=aggregate,
        cost_matrix=cost_matrix, bary_reg=bary_reg,
    )

    # Consensus spectra as (genes x k) and (k x genes), aligned to gene_list.
    consensus_kg = median_spectra.reindex(columns=gene_list).to_numpy()  # (k x genes)
    consensus_gk = consensus_kg.T                                        # (genes x k)

    # 3. Refit the usages against the FROZEN consensus spectra. Two options:
    #   * "nnls" (default): least-squares usage refit (scikit-learn
    #     ``non_negative_factorization(..., update_H=False)``), i.e. cNMF's refit.
    #     Its graded usages deconvolve genuinely mixed spots better and give the
    #     best PCC across datasets in the spOT-NMF benchmark (STARmap/seqFISH/MOB).
    #   * "ot": the fixed-spectra entropic-OT solve (:meth:`transform`). Keeps the
    #     OT usage geometry (more peaked); best RMSE/JSD on crisp panels. The
    #     usage peakiness can be tuned separately from training via ``refit_w``.
    if refit not in ("nnls", "ot"):
        raise ValueError(f"refit must be 'nnls' or 'ot', got {refit!r}")
    seed_all(seed)
    refit_losses = []
    if refit == "nnls":
        from sklearn.decomposition import non_negative_factorization

        Xhv = adata_spatial[:, adata_spatial.var["highly_variable"]].X
        Xhv = Xhv.toarray() if hasattr(Xhv, "toarray") else np.asarray(Xhv)
        Xhv = np.ascontiguousarray(Xhv, dtype=np.float64)
        # scikit-learn's NNLS requires H.dtype == X.dtype. The consensus spectra
        # inherit the replicate dtype (e.g. float32 when the replicates ran in
        # float32), so cast them to X's float64 before the refit.
        W, H_used, _ = non_negative_factorization(
            Xhv, H=np.ascontiguousarray(consensus_kg, dtype=Xhv.dtype),
            n_components=k, init="custom", update_H=False, solver="cd",
            max_iter=500, random_state=seed)
        W_out = W                       # (spots x k)
        H_out = H_used.T                # (genes x k)
    else:  # refit == "ot"
        model = spotnmf(
            factors=k,
            h_regularization=defaults["h"],
            w_regularization=refit_w if refit_w is not None else defaults["w"],
            eps=defaults["eps"],
            cost=defaults["cost"],
            pca_cost=False,
        )
        refit_adata = adata_spatial.copy()
        model.transform(
            refit_adata,
            spectra=consensus_kg,  # (k x genes)
            lr=defaults["lr"],
            optim_name=defaults["optim_name"],
            tol_inner=defaults["tol_inner"],
            tol_outer=defaults["tol_outer"],
            normalize_rows=defaults["normalize_rows"],
            max_iter=defaults["max_iter"],
            max_iter_inner=defaults["max_iter_inner"],
            impute=False,
            analytic_grad=defaults["analytic_grad"],
            readout_temp=defaults["readout_temp"],
            dtype=_resolve_dtype(defaults["dtype"]),
            device=defaults["device"],
        )
        W_out = refit_adata.obsm["W_OT"]      # (spots x k)
        H_out = refit_adata.uns["H_OT"]       # (genes x k)
        refit_losses = model.losses

    # Write the consensus factors back onto the input AnnData so callers get the
    # same obsm["W_OT"]/uns["H_OT"] contract as run_spotnmf (the replicate loop
    # above only ever mutated copies, so the original is untouched until now).
    if adata_spatial.uns is None:
        adata_spatial.uns = {}
    adata_spatial.uns["H_OT"] = H_out
    adata_spatial.obsm["W_OT"] = W_out

    df_topics_per_spot = pd.DataFrame(W_out, index=adata_spatial.obs.index)
    df_genes_per_topic = pd.DataFrame(H_out, index=gene_list)
    results = {
        "topics_per_spot": df_topics_per_spot,
        "genes_per_topic": df_genes_per_topic,
    }
    for key_matrix in results:
        results[key_matrix].columns = [f"ot_{x + 1}" for x in results[key_matrix].columns]

    if return_stability:
        return results, refit_losses, stability
    return results, refit_losses
