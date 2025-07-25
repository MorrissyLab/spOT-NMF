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
from . import utils

def seed_all(seed):
    np.random.seed(seed)  # Seed numpy (covers scipy)
    torch.manual_seed(seed)  # Seed pytorch (CPU)
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
            small values mean sparse vectors. Defaults to 1e-3.
        eps (float, optional):
            The entropy parameter for epsilon transport. Large values
            decrease importance of individual genes. Defaults to 5e-2.
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
        w_regularization: float = 1e-3,
        eps: float = 5e-2,
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
        """Initialize parameters based on input data.

        Args:
            adata_spatial:
                The input adata_spatial object.
            dtype (torch.dtype):
                The dtype to work with.
            device (torch.device):
                The device to work on.
            force_recompute (bool, optional):
                Whether to recompute the ground cost. Defaults to False.
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
        lr: float = 1,
        optim_name: str = "lbfgs",
        tol_inner: float = 1e-12,
        tol_outer: float = 1e-4,
        normalize_rows: bool = False,
        impute: bool = False,
        batch_size = 512,
    ) -> None:
        """Train the spotnmf model on an input adata object.

        Args:
            adata_spatial :
                The input adata object.
            max_iter_inner (int, optional):
                How many iterations for the inner optimization loop
                (optimizing H, or W). Defaults to 1_000.
            max_iter (int, optional):
                How many interations for the outer optimization loop (how
                many successive optimizations of H and W). Defaults to 100.
            device (torch.device, optional):
                The device to work on. Defaults to 'cpu'.
            dtype (torch.dtype, optional):
                The dtype to work with. Defaults to torch.double.
            lr (float, optional):
                The learning rate for the optimizer. The default is set
                for LBFGS and should be changed otherwise. Defaults to 1.
            optim_name (str, optional):
                The optimizer to use (`lbfgs`, `sgd` or `adam`). LBFGS
                is advised, but requires more memory. Defaults to "lbfgs".
            tol_inner (float, optional):
                The tolerance for the inner iterations before early stopping.
                Defaults to 1e-12.
            tol_outer (float, optional):
                The tolerance for the outer iterations before early stopping.
                Defaults to 1e-4.
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

    def build_optimizer(
        self, params, lr: float, optim_name: str
    ) -> torch.optim.Optimizer:
        """Generates the optimizer. The PyTorch LBGS implementation is
        parametrized following the discussion in https://discuss.pytorch.org/
        t/unclear-purpose-of-max-iter-kwarg-in-the-lbfgs-optimizer/65695.

        Args:
            params (Iterable of Tensors):
                The parameters to be optimized.
            lr (float):
                Learning rate of the optimizer.
            optim_name (str):
                Name of the optimizer, among `'lbfgs'`, `'sgd'`, `'adam'`

        Returns:
            torch.optim.Optimizer: The optimizer.
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
        """Optimize a given function.

        Args:
            loss_fn (Callable): The function to optimize.
            max_iter (int): The maximum number of iterations.
            history (List): A list to append the losses to.
            tol (float): The tolerance before early stopping.
            pbar (A tqdm progress bar): The progress bar.
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
                gpu_mem_alloc = torch.cuda.memory_allocated(device=device)

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
        """Compute the total dual loss. This is only used by the user and for,
        early stopping, not by the optimization algorithm.

        Returns:
            torch.Tensor: The loss
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
        """Computes the loss for the update of `H`.

        Returns:
            torch.Tensor: The loss.
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
        """Return the loss for the optimization of W

        Returns:
            torch.Tensor: The loss
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
    """
    Run the spotnmf model with flexible parameters provided via kwargs.

    Args:
        adata_spatial: Input spatial data in an AnnData object.
        components: Number of factorss.
        mod: The modality to use (default is 'simulated').
        **kwargs: Additional parameters for OT and training, including:
            - h: Regularization parameter for h (default 0.01).
            - w: Regularization parameter for w (default 1e-2).
            - eps: Entropic regularization (default 5e-3).
            - lr: Learning rate (default 0.001).
            - optim_name: Optimizer name (default 'adam').
            - cost: Cost function (default 'cosine').
    """
    seed_all(seed)
    
    # Define default values
    defaults = {
        "cost": "cosine",
        "optim_name": "adam",
        "tol_inner": 1e-12,
        "tol_outer": 0.00001,
        "max_iter": 100,
        "max_iter_inner": 1000
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
        device='cuda',
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
