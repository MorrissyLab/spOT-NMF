from typing import Iterable, List
import torch
from scipy.spatial.distance import cdist
import numpy as np

def clean_mixed_gene_names(genes_list, genome):
    """Clean gene names by stripping genome-specific prefixes.

    Args:
        genes_list (list): List of gene name strings that may contain genome
            prefixes (``"mm10___"`` or ``"GRCh38_"``).
        genome (str): The genome type. If ``"None"``, both prefixes are
            stripped; if ``"mm10"``, only the ``"mm10___"`` prefix is
            stripped; otherwise only the ``"GRCh38_"`` prefix is stripped.

    Returns:
        list: The cleaned list of gene name strings.
    """
    if genome == 'None':
        return [x.replace('mm10___', '').replace('GRCh38_', '') for x in genes_list]
    elif genome == 'mm10':
        return [x.replace('mm10___', '') for x in genes_list]
    else:
        return [x.replace('GRCh38_', '') for x in genes_list]
    
def reference_dataset(
    X, dtype: torch.dtype, device: torch.device
) -> torch.Tensor:
    """Transpose the dataset, densify if sparse, and convert to a Tensor.

    Args:
        X (numpy.ndarray): The input data, with observations as rows and
            variables as columns.
        dtype (torch.dtype): The dtype of the output tensor.
        device (torch.device): The device to create the tensor on.

    Returns:
        torch.Tensor: The reference dataset ``A`` (variables by observations).

    Raises:
        AssertionError: If the input data contains negative values.
    """

    # Keep only the highly variable features.
    A = X.T

    # Check that the dataset is positive.
    assert (A >= 0).all()

    # If the dataset is sparse, make it dense.
    try:
        A = A.todense()
    except:
        pass

    # Send the matrix `A` to PyTorch.
    return torch.from_numpy(A).to(device=device, dtype=dtype).contiguous()


def compute_ground_cost(
    features,
    cost: str,
    eps: float,
    force_recompute: bool,
    cost_path: str,
    dtype: torch.dtype,
    device: torch.device,
) -> torch.Tensor:
    """Compute the ground cost kernel (not lazily).

    Computes a pairwise distance matrix between features, normalizes it by
    ``eps`` times its maximum, and exponentiates it to form the kernel.

    Args:
        features (numpy.ndarray): An array with the features to compute the
            cost on.
        cost (str): The metric used to compute the cost. All SciPy ``cdist``
            metrics are allowed; ``"ones"`` uses a uniform off-diagonal cost.
        eps (float): The entropic regularization used to scale the cost.
        force_recompute (bool): Whether to recompute the cost even if a cost
            matrix is already saved at ``cost_path``.
        cost_path (str): The path to look for or save the cost matrix as a
            ``.npy`` file.
        dtype (torch.dtype): The dtype of the output tensor.
        device (torch.device): The device of the output tensor.

    Returns:
        torch.Tensor: The ground cost kernel.
    """

    # Initialize the `recomputed variable`.
    recomputed = False

    # If we force recomputing, then compute the ground cost.
    if force_recompute:
        K = cdist(features, features, metric=cost)
        recomputed = True

    # If the cost is not yet computed, try to load it or compute it.
    if not recomputed:
        try:
            K = np.load(cost_path)
        except:
            if cost == "ones":
                K = 1 - np.eye(features.shape[0])
            else:
                K = cdist(features, features, metric=cost)
            recomputed = True

    # If we did recompute the cost, save it.
    if recomputed and cost_path:
        np.save(cost_path, K)

    K = torch.from_numpy(K).to(device=device, dtype=dtype)
    K /= eps * K.max()

    # Compute the kernel K.
    K = torch.exp(-K).to(device=device, dtype=dtype)

    return K


def normalize_tensor(X: torch.Tensor) -> torch.Tensor:
    """Normalize a tensor so that each column sums to one.

    Args:
        X (torch.Tensor): The tensor to normalize.

    Returns:
        torch.Tensor: The column-normalized tensor.
    """
    return X / X.sum(0)


def entropy(
    X: torch.Tensor, min_one: bool = False, rescale: bool = False
) -> torch.Tensor:
    r"""
    Computes the entropy of a tensor, defined as :math:`E(X) = \langle X, \log X - 1 \rangle`.

    Args:
        X (torch.Tensor): 
            The input tensor for which the entropy is computed.
        min_one (bool, optional): 
            Whether to include the :math:`-1` term in the formula. Defaults to False.
        rescale (bool, optional): 
            If True, rescales the result so that the entropy lies between 0 and 1 (when min_one=False). 
            Defaults to False.

    Returns:
        torch.Tensor: 
            The computed entropy of the input tensor.
    """
    offset = 1 if min_one else 0
    scale = X.shape[1] * np.log(X.shape[0]) if rescale else 1
    return -torch.sum(X * (torch.nan_to_num(X.log()) - offset)) / scale


def entropy_dual_loss(Y: torch.Tensor) -> torch.Tensor:
    """Compute the Legendre dual of the entropy.

    Args:
        Y (torch.Tensor): The input tensor.

    Returns:
        torch.Tensor: The entropy dual loss.
    """
    return -torch.logsumexp(Y, dim=0).sum()


def ot_dual_loss(
    A: torch.Tensor, G: torch.Tensor, K: torch.Tensor, eps: float, dim=(0, 1)
) -> torch.Tensor:
    """Compute the Legendre dual of the entropic OT loss.

    Args:
        A (torch.Tensor): The reference data tensor.
        G (torch.Tensor): The dual variable tensor.
        K (torch.Tensor): The ground cost kernel tensor.
        eps (float): The entropic regularization.
        dim (tuple, optional): The dimensions over which to sum the loss.
            Defaults to (0, 1).

    Returns:
        torch.Tensor: The OT dual loss.
    """

    log_fG = G / eps

    # Compute the non stabilized product.
    scale = log_fG.max(0).values
    prod = torch.log(K @ torch.exp(log_fG - scale)) + scale

    # Compute the dot product with A.
    return eps * torch.sum(A * prod, dim=dim)


def early_stop(history: List, tol: float, nonincreasing: bool = False) -> bool:
    """Decide whether to stop early based on a loss history and tolerance.

    Args:
        history (list): The loss history.
        tol (float): The tolerance before early stopping.
        nonincreasing (bool, optional): When True, also stops early if the
            loss increases beyond the tolerance over the last iterations.
            Defaults to False.

    Returns:
        bool: Whether to stop early.

    Raises:
        ValueError: If the most recent loss value is not finite.
    """
    # If we have a nan or infinite, die.
    if len(history) > 0 and not torch.isfinite(history[-1]):
        raise ValueError("Error: Loss is not finite!")

    # If the history is too short, continue.
    if len(history) < 3:
        return False

    # If the next value is worse, stop (not normal!).
    # if history[-1] < 0.005 and history[-1]>=0:
    #     return True

    if nonincreasing and (history[-1] - history[-3]) > tol:
        return True

    # If the next value is close enough, stop.
    if abs(history[-1] - history[-2]) < tol:
        return True

    # Otherwise, keep on going.
    return False