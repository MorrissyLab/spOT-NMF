# Installation

`spOT-NMF` requires **Python ≥ 3.12**. PyTorch is installed separately so you can
pick the build (CPU or CUDA) that matches your platform.

## Recommended: uv

We recommend [**uv**](https://docs.astral.sh/uv/) for a fast, reproducible setup.

```bash
# 1. Create and activate an isolated environment (uv fetches Python 3.12 if needed)
uv venv --python 3.12
# Linux/macOS:  source .venv/bin/activate
# Windows:      .venv\Scripts\activate

# 2. Install PyTorch for your platform (see pytorch.org)
#    CPU-only:
uv pip install torch --index-url https://download.pytorch.org/whl/cpu
#    CUDA 12.x (Linux/Windows with NVIDIA GPUs):
#    uv pip install torch torchvision torchaudio --index-url https://download.pytorch.org/whl/cu121

# 3. Install spOT-NMF
uv pip install spot-nmf
```

## Alternative: pip

```bash
python -m venv .venv && source .venv/bin/activate   # Windows: .venv\Scripts\activate
pip install torch --index-url https://download.pytorch.org/whl/cpu
pip install spot-nmf
```

## GPU (CUDA) builds

Install the CUDA build of PyTorch that matches your driver **before** installing
`spOT-NMF`. See the [PyTorch install matrix](https://pytorch.org/get-started/locally/)
for the correct `--index-url` for your CUDA version, for example:

```bash
pip install torch torchvision torchaudio --index-url https://download.pytorch.org/whl/cu121
```

## Development install

To work on the package from a clone:

```bash
git clone https://github.com/MorrissyLab/spOT-NMF.git
cd spOT-NMF
uv pip install -e ".[dev]"      # or: pip install -e ".[dev]"
```

To build the documentation locally:

```bash
uv pip install -e ".[docs]"
sphinx-build -b html docs/source docs/_build/html
```

## Verifying the install

```python
import spotnmf
print(spotnmf.__version__)
```
