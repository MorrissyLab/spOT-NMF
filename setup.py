from setuptools import setup, find_packages
import os

# Read requirements.txt
requirements = []
if os.path.exists("requirements.txt"):
    with open("requirements.txt", "r") as f:
        requirements = [line.strip() for line in f if line.strip() and not line.startswith("#")]

# Read README for long description
long_description = ""
if os.path.exists("README.md"):
    with open("README.md", "r", encoding="utf-8") as f:
        long_description = f.read()

setup(
    name="spotnmf",
    version="0.1.0",
    description="Optimal Transport-based Matrix Factorization for spatial transcriptomics deconvolution.",
    long_description=long_description,
    long_description_content_type="text/markdown",
    author="Aly Abdelkareem",
    url="https://github.com/MorrissyLab/spOT-NMF",
    packages=find_packages(),
    install_requires=requirements,
    python_requires=">=3.8",
    entry_points={
        'console_scripts': [
            'spotnmf=cli:main',  # Adjust if needed to match your CLI location
        ],
    },
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
        "Operating System :: OS Independent",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
    ],
    include_package_data=True,
    zip_safe=False,
    license="GPL-3.0",
)
