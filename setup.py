import re

from pathlib import Path
from setuptools import setup, find_packages

with open("README.md") as readme:
    long_description = readme.read()


def get_version():
    """Get version number from __init__.py"""
    version_file = Path(__file__).resolve().parent / "clinker" / "__init__.py"
    version_match = re.search(
        r"^__version__ = ['\"]([^'\"]*)['\"]", version_file.read_text(), re.M
    )
    if version_match:
        return version_match.group(1)
    raise RuntimeError("Failed to find version string")


setup(
    name="clinker",
    author="Cameron Gilchrist",
    version=get_version(),
    description="",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/gamcil/clinker",
    packages=find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    install_requires=[
        "biopython>=1.78",
        "numpy>=1.13.3",
        "scipy>=1.3.3",
        "disjoint-set>=0.7.1",
        "gffutils",
    ],
    python_requires=">=3.6",
    entry_points={"console_scripts": ["clinker=clinker.main:main"]},
    include_package_data=True,
)
