from setuptools import setup
import os
from os import path
import io

this_directory = path.abspath(path.dirname(__file__))
with open(path.join(this_directory, "README.md"), encoding="utf-8") as f:
    long_description = f.read()
    
import re

VERSIONFILE = "snp2fasta/_version.py"
verstrline = open(VERSIONFILE, "rt").read()
VSRE = r"^__version__ = ['\"]([^'\"]*)['\"]"
mo = re.search(VSRE, verstrline, re.M)
if mo:
    verstr = mo.group(1)
else:
    raise RuntimeError("Unable to find version string in %s." % (VERSIONFILE,))


def _read(*parts, **kwargs):
    filepath = os.path.join(os.path.dirname(__file__), *parts)
    encoding = kwargs.pop("encoding", "utf-8")
    with io.open(filepath, encoding=encoding) as fh:
        text = fh.read()
    return text


def get_requirements(path):
    content = _read(path)
    return [
        req
        for req in content.split("\n")
        if req != "" and not (req.startswith("#") or req.startswith("-"))
    ]


setup_requires = []

INSTALL_REQUIRES = get_requirements("requirements.txt")

setup(
    name="snp2fasta",
    version=verstr,
    packages=["snp2fasta"],
    entry_points={
        "console_scripts": [
            "snp2fasta = snp2fasta.snp2fasta:main",
        ]
    },
    setup_requires=setup_requires,
    install_requires=INSTALL_REQUIRES,
    python_requires=">=3.8",
    description="SNP table to fasta of ref/alt with flanks",
    long_description=long_description,
    long_description_content_type="text/markdown",
    project_urls={
        "Source": "https://github.com/efriman/snp2fasta",
        "Issues": "https://github.com/efriman/snp2fasta/issues",
    },
    author="Elias Friman",
    author_email="elias.friman@gmail.com",
    classifiers=[
        "Programming Language :: Python",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    zip_safe=False,
)
