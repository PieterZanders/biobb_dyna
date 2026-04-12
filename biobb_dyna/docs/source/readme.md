[![](https://img.shields.io/github/v/tag/bioexcel/biobb_dyna?label=Version)](https://GitHub.com/bioexcel/biobb_dyna/tags/)
[![](https://img.shields.io/pypi/v/biobb-dyna.svg?label=Pypi)](https://pypi.python.org/pypi/biobb-dyna/)
[![](https://img.shields.io/conda/vn/bioconda/biobb_dyna?label=Conda)](https://anaconda.org/bioconda/biobb_dyna)
[![](https://img.shields.io/badge/License-Apache%202.0-blue.svg)](https://opensource.org/licenses/Apache-2.0)

[![](https://readthedocs.org/projects/biobb-dyna/badge/?version=latest&label=Docs)](https://biobb-dyna.readthedocs.io/en/latest/?badge=latest)
[![](https://img.shields.io/website?down_message=Offline&label=Biobb%20Website&up_message=Online&url=https%3A%2F%2Fmmb.irbbarcelona.org%2Fbiobb%2F)](https://mmb.irbbarcelona.org/biobb/)

# biobb_dyna

### Introduction

**biobb_dyna** is a BioExcel Biobb package for **molecular dynamics correlation analysis** and **protein network graphs**. It wraps [MDiGest](https://pypi.org/project/mdigest/) (dynamical correlation / covariance matrices), [MDAnalysis](https://www.mdanalysis.org/) for trajectories and selections, and [NetworkX](https://networkx.org/) for community detection and path analysis.

Biobb packages add a consistent Python and command-line interface on top of underlying tools so they can be composed in workflows (CWL, Galaxy, etc.).

The latest docs: [API modules](https://biobb-dyna.readthedocs.io/en/latest/modules.html) · [Command line](https://biobb-dyna.readthedocs.io/en/latest/command_line.html)

### Version

v5.0.0

### Installation

Using pip (dependencies such as `mdigest`, `MDAnalysis`, and `networkx` should be available in your environment):

```bash
pip install "biobb_dyna>=5.0.0"
```

Using conda:

```bash
conda install -c bioconda "biobb_dyna>=5.0.0"
```

A full developer environment (including test tools) can be created from `conda_env/environment.yml` at the repository root.

### Docker / Singularity

See the [GitHub repository](https://github.com/bioexcel/biobb_dyna) README for container images and example commands.

### Citation

If you use BioExcel Biobb in your research, please cite the Biobb reference publications listed on the [BioExcel website](https://bioexcel.eu/).
