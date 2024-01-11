# hiPhive module

This project introduces a set of codes for computing and analyzing high order force constants

The core technology behind this framework is based on a graph-like representation of data, based on hiPhive implementation.

This technology is suitable for different applications: from vibrational stability to thermal conductivity.

## Features

- Generating distorted structures from some reference POSCAR.
- Generating correponding force constants after energy calculations (following phonopy distribution of folders and files).
- Computing phonon lifetimes.
- Comparing with harmonic analysis.
- Analyzing convergence of clusters.
- Analyzing convergence with data points.

Please be aware that the code is under active development, bug reports are welcomed in the GitHub issues!

## Installation

To download the repository and install the dependencies:

```bash
git clone https://github.com/CibranLopez/hiPhive.git
cd hiPhive
pip3 install -r requirements.txt
```

## Execution

An user-friendly jupyter notebook has been developed, which can be run locally with common dependencies. It generates all necesserary folders and explanative files.

## Authors

This project is being developed by:

 - Cibrán López Álvarez

## Contact, questions and contributing

If you have questions, please don't hesitate to reach out at: cibran.lopez@upc.edu
