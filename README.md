# Alchemd

[![License](https://img.shields.io/badge/license-CC%20BY--NC--SA--4.0-blue.svg)](#license)
[![Platform](https://img.shields.io/badge/platform-Linux-lightgrey.svg)](#system-requirements)
[![Python](https://img.shields.io/badge/python-3.9+-blue.svg)](#installation)
[![CUDA](https://img.shields.io/badge/CUDA-12.0-green.svg)](#system-requirements)

> A computational toolkit for alchemical free energy perturbation (FEP) calculations in computer-aided drug design

## Overview

Alchemd is a Python-based toolkit for Relative Binding Free Energy (RBFE) calculations. The software provides an integrated platform for molecular preparation, dynamics simulation, and thermodynamic analysis with enhanced convergence algorithms.

## System Requirements

- **Operating System**: Linux (Ubuntu 18.04+) or Windows with WSL2
- **Hardware**: NVIDIA GPU with CUDA support (RTX 4060 Ti+ recommended)
- **Dependencies**: CUDA 12.0, Conda/Mamba, GCC 9+
- **License**: Amber license required ([apply here](https://ambermd.org/))
- **Optional dependency**: Gaussian (for NSAA non-standard amino acid module)

## Installation

### Automated Installation

**Prerequisites**: Before running the automated installation scripts, ensure you have:
- **Conda or Mamba** installed and properly configured
- **CUDA 12.0** installed and accessible in your system PATH

Execute the following installation scripts sequentially:

```bash
# 1. System environment setup
chmod +x install/setup_linux.sh
./install/setup_linux.sh

# 2. Python environment configuration
chmod +x install/setup_conda.sh
./install/setup_conda.sh

# 3. PMEMD compilation (requires Amber license)
chmod +x install/build_amber.sh
./install/build_amber.sh /path/to/pmemd24
```

### Installation Guide

For detailed installation procedures, refer to:
- [`docs/Installation_Guide_EN.md`](docs/Installation_Guide_EN.md)

## Quick Start

```bash
# Structure preparation
python core/prepare_file.py -p protein.pdb -l ligands.sdf -o workspace --cs-align -d -auto

# Molecular dynamics simulation  
python core/run_md.py -c -ln 50 -p ./workspace/prepare -m Amber -auto

# Results analysis
python core/analyze_result.py -auto -c -w ./workspace/run -m Amber

# GUI interface
python gui/gui.py
```

## Documentation

- **User Guide**: [`docs/Alchemd_User_Guide_EN.md`](docs/Alchemd_User_Guide_EN.md)
- **Installation Guide**: [`docs/Installation_Guide_EN.md`](docs/Installation_Guide_EN.md)
- **Example Data**: [`example/`](example/) directory

## Citation

If you use Alchemd in your research, please cite:

```bibtex
@software{alchemd2025,
  title={Alchemd: A toolkit for alchemical free energy perturbation calculations},
  year={2025},
  url={https://github.com/ZheLi-Lab/Alchemd-toolkit}
}
```

**Methodological References:**
- CAR method: [DOI](https://pubs.acs.org/doi/10.1021/acs.jctc.4c00939)  
- cs-fep implementation: [DOI](https://doi.org/10.1016/j.apsb.2024.06.021)

## License

Alchemd is dual-licensed:

- **Non-commercial use**: Licensed under [CC BY-NC-SA 4.0](LICENSE.md) - free for personal, academic, and non-profit use
- **Commercial use**: Requires separate commercial license - contact us for licensing inquiries

See the [LICENSE](LICENSE.md) file for complete terms and conditions.

## Support

For questions and issues, please use [GitHub Issues](../../issues).