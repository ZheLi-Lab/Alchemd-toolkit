# Alchemd Installation Guide

---
**License**: This documentation is part of Alchemd, licensed under CC BY-NC-SA 4.0 for academic use. Commercial licensing available separately. See [LICENSE.md](../LICENSE.md) for complete terms.
---

[![License](https://img.shields.io/badge/license-CC%20BY--NC--SA--4.0-blue.svg)](../LICENSE.md)
[![Platform](https://img.shields.io/badge/platform-Linux-lightgrey.svg)](#system-requirements)
[![Python](https://img.shields.io/badge/python-3.9+-blue.svg)](#python-environment)
[![CUDA](https://img.shields.io/badge/CUDA-12.0-green.svg)](#cuda-installation)

> Comprehensive installation guide for Alchemd, including automated and manual installation methods.

## Overview

This guide provides complete installation methods for Alchemd. Automated scripts are recommended for quick deployment, while manual installation is available for customized configuration based on specific requirements.

**Note**: This guide has been updated for the simplified PMEMD single-component architecture. AmberTools is no longer required.

---

## Table of Contents

- [Chapter 0: Automated Installation](#chapter-0-automated-installation)
- [Chapter 1: System Environment Setup](#chapter-1-system-environment-setup)
- [Chapter 2: Python Environment Configuration](#chapter-2-python-environment-configuration)
- [Chapter 3: External Tools Configuration](#chapter-3-external-tools-configuration)
- [Chapter 4: PMEMD Compilation](#chapter-4-pmemd-compilation)
- [Chapter 5: Alchemd Configuration](#chapter-5-alchemd-configuration)
- [Appendix: Troubleshooting](#appendix-troubleshooting)

---

## Chapter 0: Automated Installation

## Overview

Alchemd provides automated installation scripts to streamline the complex environment setup process. The automated installation sequentially executes system environment preparation, Python environment configuration, and PMEMD compilation through three scripts, suitable for rapid deployment in standard Linux environments.

**Note**: The automated scripts work for most standard configurations. For special environments or dependency conflicts, refer to the manual installation guide for customized configuration.

---

## 0.1 Prerequisites

### 0.1.1 System Requirements

- **Operating System**: Ubuntu 18.04+, CentOS 7+, Fedora, openSUSE, or Arch Linux
- **CUDA Toolkit**: CUDA 12.0 (required for GPU acceleration)
- **Amber License**: Apply for license at https://ambermd.org/

### 0.1.2 Required Software

**Conda/Mamba Package Manager**:
```bash
# Recommended: Install Miniconda
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh

# Or install Mamba (better performance)
conda install mamba -n base -c conda-forge
```

**NVIDIA Drivers and CUDA**:
- Ensure compatible NVIDIA drivers are installed
- CUDA 12.0 toolkit must be accessible via `nvcc --version`

---

## 0.2 Installation Steps

### 0.2.1 System Environment Setup

Execute the first script to install system dependencies and GCC toolchain:

```bash
chmod +x install/setup_linux.sh
./install/setup_linux.sh
```

This script automatically:
- Detects Linux distribution and package manager
- Installs base development packages (build-essential, cmake, flex, bison, etc.)
- Configures GCC 9+ compiler environment
- Verifies compiler versions and availability

### 0.2.2 Python Environment Configuration

Execute the second script to create Alchemd-specific Python environment:

```bash
chmod +x install/setup_conda.sh
./install/setup_conda.sh
```

This script automatically:
- Creates conda environment named `alchemd` (Python 3.9)
- Installs scientific computing dependencies (numpy, scipy, pandas, matplotlib, etc.)
- Installs molecular simulation tools (rdkit, openmm, parmed, lomap2, etc.)
- Downloads and configures external tools (WatVina, etc.)
- Initializes Alchemd configuration files

### 0.2.3 PMEMD Compilation

Execute the third script to compile PMEMD with CUDA support:

```bash
chmod +x install/build_amber.sh
./install/build_amber.sh /path/to/pmemd24
```

This script automatically:
- Applies Alchemd-specific patches to PMEMD source code
- Configures CUDA compilation options
- Executes parallel compilation and installation
- Sets up environment variable integration

**Parameters**:
- `/path/to/pmemd24`: Path to extracted PMEMD source directory
- Optional second parameter specifies compilation parallelism (default: 2 cores)

---

## 0.3 Installation Verification

### 0.3.1 Activate Environment

```bash
conda activate alchemd
```

### 0.3.2 Verify Python Environment

```bash
# Test core modules
python core/prepare_file.py --help

# Verify key dependencies
python -c "import numpy, rdkit, openmm, parmed; print('Core dependencies OK')"
```

### 0.3.3 Verify PMEMD Installation

```bash
# Check PMEMD CUDA version
which pmemd.cuda_SPFP
pmemd.cuda_SPFP --help

# Run environment test script
./install/test_amber_env.sh
```

### 0.3.4 Verify Configuration Files

```bash
# Check configuration file generation
ls -la core/configs.toml

# Verify path configuration
python -c "from core.analyze_tools.SettingManager import SettingManager; s=SettingManager(); print('Config loaded successfully')"
```

---

## 0.4 Troubleshooting

### 0.4.1 Common Issues

**GCC Version Too Old**:
Scripts will automatically attempt to install GCC 11. If failed, install manually or refer to the manual installation guide.

**CUDA Environment Issues**:
Ensure `nvcc --version` executes properly and displays CUDA 12.0 version.

**Package Conflicts**:
Try cleaning conda cache: `conda clean --all`, then re-run scripts.

**Network Download Failures**:
Configure mirror sources or use proxy network environment.

### 0.4.2 Reinstallation

To re-execute specific steps:

```bash
# Reconfigure Python environment (preserve existing environment)
./install/setup_conda.sh

# Update dependencies only (skip environment creation)
./install/setup_conda.sh --skip-env

# Recompile PMEMD
./install/build_amber.sh /path/to/pmemd24
```

---

## Chapter 1: System Environment Setup

### 1.1 Ubuntu System Update

Update system packages and install essential tools:

```bash
# Update package lists
sudo apt update

# Upgrade installed packages (optional but recommended)
sudo apt upgrade -y

# Install essential system tools
sudo apt install -y curl wget git vim
```

**Verification**:
```bash
# Check system version
lsb_release -a
# Should display Ubuntu 18.04+ version
```

### 1.2 Development Tools Installation

Install compilation and development tools:

```bash
# Install build essentials
sudo apt install -y build-essential

# Install compilers and toolchain
sudo apt install -y gcc g++ gfortran

# Install build systems
sudo apt install -y make cmake

# Install additional required tools
sudo apt install -y tcsh flex bison patch bc unzip

# Install development libraries
sudo apt install -y xorg-dev libz-dev libbz2-dev
```

**Verification**:
```bash
# Check if tools are installed successfully
gcc --version      # Should display GCC version
g++ --version      # Should display G++ version
gfortran --version # Should display Gfortran version
make --version     # Should display Make version
cmake --version    # Should display CMake version
```

### 1.3 GCC Compiler Installation (GCC 9+)

Check and ensure GCC version meets requirements:

```bash
# Check current GCC version
gcc --version

# Extract major version number
gcc -dumpversion | cut -d. -f1
```

**If version >= 9**:
```bash
echo "GCC version meets requirements, proceed to next step"
```

**If version < 9**:
```bash
# Install newer GCC version
sudo apt install -y gcc-11 g++-11 gfortran-11

# Set default version (optional)
sudo update-alternatives --install /usr/bin/gcc gcc /usr/bin/gcc-11 60 \
    --slave /usr/bin/g++ g++ /usr/bin/g++-11 \
    --slave /usr/bin/gfortran gfortran /usr/bin/gfortran-11

# Verify version
gcc --version
```

**Dependencies**:
- **GCC 9+**: Minimum requirement for PMEMD compilation
- **GCC 11+**: Recommended version with better optimization and compatibility
- **Gfortran**: Required for compiling Fortran code (PMEMD core)

### 1.4 CUDA 12.0 Installation

**Important**: CUDA 12.0 is a hard requirement for Alchemd and must be installed manually.

#### 1.4.1 Check NVIDIA Drivers

```bash
# Check NVIDIA drivers
nvidia-smi

# If command not found or errors, install NVIDIA drivers first
```

#### 1.4.2 Install NVIDIA Drivers (if needed)

```bash
# Auto-detect and install recommended drivers
sudo ubuntu-drivers autoinstall

# Or manually install specific version
# sudo apt install -y nvidia-driver-535

# Reboot system
sudo reboot
```

#### 1.4.3 Download and Install CUDA 12.0

**Important Note**: The CUDA version must be compatible with the GCC compiler version. CUDA 12.0 is recommended for use with GCC 9-11 versions. Avoid using newer GCC versions (like GCC 12+) to compile older CUDA versions, as this may cause compatibility issues.

```bash
# Visit NVIDIA official website to download CUDA 12.0
# Official link: https://developer.nvidia.com/cuda-downloads
# Select Linux -> x86_64 -> Ubuntu -> version -> runfile(local)
# Download to ~/Downloads directory

cd ~/Downloads
# Run installer (filename may vary by version)
# First perform normal CUDA installation
sudo sh cuda_12.0.0_*_linux.run

# In the installation interface:
# - Deselect "Driver" (if already installed)
# - Keep "CUDA Toolkit" selected
# - Confirm installation path as /usr/local/cuda-12.0

# If you need to install NVIDIA drivers additionally (after basic installation):
# sudo sh cuda_12.0.0_*_linux.run --silent --driver
```

#### 1.4.4 Configure CUDA Environment Variables

```bash
# Edit ~/.bashrc
vim ~/.bashrc

# Add the following lines to the end:
export CUDA_HOME=/usr/local/cuda-12.0
export PATH=$CUDA_HOME/bin:$PATH
export LD_LIBRARY_PATH=$CUDA_HOME/lib64:$LD_LIBRARY_PATH

# Save and reload configuration
source ~/.bashrc
```

**Verification**:
```bash
# Check CUDA installation
nvcc --version
nvidia-smi

# Both commands should run successfully and display version information
```

---

## Chapter 2: Python Environment Configuration

### 2.1 Conda/Mamba Installation

#### 2.1.1 Choose Installation Method

**Recommended: Mamba** (faster package manager):
```bash
# Visit official website to download the latest version of Mambaforge
# Official link: https://github.com/conda-forge/miniforge
# Select the version suitable for your system and download to ~/Downloads directory

# Install Mambaforge (example with x86_64 version)
cd ~/Downloads
bash Mambaforge-Linux-x86_64.sh

# Follow prompts:
# - Agree to license terms
# - Confirm installation path (default ~/mambaforge)
# - Choose shell initialization (recommend yes)
```

**Alternative: Miniconda**:
```bash
# Visit official website to download the latest version of Miniconda
# Official link: https://docs.conda.io/en/latest/miniconda.html
# Select the version suitable for your system and download to ~/Downloads directory

cd ~/Downloads
# Install Miniconda (example with x86_64 version)
bash Miniconda3-latest-Linux-x86_64.sh
```

#### 2.1.2 Restart Terminal and Verify

```bash
# Restart terminal or reload configuration
source ~/.bashrc

# Verify installation
mamba --version  # If Mamba was installed
# or
conda --version  # If Conda was installed
```

**Dependencies**:
- **Conda/Mamba**: Python package and environment management
- **Python 3.9.21**: Recommended Python version for Alchemd

### 2.2 Virtual Environment Creation

```bash
# Create Alchemd dedicated environment
mamba create -n alchemd python=3.9.21
# or use conda create -n alchemd python=3.9.21

# Activate environment
mamba activate alchemd
# or conda activate alchemd

# Verify environment
python --version  # Should display Python 3.9.21
which python      # Should display mambaforge/envs/alchemd/bin/python path
```

### 2.3 Python Dependencies Installation

#### 2.3.1 Get Alchemd Source Code

```bash
# Clone project (if not already done)
cd ~
git clone <repository-url> Alchemd_release
cd Alchemd_release
```

#### 2.3.2 Install Python Dependencies

```bash
# Ensure in alchemd environment
mamba activate alchemd

# Install dependencies (using requirements.txt)
pip install -r requirements.txt

# Verify key dependencies
python -c "import numpy; print('NumPy:', numpy.__version__)"
python -c "import rdkit; print('RDKit:', rdkit.__version__)"
python -c "import openmm; print('OpenMM:', openmm.__version__)"
```

**Dependencies**:
- **NumPy, SciPy**: Scientific computing foundations
- **RDKit**: Molecular processing and cheminformatics
- **OpenMM**: Molecular dynamics simulation engine
- **Pandas**: Data processing
- **Matplotlib**: Plotting and visualization

---

## Chapter 3: External Tools Configuration

### 3.1 WatVina Download and Configuration

```bash
# Ensure in project root directory
cd ~/Alchemd_release

# Create dependencies directory
mkdir -p dependencies/watvina

# Download WatVina
wget -O dependencies/watvina/watvina \
    https://github.com/biocheming/watvina/releases/download/v20241125/watvina

# Set executable permissions
chmod +x dependencies/watvina/watvina

# Verify
./dependencies/watvina/watvina --help
```

**Dependencies**:
- **WatVina**: Tool for molecular docking
- **Used for**: Protein-ligand interaction prediction

### 3.2 Related Project Dependencies

```bash
# Enter dependencies directory
cd dependencies

# Download Weighted_cc project
wget -O weighted_cc.zip \
    https://github.com/zlisysu/Weighted_cc/archive/refs/heads/main.zip
unzip weighted_cc.zip

# Download Pyautomd project
wget -O pyautomd.zip \
    https://github.com/phloglucinol/Pyautomd/archive/refs/heads/main.zip
unzip pyautomd.zip

# Install Pyautomd
pip install ./Pyautomd-main/

# Download and install Lomap (lomap2)
wget -O lomap.zip \
    https://github.com/OpenFreeEnergy/Lomap/archive/refs/heads/main.zip
unzip lomap.zip

# Install Lomap
pip install ./Lomap-main/

# Verify Lomap installation
python -c "import lomap; print('Lomap installed successfully')"

# Clean up download files
rm -f weighted_cc.zip pyautomd.zip lomap.zip
```

**Dependencies**:
- **Weighted_cc**: Analysis tool
- **Pyautomd**: Automated molecular preparation tool  
- **Lomap**: Perturbation map generation tool

---

## Chapter 4: PMEMD Compilation

### 4.1 Obtain Amber License

1. **Visit Amber website**: https://ambermd.org/
2. **Apply for academic license**:
   - Fill in institutional information
   - Provide supervisor/PI information
   - Wait for approval (usually 1-3 business days)
3. **Download license file**: Download PMEMD package after receiving confirmation email

### 4.2 PMEMD Source Preparation

```bash
# Create installation directory
mkdir -p ~/software
cd ~/software

# Extract PMEMD package (assuming downloaded to ~/Downloads)
tar -xjf ~/Downloads/pmemd24.tar.bz2
# This creates pmemd24 directory

# Verify directory structure
ls pmemd24/
# Should contain: build/, src/, etc.
```

### 4.3 Compilation Configuration and Execution

#### 4.3.1 Pre-compilation Setup

```bash
# Enter PMEMD directory
cd ~/software/pmemd24

# Check build directory
ls build/
# Should contain run_cmake script

# Apply Alchemd patches
cp ~/Alchemd_release/dependencies/AmberFiles/patch.sh src/pmemd/src/
cp ~/Alchemd_release/dependencies/AmberFiles/ti.F90.patch src/pmemd/src/
cd src/pmemd/src/
chmod +x patch.sh
./patch.sh
cd ../../..
```

#### 4.3.2 Configure Compilation Options

```bash
# Edit compilation script
cd build
cp run_cmake run_cmake.backup

# Modify run_cmake file
sed -i 's/DCUDA=FALSE/DCUDA=TRUE/g' run_cmake
sed -i 's/DDOWNLOAD_MINICONDA=TRUE/DDOWNLOAD_MINICONDA=FALSE/g' run_cmake

# Verify modifications
grep "DCUDA=TRUE" run_cmake
grep "DDOWNLOAD_MINICONDA=FALSE" run_cmake
```

#### 4.3.3 Execute Compilation

```bash
# Ensure in alchemd environment
mamba activate alchemd

# Set environment variables
export AMBERHOME=~/software/pmemd24
export CUDA_HOME=/usr/local/cuda-12.0

# Clean previous builds (if any)
if [ -f "./clean_build" ]; then
    echo "y" | ./clean_build
else
    echo "y" | make clean || true
fi

# Run CMake configuration
./run_cmake

# Start compilation (using 4 cores)
make install -j4

# Compilation may take 30-120 minutes depending on system performance
```

**Troubleshooting**:
- If memory insufficient, reduce parallelism: `make install -j2`
- If CUDA errors, check CUDA_HOME environment variable
- If GCC errors, confirm version >= 9

#### 4.3.4 Verify Compilation Results

```bash
# Check generated binary files
ls ~/software/pmemd24/bin/

# Should contain:
# - pmemd
# - pmemd.cuda
# - pmemd.cuda_SPFP

# Test key binary files
~/software/pmemd24/bin/pmemd.cuda_SPFP --help
```

### 4.4 Environment Variable Setup

#### 4.4.1 Automatic Environment Setup

```bash
# Use Alchemd provided environment setup script
cd ~/Alchemd_release
./install/setup_amber_env.sh ~/software/pmemd24
```

#### 4.4.2 Manual Environment Setup (Optional)

```bash
# Edit conda environment activation script
mkdir -p ~/mambaforge/envs/alchemd/etc/conda/activate.d/
cat > ~/mambaforge/envs/alchemd/etc/conda/activate.d/amber_init.sh << 'EOF'
#!/bin/bash
# PMEMD environment activation
export AMBERHOME="$HOME/software/pmemd24"
source "$AMBERHOME/amber.sh"
echo "PMEMD environment activated: $AMBERHOME"
EOF

chmod +x ~/mambaforge/envs/alchemd/etc/conda/activate.d/amber_init.sh

# Create deactivation script
mkdir -p ~/mambaforge/envs/alchemd/etc/conda/deactivate.d/
cat > ~/mambaforge/envs/alchemd/etc/conda/deactivate.d/amber_init.sh << 'EOF'
#!/bin/bash
# PMEMD environment deactivation
if [[ "$AMBERHOME" == "$HOME/software/pmemd24" ]]; then
    unset AMBERHOME
    echo "PMEMD environment deactivated"
fi
EOF

chmod +x ~/mambaforge/envs/alchemd/etc/conda/deactivate.d/amber_init.sh
```

#### 4.4.3 Test Environment Setup

```bash
# Reactivate environment
mamba deactivate
mamba activate alchemd

# Check AMBERHOME
echo $AMBERHOME
# Should display: /home/username/software/pmemd24

# Test PMEMD tools
which pmemd.cuda_SPFP
pmemd.cuda_SPFP --help
```

---

## Chapter 5: Alchemd Configuration

### 5.1 Configuration File Setup

Alchemd uses TOML format configuration files to manage paths for various tools and dependencies. We provide an automated configuration tool to simplify this process.

#### 5.1.1 Using Automatic Configuration Tool

```bash
# Ensure in alchemd environment
mamba activate alchemd
cd ~/Alchemd_release

# Run configuration initialization tool
python core/init_setup.py
```

This tool will:
- Auto-detect project paths
- Validate dependency tool availability
- Generate `core/configs.toml` configuration file
- Provide configuration integrity testing

#### 5.1.2 Interactive Configuration (Advanced Users)

If you need to customize path settings:

```bash
# Interactive configuration with customizable tool paths
python core/init_setup.py --interactive
```

#### 5.1.3 Verify Configuration

After configuration completion, check the generated configuration file:

```bash
# View configuration file content
cat core/configs.toml

# Should contain content similar to:
# [paths]
# AlchemdToolkit_Path = "/home/username/Alchemd_release"
# CAR_Path = "/home/username/Alchemd_release/dependencies/CAR-FEP"
# ...
```

### 5.2 Testing and Verification

#### 5.2.1 Basic Functionality Testing

```bash
# Ensure in alchemd environment
mamba activate alchemd
cd ~/Alchemd_release

# Test core modules
python core/prepare_file.py --help
python core/run_md.py --help
python core/analyze_result.py --help

# Test GUI (optional)
python gui/gui.py
```

#### 5.2.2 PMEMD Integration Testing

**Important Note**: Please strictly follow the sequence below to ensure proper environment setup:

```bash
# 1. Ensure you're in alchemd environment (should be automatic after compilation)
echo $CONDA_DEFAULT_ENV  # Should display alchemd

# 2. First set up PMEMD environment variable integration
./install/setup_amber_env.sh /path/to/pmemd24

# 3. Reactivate environment to load new environment variables
conda deactivate
conda activate alchemd

# 4. Run environment test script
./install/test_amber_env.sh

# 5. Check output, should display most tests passed
```

**Troubleshooting**:
- If tests fail, confirm PMEMD compilation completed without errors
- Check if files like `pmemd.cuda_SPFP` exist in `/path/to/pmemd24/bin/` directory
- Confirm CUDA environment variables are properly set

#### 5.2.3 Example Data Testing

Alchemd provides example data for testing the complete workflow:

```bash
# View example data
ls example/
# Should display:
# ligands.sdf  pairs.txt  protein.pdb

# Test preparation step using example data
python core/prepare_file.py -p example/protein.pdb -l example/ligands.sdf -o test_output --cs-align -d -auto

# Check output directory
ls test_output/
# Should generate prepared input files
```

**Example Data Description**:
- `protein.pdb`: Example protein structure file
- `ligands.sdf`: SDF file containing multiple ligand molecules
- `pairs.txt`: Predefined ligand pair file

**Note**: Complete MD simulation testing requires compiled PMEMD and typically takes significant computational time; recommended for production environments.

### 5.3 Configuration Troubleshooting

#### 5.3.1 Common Configuration Tool Issues

**Q: init_setup.py reports missing dependencies**
```bash
# Check dependencies directory
ls dependencies/

# Ensure necessary tar files are extracted
cd dependencies
tar -xf CAR-FEP.tar
tar -xf genambRBFE.tar
tar -xf AlchemConvTools.tar
```

**Q: AMBERHOME environment variable not set**
```bash
# Ensure PMEMD is compiled and environment is set
source /path/to/pmemd24/amber.sh
echo $AMBERHOME
```

**Q: Configuration file generation failed**
```bash
# Check permissions
ls -la core/

# Manually create directory (if needed)
mkdir -p core

# Re-run configuration
python core/init_setup.py --auto
```

#### 5.3.2 Configuration Validation Failures

**Q: Dependency tool tests failed**
```bash
# Manually test each tool
./dependencies/watvina/watvina --help
python dependencies/Weighted_cc-main/wcc_main.py -h
python dependencies/CAR-FEP/segmented_converge_control.py -h
```

**Q: Reconfigure system**
```bash
# Remove old configuration and regenerate
rm -f core/configs.toml
python core/init_setup.py --interactive
```

---

## Summary

Congratulations! You have completed the manual installation of Alchemd. You should now have:

- ✅ Fully configured Ubuntu system environment
- ✅ GCC 9+ compiler toolchain
- ✅ CUDA 12.0 support
- ✅ Python 3.9.21 and all required dependencies
- ✅ Fully compiled PMEMD (with CUDA support)
- ✅ Properly configured Alchemd environment

You can now start using Alchemd for molecular dynamics simulations and free energy calculations!

### Next Steps

- Read Alchemd user guides for specific usage instructions
- Prepare your protein and ligand input files
- Start your first FEP calculation project

### Recommendations

- Regularly backup your configurations and results
- Keep CUDA drivers updated
- Monitor system resource usage during calculations

---

## Chapter 6: Optional Components Configuration

### 6.1 NSAA Module Dependencies (Optional)

The NSAA (Non-Standard Amino Acid) preparation module is an optional functional component of Alchemd, specifically designed to process protein systems containing non-standard amino acids. This module requires Gaussian quantum chemistry software support.

#### 6.1.1 Gaussian Software Requirements

**Important Note**: Gaussian software must be obtained and installed independently by users. Alchemd does not provide the Gaussian software itself.

**Supported Gaussian Versions**:
- Gaussian 03 (g03)
- Gaussian 09 (g09)
- Gaussian 16 (g16) - Recommended version

**License Requirements**:
- Academic users: Valid academic license required
- Commercial users: Commercial license purchase required
- Application website: https://gaussian.com/

#### 6.1.2 Gaussian Installation Guidance

**Note**: The following is general guidance only. For specific installation steps, please refer to the official Gaussian documentation.

```bash
# Create Gaussian installation directory
sudo mkdir -p /opt/gaussian
cd /opt/gaussian

# Extract Gaussian software package (assuming license obtained and downloaded)
# tar -xvf g16.tar  # Example, actual filename may differ

# Set environment variables
export g16root="/opt/gaussian"
export GAUSS_SCRDIR="/tmp"
export GAUSS_MEMDEF="2GB"

# Load Gaussian environment
source $g16root/g16/bsd/g16.profile
```

#### 6.1.3 NSAA Configuration in Alchemd

```bash
# Enter Alchemd project directory
cd ~/Alchemd_release

# Edit NSAA configuration template
vim core/templates/nsaa_prep.txt
```

Locate the Gaussian configuration on line 11:
```
gaussian = g03
```

Modify according to your installed Gaussian version:
```
# Example configurations
gaussian = g16           # If using Gaussian 16
gaussian = g09           # If using Gaussian 09
gaussian = /opt/gaussian/g16/g16  # Using full path
```

#### 6.1.4 Verify NSAA Module Configuration

```bash
# Activate alchemd environment
mamba activate alchemd

# Test Gaussian availability
which g16  # Or your configured Gaussian command
g16 < /dev/null  # Test if Gaussian can start normally

# Test NSAA module
python -c "from core.prepare_input_tools.nsaa_preparer import NSAAProcessor; print('NSAA module available')"
```

#### 6.1.5 Usage Instructions

The NSAA module is used through the `-npd` parameter of prepare_file.py:

```bash
# Example processing protein containing non-standard amino acids
python core/prepare_file.py \
    -p protein_with_nsaa.pdb \
    -l ligands.sdf \
    -npd ./nsaa_parameters \
    -o ./nsaa_project \
    --cs-align -d -auto
```

#### 6.1.6 Troubleshooting

**Q: Gaussian command not found**
```bash
# Check environment variables
echo $g16root
echo $PATH | grep gaussian

# Reload Gaussian environment
source /opt/gaussian/g16/bsd/g16.profile
```

**Q: Gaussian license issues**
- Confirm license file location is correct
- Check if license is within validity period
- Contact Gaussian technical support

**Q: Memory insufficient**
```bash
# Adjust Gaussian memory settings
export GAUSS_MEMDEF="4GB"  # Adjust according to system configuration
```

**Important Notes**:
1. NSAA module is only needed when processing proteins containing non-standard amino acids
2. Standard amino acid FEP calculations do not require Gaussian
3. Quantum chemistry calculations are time-intensive; recommended for high-performance environments
4. Regularly check Gaussian license validity

---

## Chapter 7: Frequently Asked Questions (QA)

### 7.1 Conda/Mamba Related Issues

#### Q: Encountering CondaToSNonInteractiveError

**Problem Description**: "Terms of Service have not been accepted for the following channels. Please accept or remove them before proceeding"

**Solution**:
```bash
# Method 1: Accept terms of service for specific channels
conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/main
conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/r

# Method 2: View detailed help and handle manually
conda tos --help
```

**Alternatively**: Consult Conda official documentation to learn how to accept terms of service.

### 7.2 Compilation Related Issues

#### Q: Segment Fault during compilation in WSL environment

**Problem Description**: Random segment faults or compilation failures in WSL environment

**Solution**:
```bash
# Method 1: Adjust compilation cores (default is 4)
./install/build_amber.sh /path/to/pmemd24_src 2  # Use 2 cores
./install/build_amber.sh /path/to/pmemd24_src 4  # Use 4 cores

# Method 2: Re-run compilation
# WSL sometimes has random segment faults, re-compilation usually solves it
./install/build_amber.sh /path/to/pmemd24_src
```

**Recommendation**: In WSL environment, use fewer compilation cores (2-4 cores) to avoid memory shortage issues.

### 7.3 Installation Verification Issues

#### Q: Verification fails after build_amber.sh completion

**Problem Description**: Compilation completed but pmemd tools cannot be found or used

**Correct Verification Process**:
```bash
# 1. Stay in alchemd environment after compilation
# (build_amber.sh script will automatically keep environment activated)

# 2. Set up environment variables
./install/setup_amber_env.sh /path/to/pmemd24

# 3. Reactivate environment to load new environment variables
conda deactivate
conda activate alchemd

# 4. Verify installation
./install/test_amber_env.sh

# 5. Test PMEMD tools
which pmemd.cuda_SPFP
pmemd.cuda_SPFP --help
```

**Note**: Do not exit conda environment immediately after compilation script completion. The correct sequence is to first set up environment variables, then reactivate the environment.

---

## Appendix: Troubleshooting

### A.1 Compilation Issues

**Q: "CUDA not found" error during compilation**
```bash
# Check CUDA installation paths
ls /usr/local/cuda*

# Set correct CUDA_HOME
export CUDA_HOME=/usr/local/cuda-12.0

# Ensure path is in PATH
export PATH=$CUDA_HOME/bin:$PATH
```

**Q: GCC version too old**
```bash
# Install newer version
sudo apt install gcc-11 g++-11 gfortran-11

# Update default version
sudo update-alternatives --install /usr/bin/gcc gcc /usr/bin/gcc-11 60
```

### A.2 Runtime Issues

**Q: Cannot find pmemd.cuda_SPFP**
```bash
# Check if file exists
ls $AMBERHOME/bin/pmemd*

# If not exists, recompile PMEMD
cd ~/software/pmemd24/build
make install -j4
```

**Q: GPU cannot be used**
```bash
# Check NVIDIA drivers
nvidia-smi

# Test CUDA
cd $CUDA_HOME/samples/1_Utilities/deviceQuery
make
./deviceQuery
```

### A.3 Environment Issues

**Q: Conda environment confusion**
```bash
# Rebuild environment
mamba deactivate
mamba env remove -n alchemd
mamba create -n alchemd python=3.9.21
mamba activate alchemd
pip install -r requirements.txt
```

**Q: Permission issues**
```bash
# Check file permissions
ls -la ~/software/pmemd24/bin/

# Fix permissions (if needed)
chmod +x ~/software/pmemd24/bin/*
```

### A.4 Getting Help

If you encounter issues not covered in this guide:

1. **Check log files**: Review detailed error information from compilation or runtime
2. **Consult official documentation**: Amber official documentation and Alchemd project documentation
3. **Community help**: Get assistance through GitHub Issues or relevant forums
4. **Fallback to automated scripts**: If manual installation encounters difficulties, try using the automated installation scripts provided by the project

---

---

## License and Citation

This documentation is part of Alchemd. For license information and citation requirements, please refer to [LICENSE.md](../LICENSE.md).

### Citation Format

If you use Alchemd in your research, please cite:

```bibtex
@software{alchemd2025,
  title={Alchemd: A toolkit for alchemical free energy perturbation calculations},
  year={2025},
  url={https://github.com/ZheLi-Lab/Alchemd-toolkit}
}
```

**Methodological References**:
- CAR method: https://pubs.acs.org/doi/10.1021/acs.jctc.4c00939
- cs-fep implementation: https://doi.org/10.1016/j.apsb.2024.06.021

---

**Document Information**
- Version: 1.0
- Last Updated: December 2024
- Authors: Alchemd Development Team
- License: CC BY-NC-SA 4.0 (Academic Use)