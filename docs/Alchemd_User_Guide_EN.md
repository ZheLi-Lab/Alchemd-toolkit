# Alchemd User Guide v1.2

---
**License**: This documentation is part of Alchemd, licensed under CC BY-NC-SA 4.0 for academic use. Commercial licensing available separately. See [LICENSE.md](../LICENSE.md) for complete terms.
---

## Table of Contents
1. [Introduction](#introduction)
2. [Quick Start](#quick-start)
3. [Core Modules](#core-modules)
4. [Complete Workflow](#complete-workflow)
5. [FAQ](#faq)

## Introduction

Alchemd is a Python-based CADD (Computer-Aided Drug Design) application specialized in free energy perturbation (FEP) calculations using alchemical transformations. The toolkit provides a complete workflow from input preparation through molecular dynamics simulations to results analysis.

### Project Structure

```
Alchemd/
├── core/                    # Core modules
│   ├── prepare_file.py     # Input preparation
│   ├── run_md.py          # Molecular dynamics simulation
│   └── analyze_result.py  # Results analysis
├── gui/                    # Graphical interface
├── dependencies/           # External dependencies
└── docs/                   # Documentation
```

### Computational Methods
- **Amber**: Molecular dynamics engine
- **CAR (Convergence-Adaptive Roundtrip Method)**: Adaptive convergence method, enabled by default

## Quick Start

### Standard Workflow

```bash
# 1. Input preparation
python core/prepare_file.py -auto -p protein.pdb --cs-align -d -l ligands.sdf -o ./output

# 2. MD simulation (using CAR method)
python core/run_md.py -ln 50 -p ./output/prepare -auto

# 3. Results analysis
python core/analyze_result.py -auto -w ./output/run
```

## Core Modules

## 1. prepare_file.py - Input Preparation

### Basic Syntax
```bash
python core/prepare_file.py [options]
```

### Main Parameters

#### Required Parameters

**`-p, --protein PROTEIN`**
- **Function**: Specify protein PDB file path
- **Type**: String
- **Format**: .pdb
- **Example**: `-p protein.pdb`
- **Description**: Protein structure file to be used as receptor

**`-o, --output-dir OUTPUT_DIR`**  
- **Function**: Specify output directory
- **Type**: String
- **Example**: `-o ./output`
- **Description**: Location to save all processed files, must be specified

#### Ligand Input Methods (choose one)

**Method 1: Input ligands from SDF file**

**`-l, --ligands LIGANDS`**
- **Function**: Specify ligand SDF file path
- **Type**: String
- **Format**: .sdf
- **Example**: `-l ligands.sdf`
- **Requirements**: 
  - Must contain bond information
  - Must include explicit and implicit hydrogens
  - Each molecule requires unique name

**Method 2: Load from existing folder**

**`-af, --advance-folder ADVANCE_FOLDER`**
- **Function**: Specify pre-organized prepare folder path
- **Type**: String
- **Example**: `-af ./existing_prepare`
- **Description**: Suitable when pair molecular PDB files already exist

#### Ligand Filtering and Pairing

**`-r, --ref-ligand-list REF_LIGAND_LIST`**
- **Function**: Specify reference ligand list file
- **Type**: String
- **Format**: .lst
- **Example**: `-r ref_ligands.lst`
- **Description**: 
  - One ligand name per line
  - Names must correspond to molecules in SDF file
  - If not specified, all ligands are treated as reference ligands

**`-pl, --pair-list PAIR_LIST`**
- **Function**: Specify ligand pair list file
- **Type**: String
- **Format**: .lst
- **Example**: `-pl pairs.lst`
- **Description**: 
  - One ligand pair per line, format: `ligand1-ligand2`
  - If not specified, LOMAP2 automatically generates pairs

#### Molecular Processing Options

**`--cs-align`**
- **Function**: Use cs=align method for molecule generation
- **Type**: Boolean
- **Default**: False
- **Recommendation**: **Strongly recommended**
- **Description**: Generate 3D conformations of paired molecules from reference molecules

**`-d, --dock`**
- **Function**: Enable molecular docking
- **Type**: Boolean  
- **Default**: False
- **Recommendation**: **Recommended**
- **Description**: Perform protein-ligand docking after molecule generation

**`-npp, --not-prepare-protein`**
- **Function**: Skip protein preprocessing
- **Type**: Boolean
- **Default**: False
- **Recommendation**: Usually keep default (perform protein preprocessing)

#### Ligand Pair Processing Options

**`-sp, --strict-pair-list`**
- **Function**: Strictly follow user-defined ligand pair list
- **Type**: Boolean
- **Default**: False
- **Condition**: Only effective when using both `--cs-align` and `-p`
- **Description**: 
  - Enabled: Process pairs strictly in given order
  - Disabled: Automatically optimize pair processing order

**`-f, --flip`**
- **Function**: Set FEP alchemical transformation direction from large to small molecules
- **Type**: Boolean
- **Default**: False
- **Description**: Automatically adjust pair orientation based on van der Waals volume to ensure large-to-small transformation

**`-as, --allow-sep`**
- **Function**: Allow multiple disconnected ligand pair groups
- **Type**: Boolean
- **Default**: False
- **Description**: Handle cases with multiple independent ligand pair networks

#### Non-standard Amino Acid Support

**`-npd, --nsaa-params-dir NSAA_PARAMS_DIR`**
- **Function**: Specify directory containing pre-generated non-standard amino acid force field parameters
- **Type**: String
- **Example**: `-npd ./nsaa_parameters`
- **Description**: Path to directory containing pre-generated .frcmod and .prepi files for non-standard amino acids. Parameters must be generated beforehand using nsaa_preparer.py
- **Important**: This parameter consumes existing parameters; it does not generate them. Parameter generation is a separate preprocessing step
- **Detailed Configuration**: See [NSAA Preparation Module](#nsaa-preparation-module) section for complete configuration steps and Gaussian dependency setup

#### Automation Options

**`-auto`**
- **Function**: Automatically prepare topology and coordinate files
- **Type**: Boolean
- **Default**: False
- **Recommendation**: **Recommended**
- **Description**: Automatically execute preparation steps after pair file preparation, instead of generating a job list only.

### Output Files

Upon successful execution, the following directories and files are generated:

**`preprocess/`** - Preprocessing directory
- Stores preprocessed protein and ligand files
- Contains intermediate files for molecular docking and coordinate generation

**`prepare/`** - Preparation working directory  
- One subdirectory per ligand pair (e.g., `A-B/`)
- Each subdirectory contains:
  - `A.pdb` - Reference ligand PDB file
  - `B.pdb` - Target ligand PDB file
  - `submit.sh` - Preparation script

**`prepare_queue.lst`** - Preparation task queue
- Each line contains working directory and execution script path
- Format: `/path/to/pair/dir /path/to/submit.sh`

**`info.pickle`** - Metadata file
- Contains ligand pair information and molecular data

### Usage Examples

**Basic usage**:
```bash
python core/prepare_file.py -p protein.pdb -l ligands.sdf -o ./output --cs-align -d -auto
```

**Specify reference ligands**:
```bash
python core/prepare_file.py -p protein.pdb -l ligands.sdf -r ref_list.lst -o ./output --cs-align -d -auto
```

**Use custom ligand pairs**:
```bash
python core/prepare_file.py -p protein.pdb -l ligands.sdf -pl pairs.lst -sp -o ./output --cs-align -d -auto
```

### Input File Format Specifications

This section provides detailed format requirements and usage guidelines for three **optional important** parameters in `prepare_file.py`.

#### pairs.lst (--pair-list) Format Specification

**Purpose**: Manually specify ligand pair list, replacing automatic LOMAP2 pair generation

**File Format**:
- Plain text file, recommended `.txt` or `.lst` extension
- One ligand pair per line, format: `ligandA-ligandB`
- Use hyphen `-` to separate two ligand names
- Supports empty lines and comments (lines starting with `#` are ignored)
- File encoding: Recommend UTF-8 encoding

**Example File Content**:
```
# Example ligand pair list
mol1-mol2
mol2-mol3
mol3-mol4
mol1-mol4
mol5-mol6
```

**Important Requirements**:
- Ligand names must exactly match the `_Name` property of molecules in the SDF file
- Ligand pairs have bidirectional symmetry: `mol1-mol2` and `mol2-mol1` are treated as identical pairs
- When used with `--strict-pair-list` parameter, pairs are processed strictly in the order specified in the file
- Each ligand pair only needs to be specified once; the program automatically handles bidirectional relationships

#### ref.lst (--ref-ligand-list) Format Specification

**Purpose**: Specify the list of molecules to be used as reference ligands for coordinate alignment and structure generation

**File Format**:
- Plain text file, recommended `.txt` or `.lst` extension
- One ligand name per line
- Supports empty lines and comments (lines starting with `#` are ignored)
- File encoding: Recommend UTF-8 encoding

**Example File Content**:
```
# Example reference ligand list
mol1
mol3
mol5
```

**Important Requirements**:
- Ligand names must exactly match the `_Name` property of molecules in the SDF file
- Reference ligands provide reliable 3D conformations; recommend selecting molecules with experimental structures or high-quality conformations
- Molecules not listed in this file will be marked as active ligands and require conformation generation and optimization
- If this file is not provided, all molecules are treated as reference ligands

#### Advanced Folder Mode (--advance-folder) Format Specification

**Purpose**: Use pre-processed ligand pair PDB files, skipping automatic molecular alignment and conformation generation steps

**Directory Structure Requirements**:
```
advance_folder/
├── mol1-mol2/
│   ├── mol1.pdb
│   └── mol2.pdb
├── mol2-mol3/
│   ├── mol2.pdb
│   └── mol3.pdb
├── mol3-mol4/
│   ├── mol3.pdb
│   └── mol4.pdb
└── mol5-mol6/
    ├── mol5.pdb
    └── mol6.pdb
```

**Important Requirements**:
- Each subfolder name must follow the `ligandA-ligandB` format (must contain hyphen `-`)
- Each subfolder must contain the corresponding two PDB files
- PDB file names must exactly match the molecule names in the pair name
- PDB files must contain complete atomic coordinate information and correct residue information
- The system will automatically extract molecular information and charge information from PDB files

**Use Cases**:
- High-quality molecular conformations have been generated through external tools (such as molecular docking, quantum chemistry calculations)
- Specific molecular orientations or conformations are needed for FEP calculations
- Want to skip automatic conformation generation steps to save computational time

**Important Notes**:
- This mode bypasses all molecular preprocessing steps; users must ensure input PDB file quality
- No molecular docking or coordinate alignment will be performed; provided coordinates are used directly

#### Format Validation Recommendations

**File Encoding**: Recommend using UTF-8 encoding for all text files

**Path Usage**:
- Recommend using absolute paths to avoid path resolution issues
- Ensure all referenced files and directories exist and are accessible
- Avoid special characters (spaces, non-ASCII characters) in file paths and molecule names

**Naming Conventions**:
- Molecule names should use English letters, numbers, and underscores; avoid special symbols except hyphens in file names
- Folder names should follow system naming conventions and avoid excessively long paths

**Naming Consistency Check**:
```bash
# Verify molecule names in SDF file
grep "_Name" ligands.sdf

# Check ligand pair file format
cat pairs.lst | grep -v "^#" | grep -v "^$"

# Verify reference ligand list
cat ref.lst | grep -v "^#" | grep -v "^$"

# Check for special characters in files
file -bi pairs.lst  # Check file encoding
```

**Common Errors and Solutions**:

1. **Molecule name mismatch**: Check if `_Name` property in SDF file exactly matches names in list files
2. **File format error**: Ensure UTF-8 encoding is used and check for hidden characters
3. **Path issues**: Use absolute paths and ensure files exist
4. **Permission issues**: Ensure files have read permissions

---

## 2. run_md.py - Molecular Dynamics Simulation

### Basic Syntax
```bash
python core/run_md.py [options]
```

### Main Parameters

#### Required Parameters

**`-p, --preparation_dir PREPARATION_DIR`**
- **Function**: Specify prepare directory path
- **Type**: String
- **Example**: `-p ./output/prepare`
- **Description**: 
  - Must be a directory processed by prepare_file.py
  - Ensure all tasks in prepare_queue.lst have been completed

#### Lambda Parameter Configuration

**`-ln, --lambda_nums LAMBDA_NUMS`**
- **Function**: Set number of Lambda windows
- **Type**: Integer
- **Default**: 50
- **Example**: `-ln 50`
- **Description**: 
  - Actual Lambda points = `lambda_nums + 1`
  - Example: `-ln 50` produces 51 Lambda values (1.0, 0.98, ..., 0.02, 0.0)
  - Recommended: 50 (corresponding to Lambda interval 0.02)

#### Computational Engine Configuration

**`-m, --mode MODE`**
- **Function**: Specify molecular dynamics engine
- **Type**: String
- **Default**: `Amber`
- **Example**: `-m Amber`
- **Description**: Use Amber engine with CAR method (default setting)

#### CAR Adaptive Convergence Method Configuration

**`-ci, --CAR_input CAR_INPUT`**
- **Function**: Specify custom CAR input file
- **Type**: String
- **Format**: .txt
- **Example**: `-ci custom_car_settings.txt`
- **Description**: Use default CAR configuration template if not specified

#### Automation Options

**`-auto`**
- **Function**: Automatically run MD after preparation completion
- **Type**: Boolean
- **Default**: False
- **Recommendation**: **Recommended**
- **Description**: Automatically execute all MD tasks in local environment in **non-parallel** way.

### Output Files

Upon successful execution, generates in the parent directory of preparation_dir:

**`run/`** - MD working directory
- One subdirectory per ligand pair
- Each ligand pair contains 4 edges:
  - `lM2A/` - Alchemical transformation M→A in solution
  - `lM2B/` - Alchemical transformation M→B in solution  
  - `cM2A/` - Alchemical transformation M→A in complex
  - `cM2B/` - Alchemical transformation M→B in complex

**Each edge directory contains**:
- Topology files (`.prmtop`)
- Coordinate files (`.prmcrd`)
- Lambda configuration files (`.json`)
- Execution scripts (`submit.sh`)
- Input configuration files

**`run_all_queue.lst`** - MD task queue
- Contains all MD calculation tasks
- Format: `/path/to/edge/dir /path/to/submit.sh`

### Usage Examples

**Standard mode**:
```bash
python core/run_md.py -ln 50 -p ./output/prepare -auto
```

**Custom CAR configuration**:
```bash
python core/run_md.py -ln 50 -p ./output/prepare -ci my_car_config.txt -auto
```

---

## 3. analyze_result.py - Results Analysis

### Basic Syntax
```bash
python core/analyze_result.py [options]
```

### Main Parameters

#### Required Parameters

**`-w, --work_dir WORK_DIR`**
- **Function**: Specify MD results directory path
- **Type**: String
- **Default**: `./work`
- **Example**: `-w ./output/run`
- **Description**: Path to run directory generated by run_md.py

#### Computational Method Configuration

**`-m, --mode MODE`**
- **Function**: Specify MD engine type
- **Type**: String
- **Default**: `Amber`
- **Description**: Must be consistent with run_md.py settings

#### Automation Options

**`-auto`**
- **Function**: Automatically calculate energies and perform statistical analysis
- **Type**: Boolean
- **Default**: False
- **Recommendation**: **Recommended**
- **Description**: 
  - Automatically prepare calculation files
  - Execute energy calculations in parallel
  - Generate statistical analysis reports

### Parallel Processing Configuration

The program uses the following configuration for parallel processing:
- **Maximum worker processes**: `min(cpu_count, 4)`
- **Maximum retries**: 3
- **Task timeout**: 600 seconds

**Custom parallel configuration**: Users can modify the `MAX_WORKER_COUNT` parameter in the script to optimize performance based on computational resources. This parameter controls the maximum number of worker processes for both file preparation and calculation tasks.

### Output Files

Upon successful execution, generates in the parent directory of work_dir:

**`calculate/`** - Calculation working directory
- Contains energy calculation files for all ligand pairs
- Input and output files for each calculation unit

**`results/`** - Analysis results directory
- **`results.csv`**: Energy analysis results for all ligand pairs
  - Contains ΔΔG values, uncertainties, convergence indicators
- **`ene_detail.csv`**: Detailed energy decomposition analysis
  - Contains detailed energy data for each Lambda window
- **`summary.pdf`**: Complete analysis report
  - Includes energy plots, convergence analysis, statistical information
- **`summary_abstract.pdf`**: Simplified analysis report
  - Summary version of key results and plots

**Optional output files**:
- **`sus_pairs.lst`**: Suspicious ligand pairs list
  - Contains pairs with poor convergence or abnormal results
- **`err_pairs.lst`**: Error ligand pairs list
  - List of pairs that failed calculation
- **`calculate_queue.lst`**: Calculation task queue

### Usage Examples

**Standard analysis workflow**:
```bash
python core/analyze_result.py -auto -w ./output/run
```

**Manual analysis control**:
```bash
# Prepare calculation files only
python core/analyze_result.py -w ./output/run

# Manual execution of calculation and statistics afterwards
```

---

## NSAA Preparation Module

### Module Overview

The Non-Standard Amino Acid (NSAA) preparation module is an independent preprocessing tool within the Alchemd toolkit, designed to generate force field parameters for non-standard amino acids through quantum mechanical calculations. This module operates as a standalone preprocessing step, completely separate from the main FEP calculation workflow.

**Module Location**: `core/prepare_input_tools/nsaa_preparer.py`

**Key Characteristics**:
- Independent standalone preprocessing tool
- Requires separate execution before main workflow
- Generates quantum mechanically derived force field parameters
- Outputs structured parameter directories for consumption by prepare_file.py

### Computational Dependencies

The NSAA module requires several external computational tools:

**Required Dependencies**:
- **Gaussian** (g03/g09/g16): Quantum chemistry calculations for electronic structure
- **pyautomd**: Automated molecular dynamics parameter generation interface
- **AmberTools**: antechamber, parmchk2, tleap for force field parameter processing
  - **Important**: Unlike the main installation guide which only requires PMEMD, NSAA functionality requires a complete AmberTools compilation
  - Users must compile AmberTools separately from the main PMEMD-only installation
  - AmberTools compilation includes the necessary parameter generation utilities
- **GAUSS_SCRDIR**: Environment variable pointing to Gaussian temporary file directory

### Two-Step NSAA Workflow

The NSAA processing follows a mandatory two-step procedure that cannot be bypassed or integrated:

#### Step 1: Independent NSAA Parameter Generation

```bash
python core/prepare_input_tools/nsaa_preparer.py protein_with_nsaa.pdb -o ./nsaa_parameters
```

**Process Description**:
1. **Structure Analysis**: Scans PDB file to identify non-standard amino acid residues
2. **Residue Extraction**: Extracts NSAA and neighboring residues for context preservation
3. **Charge Calculation**: Computes formal charges using RDKit molecular descriptors
4. **Quantum Calculation**: Invokes pyautomd → Gaussian for electronic structure calculations
5. **Parameter Generation**: Produces .frcmod (force field modification) and .prepi (residue topology) files
6. **Directory Organization**: Outputs parameters to `nsaa_parameters/NSAA_NAME/` structure

**Output Structure**:
```
nsaa_parameters/
├── NSAA_NAME_1/
│   ├── NSAA_NAME_1.frcmod
│   └── NSAA_NAME_1.prepi
└── NSAA_NAME_2/
    ├── NSAA_NAME_2.frcmod
    └── NSAA_NAME_2.prepi
```

#### Step 2: Parameter Usage in Main Workflow

```bash
python core/prepare_file.py -p protein.pdb -l ligands.sdf -npd ./nsaa_parameters -o ./output
```

### Configuration Requirements

#### Environment Configuration

**Gaussian Environment Setup**:
```bash
# Set Gaussian scratch directory
export GAUSS_SCRDIR="/path/to/gaussian/scratch"

# Ensure Gaussian executable is accessible
export g16root="/path/to/gaussian"
source $g16root/g16/bsd/g16.profile
```

#### Template Modification

The NSAA preparation requires template configuration for Gaussian integration:

**Template File**: `core/templates/nsaa_prep.txt`

**Required Modifications**:
- Line 11: Specify Gaussian version (g03/g09/g16)
- Scratch directory paths for temporary file management
- Computational resource allocation parameters

### Usage Examples

#### Complete Two-Step Process

```bash
# Step 1: Generate NSAA parameters (independent preprocessing)
python core/prepare_input_tools/nsaa_preparer.py \
    protein_with_nsaa.pdb \
    -o ./nsaa_parameters \
    --no-delete  # Keep temporary files for debugging

# Step 2: Use parameters in main workflow
python core/prepare_file.py \
    -auto \
    -p protein_with_nsaa.pdb \
    -l ligands.sdf \
    -npd ./nsaa_parameters \
    --cs-align -d \
    -o ./fep_project
```

#### Multiple NSAA Systems

```bash
# Process different proteins with NSAAs
python core/prepare_input_tools/nsaa_preparer.py protein_A.pdb -o ./nsaa_params_A
python core/prepare_input_tools/nsaa_preparer.py protein_B.pdb -o ./nsaa_params_B

# Use parameters in respective workflows
python core/prepare_file.py -p protein_A.pdb -l ligands_A.sdf -npd ./nsaa_params_A -o ./project_A
python core/prepare_file.py -p protein_B.pdb -l ligands_B.sdf -npd ./nsaa_params_B -o ./project_B
```

### Troubleshooting Guide

#### Common Issues and Solutions

**Issue**: "Environment variable GAUSS_SCRDIR is not set"
- **Solution**: Configure Gaussian scratch directory: `export GAUSS_SCRDIR="/tmp/gaussian"`

**Issue**: "pyautomd command not found"
- **Solution**: Ensure pyautomd is installed and accessible in PATH

**Issue**: "Gaussian calculation timeout"
- **Solution**: Increase timeout values or optimize computational resources

**Issue**: "No non-standard amino acids found"
- **Solution**: Verify PDB file contains actual non-standard residues not in standard list

#### Validation Procedures

**Parameter Validation**:
1. Verify .frcmod files contain reasonable force constants
2. Check .prepi files for proper atomic connectivity
3. Test parameters with small tleap validation runs
4. Compare generated charges with expected chemical behavior

### Technical Limitations

**Current Limitations**:
- Manual capping required for peptide-linked NSAAs
- AM1-BCC charges used by default (RESP charges require additional QM calculations)
- Sequential processing only (no parallel NSAA parameter generation)
- Gaussian dependency limits cross-platform compatibility

**Recommendations for Advanced Users**:
- Consider RESP charge calculations for improved accuracy
- Validate parameters against experimental data when available
- Perform geometry optimization prior to parameter generation
- Use appropriate solvation models for charge calculations

---

## Complete Workflow

### Standard Workflow

**Step 1: Input preparation**
```bash
python core/prepare_file.py \
    -auto \
    -p protein.pdb \
    --cs-align \
    -d \
    -l ligands.sdf \
    -o ./my_project
```

**Step 2: MD simulation (using CAR method)**
```bash
python core/run_md.py \
    -ln 50 \
    -p ./my_project/prepare \
    -auto
```

**Step 3: Results analysis**
```bash
python core/analyze_result.py \
    -auto \
    -w ./my_project/run
```

### Advanced Workflow Options

**Using custom ligand pairs**:
```bash
# Create pairs.lst file, one ligand pair per line
echo "ligand1-ligand2" > pairs.lst
echo "ligand2-ligand3" >> pairs.lst

python core/prepare_file.py \
    -auto \
    -p protein.pdb \
    --cs-align \
    -d \
    -l ligands.sdf \
    -pl pairs.lst \
    -sp \
    -o ./custom_pairs_project
```

**Handling non-standard amino acids**:
```bash
python core/prepare_file.py \
    -auto \
    -p protein_with_nsaa.pdb \
    --cs-align \
    -d \
    -l ligands.sdf \
    -npd ./nsaa_parameters \
    -o ./nsaa_project
```

**Batch processing mode (preparation only, no auto-execution)**:
```bash
# Prepare files without auto-execution
python core/prepare_file.py \
    -p protein.pdb \
    --cs-align \
    -d \
    -l ligands.sdf \
    -o ./batch_project

# Prepare MD without auto-execution
python core/run_md.py \
    -ln 50 \
    -p ./batch_project/prepare

# Analysis without auto-calculation
python core/analyze_result.py \
    -w ./batch_project/run
```

### File Organization Structure

Complete project file structure:
```
my_project/
├── preprocess/              # Preprocessing files
│   ├── all_generate_mol.sdf # SDF file of all generated molecules
│   ├── dock.pdb            # Protein structure for molecular docking
│   ├── gen_file/           # Detailed ligand pair generation files
│   ├── gen_mol/            # Generated molecular conformation files
│   ├── given_mol/          # Reference ligand mol files
│   ├── pair/               # Ligand pair PDB files
│   ├── protein.pdb         # Preprocessed protein
│   └── rec_prepare/        # Protein preprocessing files
├── prepare/                # Preparation working directory
│   ├── ligand1-ligand2/   # Ligand pair directory
│   │   ├── ligand1.pdb
│   │   ├── ligand2.pdb
│   │   └── submit.sh
│   └── ...
├── run/                   # MD working directory
│   ├── ligand1-ligand2/   # Ligand pair calculation directory
│   │   ├── lM2A/         # Alchemical transformation M→A (solution)
│   │   ├── lM2B/         # Alchemical transformation M→B (solution)
│   │   ├── cM2A/         # Alchemical transformation M→A (complex)
│   │   └── cM2B/         # Alchemical transformation M→B (complex)
│   └── ...
├── calculate/             # Calculation analysis directory
├── results/               # Final results
│   ├── results.csv
│   ├── ene_detail.csv
│   ├── summary.pdf
│   └── summary_abstract.pdf
├── prepare_queue.lst      # Preparation task queue
├── run_all_queue.lst     # MD task queue
├── calculate_queue.lst   # Calculation task queue
└── info.pickle          # Metadata file
```

## FAQ

### Q1: How to choose appropriate Lambda numbers?
**A**: Recommended `-ln 50` (51 Lambda points) with Lambda interval 0.02. Complex transformations can increase to 100, simple ones can reduce to 25.

### Q2: How is the CAR method used?
**A**: CAR (Convergence-Adaptive Roundtrip Method) is Alchemd's core technology and is now enabled by default. It provides excellent sampling efficiency and convergence, ensuring FEP calculation reliability without requiring any manual parameters.

### Q3: How to handle failed ligand pairs?
**A**: Check `err_pairs.lst` file and logs (`run.log`, `calculate.log`). Common issues include molecular structure defects, missing hydrogens, or insufficient computational resources.

### Q4: How to adjust parallel performance?
**A**: Modify the `MAX_WORKER_COUNT` parameter in the `analyze_result.py` script to match computational resources. This parameter controls the maximum number of worker processes in the process pool.

### Q5: How to interpret results?
**A**: Focus on ΔΔG values and uncertainties in `results.csv`. Uncertainty < 1.0 kcal/mol is typically acceptable. Detailed analysis available in `summary.pdf`.

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