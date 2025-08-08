#!/bin/bash

# Alchemd Conda/Mamba Environment Setup Script
# This script sets up Python environment and installs dependencies

set -e  # Exit on any error

# Color codes for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

# Function to print colored output
print_status() {
    echo -e "${BLUE}[INFO]${NC} $1"
}

print_success() {
    echo -e "${GREEN}[SUCCESS]${NC} $1"
}

print_warning() {
    echo -e "${YELLOW}[WARNING]${NC} $1"
}

print_error() {
    echo -e "${RED}[ERROR]${NC} $1"
}

# Function to check if command exists
command_exists() {
    command -v "$1" >/dev/null 2>&1
}

# Get script directory (for relative paths)
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(dirname "$SCRIPT_DIR")"

# Parse command line arguments
SKIP_ENV=false
HELP=false

while [[ $# -gt 0 ]]; do
    case $1 in
        --skip-env)
            SKIP_ENV=true
            shift
            ;;
        -h|--help)
            HELP=true
            shift
            ;;
        *)
            print_error "Unknown option: $1"
            print_error "Use --help for usage information"
            exit 1
            ;;
    esac
done

# Show help if requested
if [[ "$HELP" == true ]]; then
    echo "Usage: $0 [OPTIONS]"
    echo ""
    echo "Options:"
    echo "  --skip-env    Skip environment creation, only install dependencies"
    echo "                (assumes alchemd environment already exists)"
    echo "  -h, --help    Show this help message"
    echo ""
    echo "Example:"
    echo "  $0              # Full setup including environment creation"
    echo "  $0 --skip-env   # Skip environment creation, only install dependencies"
    exit 0
fi

print_status "Starting Alchemd Python environment setup..."
print_status "Project root: $PROJECT_ROOT"

if [[ "$SKIP_ENV" == true ]]; then
    print_status "Skip environment creation mode enabled"
fi

# Detect conda or mamba (prefer mamba)
CONDA_CMD=""
if command_exists mamba; then
    CONDA_CMD="mamba"
    print_status "Found mamba - using mamba for faster package management"
elif command_exists conda; then
    CONDA_CMD="conda"
    print_status "Found conda - using conda for package management"
else
    print_error "Neither conda nor mamba found in PATH"
    print_error "Please install Miniconda/Anaconda or Mamba first"
    print_error "Visit: https://docs.conda.io/en/latest/miniconda.html"
    exit 1
fi

# Check conda configuration
print_status "Checking conda configuration..."
if $CONDA_CMD config --show channels | grep -q "conda-forge"; then
    print_success "conda-forge channel is configured"
else
    print_warning "conda-forge channel not found in conda configuration"
    print_warning "Some packages may not be available. Consider adding conda-forge:"
    print_warning "  conda config --add channels conda-forge"
fi

# Environment name
ENV_NAME="alchemd"
PYTHON_VERSION="3.9.21"

# Main execution logic
if [[ "$SKIP_ENV" == true ]]; then
    # Skip environment setup mode - jump directly to external tools installation
    print_status "Skipping environment setup, jumping to external tools installation..."
    print_status "Note: This assumes '$ENV_NAME' environment already exists and is properly configured"
else
    # Full environment setup mode
    print_status "Starting full environment setup..."
    
    # Check if environment already exists
    if $CONDA_CMD env list | grep -q "^$ENV_NAME "; then
        print_warning "Environment '$ENV_NAME' already exists"
        print_status "Current environment information:"
        $CONDA_CMD env list | grep "^$ENV_NAME "
        print_status ""
        print_status "Choose an option:"
        print_status "  [N] Use existing environment (default, recommended)"
        print_status "  [y] Remove and recreate environment"
        read -p "Do you want to remove and recreate it? (y/N): " -r
        if [[ $REPLY =~ ^[Yy]$ ]]; then
            print_status "Removing existing environment..."
            $CONDA_CMD env remove -n $ENV_NAME -y
            print_success "Environment removed"
        else
            print_success "Using existing environment (packages will be updated if needed)"
        fi
    fi

    # Create conda environment
    print_status "Creating conda environment '$ENV_NAME' with Python $PYTHON_VERSION..."
    if ! $CONDA_CMD create -n $ENV_NAME python=$PYTHON_VERSION -y; then
        print_error "Failed to create conda environment"
        print_warning "Trying with available Python 3.9 version..."
        if ! $CONDA_CMD create -n $ENV_NAME python=3.9 -y; then
            print_error "Failed to create conda environment with Python 3.9"
            exit 1
        fi
    fi
    print_success "Conda environment created successfully"

    # Activate environment
    print_status "Activating environment..."
    source "$($CONDA_CMD info --base)/etc/profile.d/conda.sh"
    $CONDA_CMD activate $ENV_NAME

    # Verify Python version
    PYTHON_VER=$(python --version 2>&1)
    print_status "Using Python: $PYTHON_VER"

    # Install dependencies using conda for unified package management
    REQUIREMENTS_FILE="$PROJECT_ROOT/requirements.txt"
    if [ -f "$REQUIREMENTS_FILE" ]; then
        print_status "Installing Python dependencies using conda for unified management..."
        
        # First, install packages from requirements.txt using conda with user's configured channels
        print_status "Installing packages from requirements.txt via conda..."
        if $CONDA_CMD install -n $ENV_NAME --file "$REQUIREMENTS_FILE" -y; then
            print_success "Requirements installed via conda"
        else
            print_warning "Some packages from requirements.txt failed via conda, trying individual installation..."
            
            # Read requirements and install individually
            while IFS= read -r package || [ -n "$package" ]; do
                # Skip empty lines and comments
                if [[ -z "$package" || "$package" =~ ^#.* ]]; then
                    continue
                fi
                
                print_status "Installing $package..."
                if ! $CONDA_CMD install -n $ENV_NAME "$package" -y; then
                    print_warning "Failed to install $package via conda"
                fi
            done < "$REQUIREMENTS_FILE"
        fi
        
        # Install specific version-constrained packages to override defaults
        print_status "Enforcing specific version constraints..."
        VERSIONED_PACKAGES=(
            "numpy<=1.26.0"
            "scipy<=1.12.0"
            "pandas>=2.2.3"
            "rdkit>=2024.03.5"
            "pymbar<=3.1.1"
            "pillow<=10.4.0"
            "pypdf2>3.0.0"
        )
        
        # Install all versioned packages in a single command to avoid conflicts
        print_status "Installing version-constrained packages in a single operation..."
        if $CONDA_CMD install -n $ENV_NAME "${VERSIONED_PACKAGES[@]}" -y --force-reinstall; then
            print_success "All version constraints applied successfully"
        else
            print_warning "Failed to install some version-constrained packages in batch mode"
            print_status "Falling back to individual installation..."
            
            # Fallback: install individually
            for package in "${VERSIONED_PACKAGES[@]}"; do
                print_status "Installing $package with version constraint..."
                if $CONDA_CMD install -n $ENV_NAME "$package" -y --force-reinstall; then
                    print_success "$package version constraint applied"
                else
                    print_warning "Failed to apply version constraint for $package"
                fi
            done
        fi
        
        print_success "Main dependencies installed"
    else
        print_error "requirements.txt not found at $REQUIREMENTS_FILE"
        exit 1
    fi

    # Install critical packages that often have import issues
    print_status "Installing critical packages with special handling..."
    CRITICAL_PACKAGES=("parmed" "openmm")

    for package in "${CRITICAL_PACKAGES[@]}"; do
        print_status "Installing $package via conda using user's configured channels..."
        if $CONDA_CMD install -n $ENV_NAME "$package" -y --force-reinstall; then
            print_success "$package installed successfully via conda"
        else
            print_error "Failed to install $package via conda"
            print_warning "This may cause import issues later"
            print_warning "Please ensure your conda channels include conda-forge, bioconda, or omnia"
        fi
    done
    
    print_status "Environment setup completed, proceeding to external tools..."
fi

# Create dependencies directory
DEPS_DIR="$PROJECT_ROOT/dependencies"
mkdir -p "$DEPS_DIR"

# Extract tar archives from dependencies directory
print_status "Extracting dependency tar archives..."
cd "$DEPS_DIR"

TAR_FILES=("AlchemConvTools.tar" "CAR-FEP.tar" "genambRBFE.tar")
for tar_file in "${TAR_FILES[@]}"; do
    if [ -f "$tar_file" ]; then
        print_status "Extracting $tar_file..."
        if tar -xvf "$tar_file"; then
            print_success "$tar_file extracted successfully"
        else
            print_error "Failed to extract $tar_file"
        fi
    else
        print_warning "$tar_file not found, skipping extraction"
    fi
done

cd "$PROJECT_ROOT"

# Download and setup WatVina
print_status "Setting up WatVina..."
WATVINA_DIR="$DEPS_DIR/watvina"
mkdir -p "$WATVINA_DIR"

WATVINA_URL="https://github.com/biocheming/watvina/releases/download/v20241125/watvina"
WATVINA_PATH="$WATVINA_DIR/watvina"

if [ ! -f "$WATVINA_PATH" ]; then
    print_status "Downloading WatVina..."
    if wget -O "$WATVINA_PATH" "$WATVINA_URL"; then
        print_success "WatVina downloaded successfully"
    else
        print_error "Failed to download WatVina"
        print_warning "Please download manually from: $WATVINA_URL"
    fi
fi

# Always ensure WatVina has correct permissions
if [ -f "$WATVINA_PATH" ]; then
    print_status "Setting executable permissions for WatVina..."
    chmod +x "$WATVINA_PATH"
    if [ -x "$WATVINA_PATH" ]; then
        print_success "WatVina is now executable"
    else
        print_error "Failed to set executable permissions for WatVina"
    fi
else
    print_warning "WatVina file not found, skipping permissions setup"
fi

# Download Weighted_cc
print_status "Setting up Weighted_cc..."
WEIGHTED_CC_ZIP="$DEPS_DIR/weighted_cc.zip"
WEIGHTED_CC_URL="https://github.com/zlisysu/Weighted_cc/archive/refs/heads/main.zip"

if [ ! -d "$DEPS_DIR/Weighted_cc-main" ]; then
    print_status "Downloading Weighted_cc..."
    if wget -O "$WEIGHTED_CC_ZIP" "$WEIGHTED_CC_URL"; then
        cd "$DEPS_DIR"
        if unzip "$WEIGHTED_CC_ZIP"; then
            rm -f "$WEIGHTED_CC_ZIP"
            print_success "Weighted_cc downloaded and extracted"
        else
            print_error "Failed to extract Weighted_cc"
        fi
        cd "$PROJECT_ROOT"
    else
        print_error "Failed to download Weighted_cc"
        print_warning "Please download manually from: $WEIGHTED_CC_URL"
    fi
else
    print_success "Weighted_cc already exists"
fi

# Download and install Pyautomd
print_status "Setting up Pyautomd..."
PYAUTOMD_ZIP="$DEPS_DIR/pyautomd.zip"
PYAUTOMD_URL="https://github.com/phloglucinol/Pyautomd/archive/refs/heads/main.zip"

if [ ! -d "$DEPS_DIR/Pyautomd-main" ]; then
    print_status "Downloading Pyautomd..."
    if wget -O "$PYAUTOMD_ZIP" "$PYAUTOMD_URL"; then
        cd "$DEPS_DIR"
        if unzip "$PYAUTOMD_ZIP"; then
            rm -f "$PYAUTOMD_ZIP"
            print_success "Pyautomd downloaded and extracted"
        else
            print_error "Failed to extract Pyautomd"
        fi
        cd "$PROJECT_ROOT"
    else
        print_error "Failed to download Pyautomd"
        print_warning "Please download manually from: $PYAUTOMD_URL"
    fi
fi

# Install Pyautomd package
PYAUTOMD_PATH="$DEPS_DIR/Pyautomd-main"
if [ -d "$PYAUTOMD_PATH" ]; then
    print_status "Installing Pyautomd package..."
    if pip install "$PYAUTOMD_PATH/"; then
        print_success "Pyautomd package installed successfully"
    else
        print_error "Failed to install Pyautomd package"
        print_warning "You may need to install it manually: pip install $PYAUTOMD_PATH/"
    fi
else
    print_warning "Pyautomd directory not found, skipping installation"
fi

# Install Lomap via conda
print_status "Setting up Lomap2 via conda..."
if [[ "$SKIP_ENV" != true ]]; then
    # Only install if we're not skipping environment setup
    if $CONDA_CMD install -n $ENV_NAME conda-forge::lomap2 -y; then
        print_success "Lomap2 installed successfully via conda"
    else
        print_error "Failed to install lomap2 via conda"
        print_warning "Please ensure conda-forge channel is available"
        print_warning "You can try manually: conda install -c conda-forge lomap2"
    fi
else
    print_status "Skipping lomap2 installation (environment setup skipped)"
    print_warning "Please ensure lomap2 is installed in your environment manually if needed"
fi

# Verify installation
print_status "Verifying installation..."

# Test Python imports
print_status "Testing critical Python packages..."
CRITICAL_PACKAGES=("numpy" "pandas" "rdkit" "matplotlib" "parmed" "openmm")
FAILED_IMPORTS=()

for package in "${CRITICAL_PACKAGES[@]}"; do
    print_status "Testing import: $package"
    if python -c "import $package; print(f'✓ {$package.__version__ if hasattr($package, \"__version__\") else \"imported successfully\"}')" 2>/dev/null; then
        print_success "$package import successful"
    else
        FAILED_IMPORTS+=("$package")
        print_warning "$package import failed"
    fi
done

# Special test for problematic packages
print_status "Running additional tests for critical packages..."

# Test parmed specifically
if python -c "import parmed; print('parmed version:', parmed.__version__)" 2>/dev/null; then
    print_success "parmed detailed test passed"
else
    print_error "parmed detailed test failed - this may cause issues"
fi

# Test openmm specifically  
if python -c "import openmm; print('OpenMM version:', openmm.__version__)" 2>/dev/null; then
    print_success "OpenMM detailed test passed"
else
    print_error "OpenMM detailed test failed - this may cause issues"
fi

# Test lomap specifically (conda installed version)
if python -c "import lomap; print('Lomap version:', getattr(lomap, '__version__', 'version not available'))" 2>/dev/null; then
    print_success "Lomap detailed test passed"
else
    print_error "Lomap detailed test failed - this may cause issues"
    print_warning "If using --skip-env, ensure lomap2 is installed: conda install -c conda-forge lomap2"
fi


if [ ${#FAILED_IMPORTS[@]} -gt 0 ]; then
    print_warning "Some packages failed to import: ${FAILED_IMPORTS[*]}"
    print_warning "Please check these packages manually"
    print_warning "Critical packages (parmed, openmm, lomap) are essential for Alchemd functionality"
fi

# Check external tools
print_status "Checking external tools..."
if [ -x "$WATVINA_PATH" ]; then
    print_success "WatVina is executable"
else
    print_warning "WatVina is not executable or not found"
fi

print_success "Python environment setup completed!"
print_status ""
print_status "IMPORTANT: This setup is for Alchemd - FEP calculations toolkit"
print_status ""
print_status "To activate the environment, run:"
print_status "  conda activate $ENV_NAME"
print_status ""
print_status "Usage options for this script:"
print_status "  ./install/setup_conda.sh             # Full setup (default)"
print_status "  ./install/setup_conda.sh --skip-env  # Skip environment creation,"
print_status "                                       # only install dependencies"
print_status ""
print_status "To test the installation, try:"
print_status "  conda activate $ENV_NAME"
print_status "  python core/prepare_file.py --help"
print_status ""
print_status "Next step: Build PMEMD (REQUIRED for Alchemd):"
print_status "  chmod +x install/build_amber.sh"
print_status "  ./install/build_amber.sh /path/to/pmemd24"
print_status ""
print_warning "Note: Amber license is required. Apply at https://ambermd.org/"
print_warning "Important: Only PMEMD with CUDA support is required"

# Initialize Alchemd configuration
print_status ""
print_status "Initializing Alchemd configuration..."
cd "$PROJECT_ROOT"

# Run configuration initialization in auto mode
if python core/init_setup.py --auto; then
    print_success "Configuration initialized successfully"
    print_status "Configuration file created: core/configs.toml"
    print_status ""
    print_status "You can reconfigure anytime by running:"
    print_status "  python core/init_setup.py --interactive"
else
    print_warning "Configuration initialization failed"
    print_warning "You can run it manually later:"
    print_warning "  conda activate $ENV_NAME"
    print_warning "  python core/init_setup.py"
fi

print_status ""
print_success "Python environment and configuration setup completed!"
print_status ""
print_status "Environment: $ENV_NAME"
print_status "Project root: $PROJECT_ROOT"
print_status "Configuration: core/configs.toml"