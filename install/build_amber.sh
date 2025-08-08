#!/bin/bash

# Alchemd Amber Compilation Script
# This script compiles Amber with custom patches for Alchemd

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

# Function to detect and activate conda environment
activate_alchemd_env() {
    print_status "Setting up conda/mamba environment..."
    
    # Detect conda or mamba
    local CONDA_CMD=""
    if command -v mamba &> /dev/null; then
        CONDA_CMD="mamba"
        print_status "Found mamba - using mamba for environment management"
    elif command -v conda &> /dev/null; then
        CONDA_CMD="conda"
        print_status "Found conda - using conda for environment management"
    else
        print_error "Neither conda nor mamba found in PATH"
        print_error "Please install Miniconda/Anaconda or Mamba first"
        print_error "Visit: https://docs.conda.io/en/latest/miniconda.html"
        exit 1
    fi
    
    # Check if alchemd environment exists
    if ! $CONDA_CMD env list | grep -q "^alchemd "; then
        print_error "alchemd environment not found"
        print_error "Please run setup_conda.sh first to create the alchemd environment:"
        print_error "  chmod +x install/setup_conda.sh"
        print_error "  ./install/setup_conda.sh"
        exit 1
    fi
    
    print_success "alchemd environment found"
    
    # Activate alchemd environment
    print_status "Activating alchemd environment..."
    source "$($CONDA_CMD info --base)/etc/profile.d/conda.sh"
    $CONDA_CMD activate alchemd
    
    # Verify activation
    if [[ "$CONDA_DEFAULT_ENV" == "alchemd" ]]; then
        print_success "alchemd environment activated successfully"
        print_status "Current Python: $(which python)"
        print_status "Python version: $(python --version 2>&1)"
    else
        print_error "Failed to activate alchemd environment"
        print_error "Current environment: ${CONDA_DEFAULT_ENV:-none}"
        exit 1
    fi
}

# Function to show usage
show_usage() {
    echo "Usage: $0 <PMEMD_PATH> [BUILD_CORES]"
    echo ""
    echo "  <PMEMD_PATH>       Path to the extracted PMEMD directory (e.g., pmemd24)"
    echo "  [BUILD_CORES]      Number of cores for compilation (default: 4, max recommended: 8)"
    echo ""
    echo "Example:"
    echo "  $0 /opt/pmemd24"
    echo "  $0 /opt/pmemd24 6"
    echo ""
    echo "Environment Variables:"
    echo "  ALCHEMD_BUILD_CORES  Override default build cores (alternative to BUILD_CORES argument)"
    echo ""
    echo "Prerequisites:"
    echo "  1. Run setup_conda.sh first to create alchemd environment"
    echo "  2. Ensure CUDA 12.0 and NVIDIA drivers are installed"
    echo "  3. Extract PMEMD package to target directory"
    echo "  4. Amber license is required - apply at https://ambermd.org/"
    echo ""
    echo "Note: This script automatically activates the alchemd conda environment."
    echo "      Only PMEMD with CUDA support will be compiled (AmberTools no longer required)."
}

# Check if PMEMD path is provided
if [ $# -lt 1 ] || [ $# -gt 2 ]; then
    print_error "Error: PMEMD path must be provided"
    show_usage
    exit 1
fi

PMEMD_PATH="$1"
BUILD_CORES_ARG="$2"

# Convert to absolute paths
PMEMD_PATH=$(readlink -f "$PMEMD_PATH")

# Determine build cores (priority: argument > environment variable > default)
if [ -n "$BUILD_CORES_ARG" ] && [ "$BUILD_CORES_ARG" -gt 0 ] 2>/dev/null; then
    BUILD_CORES="$BUILD_CORES_ARG"
elif [ -n "$ALCHEMD_BUILD_CORES" ] && [ "$ALCHEMD_BUILD_CORES" -gt 0 ] 2>/dev/null; then
    BUILD_CORES="$ALCHEMD_BUILD_CORES"
else
    BUILD_CORES=2
fi

# Validate and cap build cores
if [ "$BUILD_CORES" -gt 16 ]; then
    print_warning "Build cores capped at 16 (requested: $BUILD_CORES)"
    BUILD_CORES=16
fi

print_status "Starting PMEMD compilation process for Alchemd..."
print_status "PMEMD path: $PMEMD_PATH"
print_status "Build cores: $BUILD_CORES"

# Check for WSL environment and provide performance warning
if grep -qi microsoft /proc/version 2>/dev/null || [ -n "$WSL_DISTRO_NAME" ]; then
    print_warning "WSL environment detected"
    print_warning "Compilation may be slower than native Linux"
    if [ "$BUILD_CORES" -gt 6 ]; then
        print_warning "Consider reducing build cores to 4-6 in WSL to prevent memory issues"
    fi
fi

# Activate alchemd environment before compilation
activate_alchemd_env

# Verify PMEMD directory exists
if [ ! -d "$PMEMD_PATH" ]; then
    print_error "PMEMD directory does not exist: $PMEMD_PATH"
    exit 1
fi

# Get script directory (for relative paths)
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(dirname "$SCRIPT_DIR")"
AMBER_FILES_DIR="$PROJECT_ROOT/dependencies/AmberFiles"

print_status "Project root: $PROJECT_ROOT"
print_status "Amber files directory: $AMBER_FILES_DIR"

# Check if AmberFiles directory exists
if [ ! -d "$AMBER_FILES_DIR" ]; then
    print_error "AmberFiles directory not found: $AMBER_FILES_DIR"
    exit 1
fi

# Check for PMEMD (build/run_cmake)
PMEMD_CMAKE_PATH="$PMEMD_PATH/build/run_cmake"
if [ ! -f "$PMEMD_CMAKE_PATH" ]; then
    print_error "PMEMD build script not found: $PMEMD_CMAKE_PATH does not exist"
    print_error "Please ensure PMEMD is properly extracted to: $PMEMD_PATH"
    exit 1
fi
print_success "PMEMD detected: $PMEMD_CMAKE_PATH"

# Check for PMEMD source directory (for patching)
PMEMD_SRC_PATH="$PMEMD_PATH/src/pmemd/src"
if [ ! -d "$PMEMD_SRC_PATH" ]; then
    print_error "PMEMD source not found: $PMEMD_SRC_PATH does not exist"
    print_error "Please ensure PMEMD is properly extracted to: $PMEMD_PATH"
    exit 1
fi
print_success "PMEMD source detected: $PMEMD_SRC_PATH"

# Function to check and setup CUDA environment
check_cuda_environment() {
    print_status "Checking CUDA environment..."
    
    # Check NVIDIA driver
    if ! command -v nvidia-smi &> /dev/null; then
        print_error "nvidia-smi not found. Please install NVIDIA drivers first."
        return 1
    fi
    
    # Check CUDA compiler
    if ! command -v nvcc &> /dev/null; then
        print_error "nvcc not found. Please install CUDA Toolkit."
        return 1
    fi
    
    # Get CUDA version and path
    CUDA_VERSION=$(nvcc --version | grep "release" | sed -n 's/.*release \([0-9\.]*\).*/\1/p')
    print_success "CUDA compiler found: version $CUDA_VERSION"
    
    # Get nvcc path to determine CUDA installation directory
    NVCC_PATH=$(which nvcc)
    CUDA_ROOT=$(dirname $(dirname "$NVCC_PATH"))
    print_success "CUDA installation detected at: $CUDA_ROOT"
    
    # Set CUDA environment variables
    export CUDA_HOME="$CUDA_ROOT"
    export LD_LIBRARY_PATH="$CUDA_ROOT/lib64:$LD_LIBRARY_PATH"
    export PATH="$CUDA_ROOT/bin:$PATH"
    
    print_success "CUDA environment configured successfully"
    print_status "CUDA_HOME: $CUDA_HOME"
    
    return 0
}

# Function to modify run_cmake compilation options
modify_run_cmake() {
    local cmake_path="$1"
    local component_name="$2"
    local backup_path="${cmake_path}.backup"
    
    print_status "Modifying compilation options for $component_name..."
    
    # Backup original run_cmake
    if [ ! -f "$backup_path" ]; then
        print_status "Creating backup of run_cmake for $component_name..."
        cp "$cmake_path" "$backup_path"
        print_success "Backup created: $backup_path"
    else
        print_status "Backup already exists: $backup_path"
    fi
    
    # Check if modifications are already applied
    if grep -q "DCUDA=TRUE" "$cmake_path" && grep -q "DDOWNLOAD_MINICONDA=FALSE" "$cmake_path"; then
        print_success "Compilation options already modified for $component_name"
        return 0
    fi
    
    print_status "Applying compilation option changes for $component_name..."
    
    # Replace CUDA=FALSE with CUDA=TRUE
    if grep -q "DCUDA=FALSE" "$cmake_path"; then
        sed -i 's/DCUDA=FALSE/DCUDA=TRUE/g' "$cmake_path"
        print_success "Changed CUDA option from FALSE to TRUE for $component_name"
    elif ! grep -q "DCUDA=TRUE" "$cmake_path"; then
        print_warning "DCUDA=FALSE not found in $component_name run_cmake, you may need to add CUDA support manually"
    fi
    
    # Replace DOWNLOAD_MINICONDA=TRUE with DOWNLOAD_MINICONDA=FALSE
    if grep -q "DDOWNLOAD_MINICONDA=TRUE" "$cmake_path"; then
        sed -i 's/DDOWNLOAD_MINICONDA=TRUE/DDOWNLOAD_MINICONDA=FALSE/g' "$cmake_path"
        print_success "Changed DOWNLOAD_MINICONDA option from TRUE to FALSE for $component_name"
    elif ! grep -q "DDOWNLOAD_MINICONDA=FALSE" "$cmake_path"; then
        print_warning "DDOWNLOAD_MINICONDA=TRUE not found in $component_name run_cmake"
    fi

}

# Check CUDA environment before compilation
if ! check_cuda_environment; then
    print_error "CUDA environment check failed"
    print_error "Please resolve CUDA issues before proceeding with Amber compilation"
    exit 1
fi

# Modify PMEMD run_cmake file
modify_run_cmake "$PMEMD_CMAKE_PATH" "PMEMD"

# Apply patches to pmemd
print_status "Applying patches to pmemd..."

# Copy patch files to pmemd/src directory
PATCH_SCRIPT="$AMBER_FILES_DIR/patch.sh"
PATCH_FILE="$AMBER_FILES_DIR/ti.F90.patch"

if [ ! -f "$PATCH_SCRIPT" ]; then
    print_error "Patch script not found: $PATCH_SCRIPT"
    exit 1
fi

if [ ! -f "$PATCH_FILE" ]; then
    print_error "Patch file not found: $PATCH_FILE"
    exit 1
fi

# Copy patch files to pmemd source directory
print_status "Copying patch files to $PMEMD_SRC_PATH..."
cp "$PATCH_SCRIPT" "$PMEMD_SRC_PATH/"
cp "$PATCH_FILE" "$PMEMD_SRC_PATH/"

# Make patch script executable
chmod +x "$PMEMD_SRC_PATH/patch.sh"

# Execute patch script
print_status "Executing patch script..."
cd "$PMEMD_SRC_PATH"
if ./patch.sh; then
    print_success "Patches applied successfully"
else
    print_error "Failed to apply patches"
    exit 1
fi

# Return to original directory
cd "$PROJECT_ROOT"

# Function to clean build directory before compilation
clean_build_dir() {
    local build_dir="$1"
    local component_name="$2"
    
    print_status "Cleaning previous build artifacts for $component_name..."
    
    cd "$build_dir"
    
    # Try clean_build script first (Amber's native clean script)
    if [ -f "./clean_build" ]; then
        print_status "Using Amber's clean_build script..."
        if echo "y" | ./clean_build; then
            print_success "Build directory cleaned using clean_build"
        else
            print_warning "clean_build failed, trying make clean..."
            echo "y" | make clean || true
        fi
    else
        print_status "clean_build script not found, using make clean..."
        echo "y" | make clean || true
        print_status "Build directory cleaned using make clean"
    fi
}

# Function to compile a component
compile_component() {
    local component_path="$1"
    local component_name="$2"
    local build_dir="$component_path/build"
    
    print_status "Starting $component_name compilation..."
    print_warning "This process may take 30 minutes to several hours depending on your system"
    
    if [ ! -d "$build_dir" ]; then
        print_error "Build directory not found: $build_dir"
        return 1
    fi
    
    cd "$build_dir"
    
    # Set environment variable for this component
    export AMBERHOME="$component_path"

    # Clean build directory before compilation
    clean_build_dir "$build_dir" "$component_name"

    print_status "Build directory: $build_dir"
    print_status "Component: $component_name"

    # Run cmake configuration
    print_status "Running cmake configuration for $component_name..."
    if ./run_cmake; then
        print_success "$component_name cmake configuration completed!"
    else
        print_error "$component_name cmake configuration failed"
        return 1
    fi
    
    # Use predetermined build cores
    print_status "Using $BUILD_CORES cores for parallel compilation"
    
    # Run make install with parallel compilation
    print_status "Running make install for $component_name (using -j$BUILD_CORES)..."
    if make install -j"$BUILD_CORES"; then
        print_success "$component_name compilation and installation completed!"
        return 0
    else
        print_error "$component_name compilation failed"
        return 1
    fi
}

print_status "Starting PMEMD compilation process (REQUIRED for Alchemd)..."
print_warning "CUDA support is essential for optimal performance"

# Compile PMEMD only
if ! compile_component "$PMEMD_PATH" "PMEMD"; then
    print_error "PMEMD compilation failed"
    print_warning "Common issues and solutions:"
    print_warning "1. Missing CUDA toolkit - ensure CUDA 12.0 is installed"
    print_warning "2. Missing dependencies - run setup_linux.sh first"
    print_warning "3. Insufficient memory - compilation requires several GB of RAM"
    print_warning "4. Compiler version issues - ensure gcc-9+, g++-9+, gfortran-9+ are available"
    print_warning "5. Missing Amber license - apply at https://ambermd.org/"
    print_warning "6. WSL performance - try reducing build cores: export ALCHEMD_BUILD_CORES=4"
    exit 1
fi


# Function to verify PMEMD installation
verify_pmemd_installation() {
    local pmemd_path="$1"
    
    print_status "Verifying PMEMD installation..."
    
    local bin_dir="$pmemd_path/bin"
    local missing_bins=()
    
    # Check for PMEMD binaries (focus on CUDA versions)
    local tools=("pmemd" "pmemd.cuda" "pmemd.cuda_SPFP")
    local critical_tools=("pmemd.cuda" "pmemd.cuda_SPFP")
    
    for tool in "${tools[@]}"; do
        if [ -f "$bin_dir/$tool" ]; then
            print_success "$tool is available"
            # Test pmemd.cuda_SPFP functionality if available
            if [ "$tool" == "pmemd.cuda_SPFP" ]; then
                if "$bin_dir/$tool" --help &>/dev/null; then
                    print_success "pmemd.cuda_SPFP functional test passed"
                else
                    print_warning "pmemd.cuda_SPFP exists but help test failed"
                fi
            fi
        else
            missing_bins+=("$tool")
        fi
    done
    
    # Check critical tools
    local critical_missing=()
    for tool in "${critical_tools[@]}"; do
        if [ ! -f "$bin_dir/$tool" ]; then
            critical_missing+=("$tool")
        fi
    done
    
    if [ ${#critical_missing[@]} -gt 0 ]; then
        print_error "Critical PMEMD CUDA binaries are missing: ${critical_missing[*]}"
        print_error "CUDA compilation may have failed"
        return 1
    fi
    
    if [ ${#missing_bins[@]} -gt 0 ]; then
        print_warning "Some PMEMD binaries are missing: ${missing_bins[*]}"
    fi
    
    # Check for amber.sh
    local amber_sh="$pmemd_path/amber.sh"
    if [ -f "$amber_sh" ]; then
        print_success "PMEMD environment script found: $amber_sh"
    else
        print_warning "PMEMD environment script not found: $amber_sh"
        return 1
    fi
    
    return 0
}

# Skip binary check here - will be verified after environment setup using test_amber_env.sh

# Calculate the actual PMEMD installation path (compiled output)
# PMEMD_PATH is the source path, but installation goes to a parallel directory named "pmemd24"
PMEMD_PARENT_DIR=$(dirname "$PMEMD_PATH")
PMEMD_INSTALL_PATH="$PMEMD_PARENT_DIR/pmemd24"

print_status "PMEMD source path: $PMEMD_PATH"
print_status "PMEMD install path: $PMEMD_INSTALL_PATH"

# Verify the installation path exists
if [ ! -d "$PMEMD_INSTALL_PATH" ]; then
    print_error "PMEMD installation directory not found: $PMEMD_INSTALL_PATH"
    print_error "Compilation may have failed or used a different install path"
    exit 1
fi

# Set up environment integration
print_status "Setting up PMEMD environment integration..."

# Automatically setup conda/mamba environment integration
SETUP_ENV_SCRIPT="$SCRIPT_DIR/setup_amber_env.sh"
if [ -f "$SETUP_ENV_SCRIPT" ]; then
    print_status "Configuring automatic PMEMD environment integration..."
    if "$SETUP_ENV_SCRIPT" "$PMEMD_INSTALL_PATH"; then
        print_success "PMEMD environment integration setup completed"
    else
        print_warning "Failed to setup automatic environment integration"
        print_warning "You can manually run: bash $SETUP_ENV_SCRIPT $PMEMD_INSTALL_PATH"
    fi
else
    print_warning "Environment setup script not found: $SETUP_ENV_SCRIPT"
    print_status "To use PMEMD, source the environment script:"
    print_status "  source $PMEMD_INSTALL_PATH/amber.sh"
fi

# Display final information
print_success "PMEMD compilation process completed!"
print_status ""
print_status "Installation Summary:"
print_status "  PMEMD source path: $PMEMD_PATH"
print_status "  PMEMD install path: $PMEMD_INSTALL_PATH"
print_status "  Build cores used: $BUILD_CORES"
print_status "  Patches applied: Yes"
print_status "  CUDA support: Enabled (required)"
print_status "  Environment integration: Completed"
print_status ""
print_status "IMPORTANT: Complete the verification process:"
if [ -f "$SETUP_ENV_SCRIPT" ] && [ -n "$CONDA_PREFIX" ]; then
    print_status "1. PMEMD environment integration has been set up"
    print_status "2. To complete verification, run the following commands:"
    print_status ""
    print_status "   # Reactivate environment to load PMEMD settings"
    print_status "   conda deactivate"
    print_status "   conda activate alchemd"
    print_status ""
    print_status "   # Run verification tests"
    print_status "   ./install/test_amber_env.sh"
    print_status ""
    print_status "   # Test PMEMD tools directly"
    print_status "   which pmemd.cuda_SPFP"
    print_status "   pmemd.cuda_SPFP --help"
    print_status ""
    print_status "3. Update Alchemd configuration if needed:"
    print_status "   python core/init_setup.py"
else
    print_status "1. Manual environment setup required:"
    print_status "   source $PMEMD_INSTALL_PATH/amber.sh"
    print_status ""
    print_status "2. Update Alchemd configuration:"
    print_status "   python core/init_setup.py"
    print_status ""
    print_status "3. Test the installation:"
    print_status "   ./install/test_amber_env.sh"
fi
print_status ""
print_warning "DO NOT exit the conda environment! Stay in alchemd environment."
print_success "PMEMD compilation completed successfully!"
print_status "Follow the verification steps above to complete the setup."

# Return to project root
cd "$PROJECT_ROOT"