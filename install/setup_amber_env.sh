#!/bin/bash

# Alchemd Amber Environment Setup Script
# This script automatically adds Amber initialization to conda/mamba environment

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

# Function to detect environment manager
detect_env_manager() {
    if command -v mamba &> /dev/null; then
        echo "mamba"
    elif command -v conda &> /dev/null; then
        echo "conda"
    else
        echo "none"
    fi
}

# Function to get current environment path
get_current_env_path() {
    local env_manager="$1"
    
    # Try to get from environment variable first
    if [[ -n "$CONDA_PREFIX" ]]; then
        echo "$CONDA_PREFIX"
        return 0
    fi
    
    # Fallback: try to detect from manager commands
    case "$env_manager" in
        "mamba")
            if command -v mamba &> /dev/null; then
                mamba info --base 2>/dev/null || echo ""
            fi
            ;;
        "conda")
            if command -v conda &> /dev/null; then
                conda info --base 2>/dev/null || echo ""
            fi
            ;;
        *)
            echo ""
            ;;
    esac
}

# Function to get current environment name
get_current_env_name() {
    local env_manager="$1"
    
    # Try environment variable first
    if [[ -n "$CONDA_DEFAULT_ENV" ]]; then
        echo "$CONDA_DEFAULT_ENV"
        return 0
    fi
    
    # Try to extract from CONDA_PREFIX
    if [[ -n "$CONDA_PREFIX" ]]; then
        basename "$CONDA_PREFIX"
        return 0
    fi
    
    echo "base"
}

# Function to setup environment activation script
setup_env_activation() {
    local pmemd_path="$1"
    local env_manager="$2"
    local env_path="$3"
    local env_name="$4"
    
    print_status "Setting up PMEMD environment activation for $env_manager environment: $env_name"
    
    # Create activate.d directory if it doesn't exist
    local activate_dir="$env_path/etc/conda/activate.d"
    local deactivate_dir="$env_path/etc/conda/deactivate.d"
    
    if [[ ! -d "$activate_dir" ]]; then
        print_status "Creating activation directory: $activate_dir"
        mkdir -p "$activate_dir"
    fi
    
    if [[ ! -d "$deactivate_dir" ]]; then
        print_status "Creating deactivation directory: $deactivate_dir"
        mkdir -p "$deactivate_dir"
    fi
    
    # Create activation script
    local activate_script="$activate_dir/amber_init.sh"
    local deactivate_script="$deactivate_dir/amber_init.sh"
    
    print_status "Creating Amber25 activation script: $activate_script"
    
    cat > "$activate_script" << EOF
#!/bin/bash
# Auto-generated PMEMD initialization script for Alchemd
# This script is automatically sourced when the conda/mamba environment is activated

# Set AMBERHOME environment variable to PMEMD path
export AMBERHOME="$pmemd_path"

# Source PMEMD environment script if it exists
if [[ -f "$pmemd_path/amber.sh" ]]; then
    source "$pmemd_path/amber.sh"
    echo -e "\\033[0;32m[Alchemd]\\033[0m PMEMD environment activated: $pmemd_path"
else
    echo -e "\\033[0;33m[Alchemd WARNING]\\033[0m PMEMD environment script not found: $pmemd_path/amber.sh"
fi

echo -e "\\033[0;32m[Alchemd]\\033[0m PMEMD environment ready for Alchemd"
EOF

    # Make activation script executable
    chmod +x "$activate_script"
    
    # Create deactivation script
    print_status "Creating Amber25 deactivation script: $deactivate_script"
    
    cat > "$deactivate_script" << EOF
#!/bin/bash
# Auto-generated PMEMD deactivation script for Alchemd
# This script is automatically sourced when the conda/mamba environment is deactivated

# Unset AMBERHOME if it was set by our activation script
if [[ "\$AMBERHOME" == "$pmemd_path" ]]; then
    unset AMBERHOME
    echo -e "\\033[0;34m[Alchemd]\\033[0m PMEMD environment deactivated"
fi
EOF

    # Make deactivation script executable
    chmod +x "$deactivate_script"
    
    print_success "PMEMD environment scripts created successfully"
    print_status "Activation script: $activate_script"
    print_status "Deactivation script: $deactivate_script"
    
    return 0
}

# Function to verify setup
verify_setup() {
    local env_name="$1"
    local pmemd_path="$2"
    
    print_status "Verifying PMEMD environment setup..."
    
    # Check if amber.sh exists for PMEMD
    if [[ -f "$pmemd_path/amber.sh" ]]; then
        print_success "PMEMD environment script exists: $pmemd_path/amber.sh"
    else
        print_warning "PMEMD environment script not found: $pmemd_path/amber.sh"
        print_warning "Make sure PMEMD compilation completed successfully"
    fi
    
    print_status "Environment setup verification completed"
    print_status ""
    print_status "To test the setup:"
    print_status "1. Deactivate and reactivate your $env_name environment:"
    print_status "   conda deactivate && conda activate $env_name"
    print_status "   # or for mamba:"
    print_status "   mamba deactivate && mamba activate $env_name"
    print_status ""
    print_status "2. Check if AMBERHOME is set correctly:"
    print_status "   echo \$AMBERHOME"
    print_status ""
    print_status "3. Test PMEMD tools:"
    print_status "   which pmemd"
    print_status "   which pmemd.cuda"
    print_status "   which pmemd.cuda_SPFP"
    print_status "   pmemd.cuda_SPFP --help"
}

# Main function
main() {
    if [[ $# -ne 1 ]]; then
        print_error "Usage: $0 <PMEMD_PATH>"
        print_error ""
        print_error "  <PMEMD_PATH>       Path to the compiled PMEMD installation"
        print_error ""
        print_error "Example:"
        print_error "  $0 /opt/pmemd24"
        exit 1
    fi
    
    local pmemd_path="$1"
    
    # Convert to absolute paths
    pmemd_path=$(readlink -f "$pmemd_path")
    
    print_status "Setting up PMEMD environment integration for Alchemd"
    print_status "PMEMD path: $pmemd_path"
    
    # Verify directory exists
    if [[ ! -d "$pmemd_path" ]]; then
        print_error "PMEMD directory does not exist: $pmemd_path"
        exit 1
    fi
    
    # Detect environment manager
    local env_manager
    env_manager=$(detect_env_manager)
    
    if [[ "$env_manager" == "none" ]]; then
        print_error "Neither conda nor mamba was found in PATH"
        print_error "This script requires conda or mamba to manage environments"
        exit 1
    fi
    
    print_success "Detected environment manager: $env_manager"
    
    # Get current environment information
    local env_path
    env_path=$(get_current_env_path "$env_manager")
    
    if [[ -z "$env_path" || ! -d "$env_path" ]]; then
        print_error "Could not detect current conda/mamba environment"
        print_error "Make sure you have activated a conda/mamba environment"
        print_error "Current CONDA_PREFIX: ${CONDA_PREFIX:-'not set'}"
        exit 1
    fi
    
    local env_name
    env_name=$(get_current_env_name "$env_manager")
    
    print_success "Current environment: $env_name"
    print_success "Environment path: $env_path"
    
    # Setup environment activation
    if setup_env_activation "$pmemd_path" "$env_manager" "$env_path" "$env_name"; then
        print_success "PMEMD environment integration completed successfully!"
        verify_setup "$env_name" "$pmemd_path"
    else
        print_error "Failed to setup environment activation"
        exit 1
    fi
}

# Run main function with all arguments
main "$@"