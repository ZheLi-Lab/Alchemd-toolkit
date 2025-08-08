#!/bin/bash

# Alchemd Linux Environment Setup Script
# This script installs required system packages for Alchemd

set -e  # Exit on any error

# Ensure we're running with bash
if [ -z "$BASH_VERSION" ]; then
    echo "Error: This script requires bash to run properly"
    echo "Please run with: bash $0"
    exit 1
fi

# Enable debug mode if requested
if [[ "$1" == "--debug" || "$DEBUG" == "1" ]]; then
    set -x
    print_status "Debug mode enabled"
fi

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

print_status "Starting Alchemd Linux environment setup..."

# Detect Linux distribution and package manager
if [ -f /etc/os-release ]; then
    . /etc/os-release
    DISTRO="$ID"
    VERSION="$VERSION_ID"
else
    print_error "Cannot detect Linux distribution"
    exit 1
fi

print_status "Detected Linux distribution: $DISTRO $VERSION"

# Detect package manager
PKG_MANAGER=""
if command_exists apt-get; then
    PKG_MANAGER="apt"
    UPDATE_CMD="sudo apt-get update"
    INSTALL_CMD="sudo apt-get install -y"
elif command_exists dnf; then
    PKG_MANAGER="dnf"
    UPDATE_CMD="sudo dnf update -y"
    INSTALL_CMD="sudo dnf install -y"
elif command_exists yum; then
    PKG_MANAGER="yum"
    UPDATE_CMD="sudo yum update -y"
    INSTALL_CMD="sudo yum install -y"
elif command_exists zypper; then
    PKG_MANAGER="zypper"
    UPDATE_CMD="sudo zypper refresh"
    INSTALL_CMD="sudo zypper install -y"
elif command_exists pacman; then
    PKG_MANAGER="pacman"
    UPDATE_CMD="sudo pacman -Sy"
    INSTALL_CMD="sudo pacman -S --noconfirm"
else
    print_error "No supported package manager found"
    exit 1
fi

print_status "Detected package manager: $PKG_MANAGER"

# Update package lists
print_status "Updating package lists..."
if ! $UPDATE_CMD 2>/dev/null; then
    print_warning "Failed to update package lists, continuing anyway..."
    print_status "This might be due to network issues or repository problems"
else
    print_success "Package update completed"
fi

# Function to install GCC toolchain (gcc9+ compatible)
install_gcc_toolchain() {
    print_status "Installing GCC toolchain (gcc9+ required)..."
    
    case "$PKG_MANAGER" in
        "apt")
            # Ubuntu/Debian
            # First try to install default GCC and check version
            $INSTALL_CMD gcc g++ gfortran
            if command -v gcc &> /dev/null; then
                local gcc_version=$(gcc -dumpversion | cut -d. -f1)
                if [ "$gcc_version" -ge 9 ] 2>/dev/null; then
                    print_success "System GCC $gcc_version is suitable (>= 9), using default toolchain"
                    return 0
                fi
            fi
            
            # If default GCC is too old, try to install newer version
            print_status "System GCC is too old, attempting to install GCC-11..."
            if $INSTALL_CMD gcc-11 g++-11 gfortran-11; then
                # Check if gcc-11 binaries exist before setting alternatives
                if [ -f /usr/bin/gcc-11 ] && [ -f /usr/bin/g++-11 ] && [ -f /usr/bin/gfortran-11 ]; then
                    # Set alternatives to make gcc-11 default
                    sudo update-alternatives --install /usr/bin/gcc gcc /usr/bin/gcc-11 60 \
                        --slave /usr/bin/g++ g++ /usr/bin/g++-11 \
                        --slave /usr/bin/gfortran gfortran /usr/bin/gfortran-11
                    print_success "GCC-11 installed and set as default via update-alternatives"
                else
                    print_warning "GCC-11 binaries not found in expected locations, keeping system default"
                fi
            else
                print_warning "Could not install GCC-11, keeping system default (may cause compilation issues)"
            fi
            ;;
            
        "dnf")
            # RHEL/CentOS/Fedora with dnf
            if [[ "$DISTRO" == "rhel" || "$DISTRO" == "centos" || "$DISTRO" == "rocky" || "$DISTRO" == "almalinux" ]]; then
                # First try default GCC
                $INSTALL_CMD gcc gcc-c++ gcc-gfortran
                if command -v gcc &> /dev/null; then
                    local gcc_version=$(gcc -dumpversion | cut -d. -f1)
                    if [ "$gcc_version" -ge 9 ] 2>/dev/null; then
                        print_success "System GCC $gcc_version is suitable (>= 9), using default toolchain"
                        return 0
                    fi
                fi
                
                # If default GCC is too old, try gcc-toolset
                print_status "System GCC is too old, attempting to install gcc-toolset-11..."
                if $INSTALL_CMD gcc-toolset-11-gcc gcc-toolset-11-gcc-c++ gcc-toolset-11-gcc-gfortran; then
                    # Create activation script for gcc-toolset-11
                    if [ -f /opt/rh/gcc-toolset-11/enable ]; then
                        GCC_ACTIVATION_SCRIPT="/etc/profile.d/gcc-toolset-11.sh"
                        sudo tee "$GCC_ACTIVATION_SCRIPT" > /dev/null << 'EOF'
#!/bin/bash
# Auto-activate gcc-toolset-11 for all users
if [ -f /opt/rh/gcc-toolset-11/enable ]; then
    source /opt/rh/gcc-toolset-11/enable
fi
EOF
                        sudo chmod +x "$GCC_ACTIVATION_SCRIPT"
                        
                        # Source it for current session
                        source /opt/rh/gcc-toolset-11/enable
                        print_success "GCC-toolset-11 installed and activated"
                        print_status "GCC-toolset-11 will be automatically activated for all future sessions"
                    else
                        print_warning "GCC-toolset-11 enable script not found, installation may be incomplete"
                        print_warning "Using system default GCC (may cause compilation issues if version < 9)"
                    fi
                else
                    print_warning "GCC-toolset-11 installation failed, using system default GCC"
                    print_warning "May cause compilation issues if GCC version < 9"
                fi
            else
                # Fedora - typically has recent GCC by default
                $INSTALL_CMD gcc gcc-c++ gcc-gfortran
                print_success "GCC toolchain installed (Fedora typically uses recent GCC by default)"
            fi
            ;;
            
        "yum")
            # RHEL/CentOS with yum
            # First try default GCC
            $INSTALL_CMD gcc gcc-c++ gcc-gfortran
            if command -v gcc &> /dev/null; then
                local gcc_version=$(gcc -dumpversion | cut -d. -f1)
                if [ "$gcc_version" -ge 9 ] 2>/dev/null; then
                    print_success "System GCC $gcc_version is suitable (>= 9), using default toolchain"
                    return 0
                fi
            fi
            
            # If default GCC is too old and SCL is available, try devtoolset
            if command_exists scl; then
                print_status "System GCC is too old, attempting to install devtoolset-11..."
                if $INSTALL_CMD centos-release-scl && $INSTALL_CMD devtoolset-11-gcc devtoolset-11-gcc-c++ devtoolset-11-gcc-gfortran; then
                    # Create activation script for devtoolset-11
                    if [ -f /opt/rh/devtoolset-11/enable ]; then
                        GCC_ACTIVATION_SCRIPT="/etc/profile.d/devtoolset-11.sh"
                        sudo tee "$GCC_ACTIVATION_SCRIPT" > /dev/null << 'EOF'
#!/bin/bash
# Auto-activate devtoolset-11 for all users
if [ -f /opt/rh/devtoolset-11/enable ]; then
    source /opt/rh/devtoolset-11/enable
fi
EOF
                        sudo chmod +x "$GCC_ACTIVATION_SCRIPT"
                        
                        # Source it for current session
                        source /opt/rh/devtoolset-11/enable
                        print_success "Devtoolset-11 installed and activated"
                        print_status "Devtoolset-11 will be automatically activated for all future sessions"
                    else
                        print_warning "Devtoolset-11 enable script not found, using system default GCC"
                        print_warning "May cause compilation issues if GCC version < 9"
                    fi
                else
                    print_warning "Devtoolset-11 installation failed, using system default GCC"
                    print_warning "May cause compilation issues if GCC version < 9"
                fi
            else
                print_warning "SCL not available, using system default GCC"
                print_warning "May cause compilation issues if GCC version < 9"
            fi
            ;;
            
        "zypper")
            # openSUSE
            # First try default GCC
            $INSTALL_CMD gcc gcc-c++ gcc-fortran
            if command -v gcc &> /dev/null; then
                local gcc_version=$(gcc -dumpversion | cut -d. -f1)
                if [ "$gcc_version" -ge 9 ] 2>/dev/null; then
                    print_success "System GCC $gcc_version is suitable (>= 9), using default toolchain"
                    return 0
                fi
            fi
            
            # If default GCC is too old, try to install gcc11
            print_status "System GCC is too old, attempting to install GCC-11..."
            if $INSTALL_CMD gcc11 gcc11-c++ gcc11-fortran; then
            
            # Find gcc-11 installation path
            GCC11_PATH=""
            for path in /usr/bin/gcc-11 /usr/bin/gcc11; do
                if [ -f "$path" ]; then
                    GCC11_PATH="$path"
                    break
                fi
            done
            
            if [ -n "$GCC11_PATH" ]; then
                G11_PATH="${GCC11_PATH/gcc/g++}"
                GFORTRAN11_PATH="${GCC11_PATH/gcc/gfortran}"
                
                if [ -f "$G11_PATH" ] && [ -f "$GFORTRAN11_PATH" ]; then
                    sudo update-alternatives --install /usr/bin/gcc gcc "$GCC11_PATH" 60 \
                        --slave /usr/bin/g++ g++ "$G11_PATH" \
                        --slave /usr/bin/gfortran gfortran "$GFORTRAN11_PATH"
                    print_success "GCC-11 installed and set as default"
                else
                    print_warning "GCC-11 components not fully available, using system default"
                    print_warning "May cause compilation issues if GCC version < 9"
                fi
            else
                print_warning "GCC-11 not found in expected locations, using system default"
                print_warning "May cause compilation issues if GCC version < 9"
            fi
            else
                print_warning "Could not install GCC-11, using system default"
                print_warning "May cause compilation issues if GCC version < 9"
            fi
            ;;
            
        "pacman")
            # Arch Linux
            $INSTALL_CMD gcc gcc-fortran
            print_success "GCC toolchain installed (Arch uses recent GCC by default)"
            ;;
            
        *)
            print_error "Unsupported package manager for GCC-11 installation"
            return 1
            ;;
    esac
    
    # Verify GCC installation and version
    if command_exists gcc; then
        GCC_VERSION=$(gcc -dumpversion 2>/dev/null | cut -d. -f1)
        if [ -n "$GCC_VERSION" ] && [ "$GCC_VERSION" -ge 9 ]; then
            print_success "GCC $GCC_VERSION is installed and should work with Alchemd"
        elif [ -n "$GCC_VERSION" ] && [ "$GCC_VERSION" -ge 7 ]; then
            print_warning "GCC $GCC_VERSION may work but GCC 9+ is recommended for best compatibility"
        else
            print_warning "GCC version could not be determined or may be too old (< 7)"
        fi
    else
        print_error "GCC installation verification failed"
        return 1
    fi
    
    return 0
}

# Function to install base packages
install_base_packages() {
    print_status "Installing base development packages..."
    
    case "$PKG_MANAGER" in
        "apt")
            $INSTALL_CMD tcsh make flex bison patch bc wget cmake unzip \
                xorg-dev libz-dev libbz2-dev build-essential
            ;;
        "dnf"|"yum")
            $INSTALL_CMD tcsh make flex bison patch bc wget cmake unzip \
                libX11-devel zlib-devel bzip2-devel
            ;;
        "zypper")
            $INSTALL_CMD tcsh make flex bison patch bc wget cmake unzip \
                libX11-devel zlib-devel libbz2-devel
            ;;
        "pacman")
            $INSTALL_CMD tcsh make flex bison patch bc wget cmake unzip \
                libx11 zlib bzip2
            ;;
    esac
    
    print_success "Base packages installed"
}

# Install base packages
install_base_packages

# Install GCC toolchain
if ! install_gcc_toolchain; then
    print_error "Failed to install GCC toolchain"
    exit 1
fi

# Verify installation
print_status "Verifying GCC installation..."

# Check if compilers are available
MISSING_TOOLS=""
for tool in gcc g++ gfortran make cmake; do
    if ! command_exists "$tool"; then
        if [ -z "$MISSING_TOOLS" ]; then
            MISSING_TOOLS="$tool"
        else
            MISSING_TOOLS="$MISSING_TOOLS $tool"
        fi
    fi
done

if [ -n "$MISSING_TOOLS" ]; then
    print_error "Critical tools are missing: $MISSING_TOOLS"
    print_error "Please check the installation and try again"
    exit 1
fi

# Display compiler information
print_status "Compiler information:"
echo "GCC: $(gcc --version | head -1)"
echo "G++: $(g++ --version | head -1)"
echo "Gfortran: $(gfortran --version | head -1)"
echo "Make: $(make --version | head -1)"
echo "CMake: $(cmake --version | head -1)"

# Check GCC version
GCC_VERSION=$(gcc -dumpversion | cut -d. -f1)
if [ "$GCC_VERSION" -ge 11 ]; then
    print_success "GCC version $GCC_VERSION is excellent for Alchemd"
elif [ "$GCC_VERSION" -ge 9 ]; then
    print_success "GCC version $GCC_VERSION is suitable for Alchemd"
else
    print_warning "GCC version $GCC_VERSION may be too old for optimal Alchemd compilation"
    print_warning "GCC 9 or newer is recommended, GCC 11+ is preferred"
fi

print_success "Linux environment setup completed successfully!"
print_status ""
print_status "Important notes:"
print_status "1. CUDA 12.0 support is required for Alchemd GPU acceleration"
print_status "2. Please ensure NVIDIA drivers and CUDA toolkit are installed separately"
print_status "3. GCC 9+ is required for compilation; GCC 11+ is preferred"
if [[ "$PKG_MANAGER" == "dnf" && ("$DISTRO" == "rhel" || "$DISTRO" == "centos" || "$DISTRO" == "rocky" || "$DISTRO" == "almalinux") ]]; then
    if [ -f /opt/rh/gcc-toolset-11/enable ]; then
        print_status "4. GCC-toolset-11 has been configured to activate automatically"
        print_status "   If you need to manually activate it: source /opt/rh/gcc-toolset-11/enable"
    fi
elif [[ "$PKG_MANAGER" == "yum" ]]; then
    if [ -f /opt/rh/devtoolset-11/enable ]; then
        print_status "4. Devtoolset-11 has been configured to activate automatically"
        print_status "   If you need to manually activate it: source /opt/rh/devtoolset-11/enable"
    fi
fi
print_status ""
print_status "Next steps:"
print_status "1. Install CUDA 12.0 (if not already installed)"
print_status "2. Setup Python environment: chmod +x install/setup_conda.sh && ./install/setup_conda.sh"
print_status "3. Build PMEMD: chmod +x install/build_amber.sh && ./install/build_amber.sh <pmemd_path>"