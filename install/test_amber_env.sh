#!/bin/bash

# Test script for Amber environment auto-setup
# This script helps verify that the automatic Amber environment integration is working

set -e

# Color codes for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

print_status() {
    echo -e "${BLUE}[TEST]${NC} $1"
}

print_success() {
    echo -e "${GREEN}[PASS]${NC} $1"
}

print_warning() {
    echo -e "${YELLOW}[WARN]${NC} $1"
}

print_error() {
    echo -e "${RED}[FAIL]${NC} $1"
}

print_header() {
    echo ""
    echo "=================================="
    echo "$1"
    echo "=================================="
}

# Test environment detection
test_environment_detection() {
    print_header "Testing Environment Detection"
    
    # Check for conda/mamba
    if command -v mamba &> /dev/null; then
        print_success "mamba found: $(which mamba)"
    elif command -v conda &> /dev/null; then
        print_success "conda found: $(which conda)"
    else
        print_error "Neither conda nor mamba found"
        return 1
    fi
    
    # Check current environment
    if [[ -n "$CONDA_PREFIX" ]]; then
        print_success "Current environment: ${CONDA_DEFAULT_ENV:-$(basename "$CONDA_PREFIX")}"
        print_success "Environment path: $CONDA_PREFIX"
    else
        print_error "No conda/mamba environment currently activated"
        return 1
    fi
    
    return 0
}

# Test activation scripts
test_activation_scripts() {
    print_header "Testing Activation Scripts"
    
    if [[ -z "$CONDA_PREFIX" ]]; then
        print_error "No environment activated"
        return 1
    fi
    
    local activate_script="$CONDA_PREFIX/etc/conda/activate.d/amber_init.sh"
    local deactivate_script="$CONDA_PREFIX/etc/conda/deactivate.d/amber_init.sh"
    
    # Check activation script
    if [[ -f "$activate_script" ]]; then
        print_success "Activation script exists: $activate_script"
        if [[ -x "$activate_script" ]]; then
            print_success "Activation script is executable"
        else
            print_warning "Activation script is not executable"
        fi
    else
        print_error "Activation script not found: $activate_script"
        return 1
    fi
    
    # Check deactivation script
    if [[ -f "$deactivate_script" ]]; then
        print_success "Deactivation script exists: $deactivate_script"
        if [[ -x "$deactivate_script" ]]; then
            print_success "Deactivation script is executable"
        else
            print_warning "Deactivation script is not executable"
        fi
    else
        print_error "Deactivation script not found: $deactivate_script"
        return 1
    fi
    
    return 0
}

# Test AMBERHOME setting
test_amberhome() {
    print_header "Testing AMBERHOME Environment Variable"
    
    if [[ -n "$AMBERHOME" ]]; then
        print_success "AMBERHOME is set: $AMBERHOME"
        
        if [[ -d "$AMBERHOME" ]]; then
            print_success "AMBERHOME directory exists (primary AmberTools location)"
        else
            print_error "AMBERHOME directory does not exist: $AMBERHOME"
            return 1
        fi
        
        # Check for amber.sh in primary location
        if [[ -f "$AMBERHOME/amber.sh" ]]; then
            print_success "AmberTools environment script found: $AMBERHOME/amber.sh"
        else
            print_warning "AmberTools environment script not found: $AMBERHOME/amber.sh"
        fi
        
    else
        print_error "AMBERHOME is not set"
        return 1
    fi
    
    return 0
}

# Test Amber25 dual-component setup
test_amber25_components() {
    print_header "Testing Amber25 Dual-Component Setup"
    
    local components_found=0
    
    # Try to detect AmberTools and PMEMD installations from activation script
    local activate_script=""
    if [[ -n "$CONDA_PREFIX" ]]; then
        activate_script="$CONDA_PREFIX/etc/conda/activate.d/amber_init.sh"
    fi
    
    if [[ -f "$activate_script" ]]; then
        print_success "Amber25 activation script found: $activate_script"
        
        # Extract paths from activation script
        local ambertools_path=$(grep "export AMBERHOME=" "$activate_script" | cut -d'"' -f2)
        local pmemd_paths=$(grep -o '"/[^"]*pmemd[^"]*"' "$activate_script" | tr -d '"')
        
        if [[ -n "$ambertools_path" && -d "$ambertools_path" ]]; then
            print_success "AmberTools component detected: $ambertools_path"
            ((components_found++))
            
            if [[ -f "$ambertools_path/amber.sh" ]]; then
                print_success "AmberTools environment script exists"
            else
                print_warning "AmberTools environment script missing"
            fi
        fi
        
        for pmemd_path in $pmemd_paths; do
            if [[ -d "$pmemd_path" ]]; then
                print_success "PMEMD component detected: $pmemd_path"
                ((components_found++))
                
                if [[ -f "$pmemd_path/amber.sh" ]]; then
                    print_success "PMEMD environment script exists"
                else
                    print_warning "PMEMD environment script missing"
                fi
                break
            fi
        done
    else
        print_warning "Amber25 activation script not found: $activate_script"
    fi
    
    if [[ $components_found -ge 2 ]]; then
        print_success "Both Amber25 components detected correctly"
        return 0
    else
        print_warning "Incomplete Amber25 dual-component setup (found $components_found/2 components)"
        return 1
    fi
}

# Test Amber tools
test_amber_tools() {
    print_header "Testing Amber Tools Availability"
    
    # AmberTools components
    local ambertools_bins=("sander" "antechamber" "tleap")
    local ambertools_found=0
    
    print_status "Checking AmberTools binaries..."
    for tool in "${ambertools_bins[@]}"; do
        if command -v "$tool" &> /dev/null; then
            print_success "$tool found: $(which "$tool")"
            ((ambertools_found++))
        else
            print_warning "$tool not found in PATH"
        fi
    done
    
    # PMEMD components
    local pmemd_bins=("pmemd" "pmemd.cuda")
    local pmemd_found=0
    
    print_status "Checking PMEMD binaries..."
    for tool in "${pmemd_bins[@]}"; do
        if command -v "$tool" &> /dev/null; then
            print_success "$tool found: $(which "$tool")"
            ((pmemd_found++))
        else
            if [[ "$tool" == "pmemd.cuda" ]]; then
                print_warning "$tool not found (GPU support may be unavailable)"
            else
                print_warning "$tool not found in PATH"
            fi
        fi
    done
    
    local total_found=$((ambertools_found + pmemd_found))
    local total_expected=$((${#ambertools_bins[@]} + ${#pmemd_bins[@]}))
    
    print_status "Found $total_found out of $total_expected Amber tools"
    print_status "  AmberTools: $ambertools_found/${#ambertools_bins[@]}"
    print_status "  PMEMD: $pmemd_found/${#pmemd_bins[@]}"
    
    # Consider test passed if we have most tools (at least 3 out of 5)
    if [[ $total_found -ge 3 ]]; then
        return 0
    else
        return 1
    fi
}

# Test Alchemd integration
test_alchemd_integration() {
    print_header "Testing Alchemd Integration"
    
    # Check if we're in the Alchemd directory
    if [[ -f "core/run_md.py" ]]; then
        print_success "Alchemd core directory found"
        
        # Test if run_md.py can show help (basic functionality test)
        if python core/run_md.py --help &> /dev/null; then
            print_success "run_md.py is functional"
        else
            print_warning "run_md.py may have issues (this might be normal if dependencies are missing)"
        fi
    else
        print_warning "Not in Alchemd root directory or core/run_md.py not found"
    fi
    
    return 0
}

# Main test function
main() {
    print_header "Alchemd PMEMD Environment Integration Test"
    print_status "This script tests if the automatic PMEMD environment setup is working correctly"
    
    local tests_passed=0
    local total_tests=6
    
    # Run tests
    if test_environment_detection; then
        ((tests_passed++))
    fi
    
    if test_activation_scripts; then
        ((tests_passed++))
    fi
    
    if test_amberhome; then
        ((tests_passed++))
    fi
    
    if test_amber25_components; then
        ((tests_passed++))
    fi
    
    if test_amber_tools; then
        ((tests_passed++))
    fi
    
    if test_alchemd_integration; then
        ((tests_passed++))
    fi
    
    # Final results
    print_header "Test Results"
    print_status "Passed: $tests_passed / $total_tests tests"
    
    if [[ $tests_passed -eq $total_tests ]]; then
        print_success "All tests passed! Amber25 dual-component environment integration is working correctly."
        echo ""
        print_status "You can now use Alchemd directly without manual Amber setup:"
        print_status "  python core/run_md.py --help"
        print_status "  python core/prepare_file.py --help"
        print_status "  python core/analyze_result.py --help"
        return 0
    elif [[ $tests_passed -ge 4 ]]; then
        print_warning "Most tests passed. There might be minor issues but basic functionality should work."
        return 0
    else
        print_error "Multiple tests failed. Please check your Amber25 installation and environment setup."
        echo ""
        print_status "Troubleshooting suggestions:"
        print_status "1. Make sure you're in a conda/mamba environment"
        print_status "2. Re-run the setup: bash install/setup_amber_env.sh /path/to/ambertools25 /path/to/pmemd24"
        print_status "3. Check if both Amber components were compiled successfully"
        print_status "4. Deactivate and reactivate your environment: conda deactivate && conda activate your_env"
        print_status "5. Verify both AmberTools and PMEMD directories exist and contain amber.sh files"
        return 1
    fi
}

# Run main function
main "$@"