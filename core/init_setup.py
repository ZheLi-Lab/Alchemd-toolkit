#!/usr/bin/env python3
"""
Alchemd Configuration Initialization Tool

This script provides a command-line interface to initialize Alchemd configuration
without requiring the GUI. It automatically detects project paths, validates 
dependencies, and generates the configs.toml file.

Based on the GUI SettingPage logic from gui/gui.py.

Author: Alchemd Development Team
"""

import os
import sys
import argparse
import toml
from pathlib import Path
from typing import Dict, List, Optional, Tuple
import subprocess


class Colors:
    """ANSI color codes for terminal output."""
    RED = '\033[0;31m'
    GREEN = '\033[0;32m'
    YELLOW = '\033[1;33m'
    BLUE = '\033[0;34m'
    CYAN = '\033[0;36m'
    NC = '\033[0m'  # No Color


def print_status(message: str) -> None:
    """Print a status message in blue."""
    print(f"{Colors.BLUE}[INFO]{Colors.NC} {message}")


def print_success(message: str) -> None:
    """Print a success message in green."""
    print(f"{Colors.GREEN}[SUCCESS]{Colors.NC} {message}")


def print_warning(message: str) -> None:
    """Print a warning message in yellow."""
    print(f"{Colors.YELLOW}[WARNING]{Colors.NC} {message}")


def print_error(message: str) -> None:
    """Print an error message in red."""
    print(f"{Colors.RED}[ERROR]{Colors.NC} {message}")


def print_header(title: str) -> None:
    """Print a section header."""
    print(f"\n{Colors.CYAN}{'='*60}{Colors.NC}")
    print(f"{Colors.CYAN}{title.center(60)}{Colors.NC}")
    print(f"{Colors.CYAN}{'='*60}{Colors.NC}\n")


class ConfigInitializer:
    """
    Configuration initializer for Alchemd.
    
    Based on the GUI SettingPage class logic, this provides command-line
    configuration setup functionality.
    """
    
    def __init__(self, project_root: Optional[str] = None):
        """
        Initialize the configuration setup.
        
        Args:
            project_root: Optional project root path. If None, auto-detected.
        """
        self.project_root = Path(project_root) if project_root else self._detect_project_root()
        self.config_file = self.project_root / "core" / "configs.toml"
        self.dependencies_dir = self.project_root / "dependencies"
        
        # Configuration dictionary - matches GUI SettingPage structure
        self.config_dict: Dict[str, str] = {}
        self.toml_dict: Dict[str, Dict[str, str]] = {}
        
    def _detect_project_root(self) -> Path:
        """
        Auto-detect the Alchemd project root directory.
        
        Returns:
            Path to the project root directory.
            
        Raises:
            RuntimeError: If project root cannot be detected.
        """
        # Start from the script's directory and work upward
        current_dir = Path(__file__).parent.parent.absolute()
        
        # Look for characteristic files/directories
        markers = ["core", "gui", "dependencies", "install"]
        
        for potential_root in [current_dir] + list(current_dir.parents):
            if all((potential_root / marker).exists() for marker in markers):
                return potential_root
                
        # Fallback: assume script is in core/ directory
        return Path(__file__).parent.parent.absolute()
    
    def setup_default_paths(self) -> None:
        """
        Set up default paths based on project structure.
        Mirrors the setup_paths() method from GUI SettingPage.
        """
        print_status("Setting up default paths...")
        
        # Base paths
        toolkit_path = str(self.project_root)
        deps_path = str(self.dependencies_dir)
        
        # Default configuration - matches GUI SettingPage defaults
        self.config_dict = {
            'AlchemdToolkit_Path': toolkit_path,
            'CAR_Path': os.path.join(deps_path, 'CAR-FEP'),
            'GenambRBFE_Path': os.path.join(deps_path, 'genambRBFE'),
            'AlchemdConvTools_Path': os.path.join(deps_path, 'AlchemConvTools'),
            'Wcc_Path': os.path.join(deps_path, 'Weighted_cc-main'),
            'Watvina_Path': os.path.join(deps_path, 'watvina', 'watvina'),
            'Alchemd_Path': os.path.join(deps_path, 'Alchemd'),
        }
        
        print_success("Default paths configured")
        
    def interactive_setup(self) -> None:
        """
        Interactive configuration setup for advanced users.
        """
        print_header("Interactive Configuration Setup")
        print("You can customize the paths or press Enter to use defaults.\n")
        
        # Display and optionally modify each path
        for key, default_value in self.config_dict.items():
            display_name = key.replace('_', ' ').replace('Path', '').strip()
            
            print(f"{Colors.CYAN}{display_name}:{Colors.NC}")
            print(f"  Default: {default_value}")
            
            user_input = input(f"  Custom path (or Enter for default): ").strip()
            
            if user_input:
                self.config_dict[key] = os.path.abspath(user_input)
                print_success(f"Updated {display_name}")
            else:
                print_status(f"Using default for {display_name}")
            print()
    
    def validate_dependencies(self) -> List[str]:
        """
        Validate that dependencies are available and functional.
        Based on the test_dependencies() method from GUI SettingPage.
        
        Returns:
            List of failed dependencies.
        """
        print_status("Validating dependencies...")
        failed_dependencies = []
        
        def test_program_executable(path: str, program: str, is_binary: bool = False) -> bool:
            """Test if a program is executable."""
            if is_binary:
                exec_path = os.path.join(path, program)
                if not os.path.isfile(exec_path):
                    return False
                if not os.access(exec_path, os.X_OK):
                    return False
                # Try to run with --help or -h
                try:
                    result = subprocess.run([exec_path, '--help'], 
                                          capture_output=True, timeout=10)
                    return result.returncode == 0
                except (subprocess.TimeoutExpired, subprocess.SubprocessError, FileNotFoundError):
                    try:
                        result = subprocess.run([exec_path, '-h'], 
                                              capture_output=True, timeout=10)
                        return result.returncode == 0
                    except (subprocess.TimeoutExpired, subprocess.SubprocessError, FileNotFoundError):
                        return False
            else:
                # Python script test
                script_path = os.path.join(path, program)
                if not os.path.isfile(script_path):
                    return False
                try:
                    result = subprocess.run([sys.executable, script_path, '-h'], 
                                          capture_output=True, timeout=10)
                    return result.returncode == 0
                except (subprocess.TimeoutExpired, subprocess.SubprocessError):
                    return False
        
        # Test GenambRBFE
        genambrbfe_path = self.config_dict['GenambRBFE_Path']
        if not os.path.exists(os.path.join(genambrbfe_path, 'genambRBFE')):
            failed_dependencies.append('GenambRBFE')
            print_warning("GenambRBFE not found")
        else:
            print_success("GenambRBFE found")
        
        # Test WatVina
        watvina_path = os.path.dirname(self.config_dict['Watvina_Path'])
        watvina_exec = os.path.basename(self.config_dict['Watvina_Path'])
        if not test_program_executable(watvina_path, watvina_exec, is_binary=True):
            failed_dependencies.append('WatVina')
            print_warning("WatVina not executable or not found")
        else:
            print_success("WatVina is functional")
        
        # Test Weighted_cc
        wcc_path = self.config_dict['Wcc_Path']
        if not test_program_executable(wcc_path, 'wcc_main.py'):
            failed_dependencies.append('Weighted_cc')
            print_warning("Weighted_cc not functional")
        else:
            print_success("Weighted_cc is functional")
        
        # Test CAR
        car_path = self.config_dict['CAR_Path']
        if not test_program_executable(car_path, 'segmented_converge_control.py'):
            failed_dependencies.append('CAR')
            print_warning("CAR not functional")
        else:
            print_success("CAR is functional")
        
        # Test AlchemConvTools
        alchemconv_path = self.config_dict['AlchemdConvTools_Path']
        if not test_program_executable(alchemconv_path, 'one_end_fe_aly.py'):
            failed_dependencies.append('AlchemConvTools')
            print_warning("AlchemConvTools not functional")
        else:
            print_success("AlchemConvTools is functional")
        
        return failed_dependencies
    
    def check_amber_environment(self) -> bool:
        """
        Check if AMBER environment is properly configured.
        
        Returns:
            True if AMBER is properly configured, False otherwise.
        """
        print_status("Checking AMBER environment...")
        
        amberhome = os.environ.get('AMBERHOME')
        if not amberhome:
            print_warning("AMBERHOME environment variable not set")
            print_warning("Make sure to run: source /path/to/pmemd/amber.sh")
            return False
        
        print_success(f"AMBERHOME found: {amberhome}")
        
        # Check for critical AMBER binaries
        amber_binaries = ['pmemd.cuda', 'pmemd.cuda_SPFP']
        missing_binaries = []
        
        for binary in amber_binaries:
            binary_path = os.path.join(amberhome, 'bin', binary)
            if os.path.isfile(binary_path) and os.access(binary_path, os.X_OK):
                print_success(f"{binary} found and executable")
            else:
                missing_binaries.append(binary)
                print_warning(f"{binary} not found or not executable")
        
        if missing_binaries:
            print_warning("Some critical AMBER binaries are missing")
            print_warning("Please ensure PMEMD compilation completed successfully")
            return False
        
        return True
    
    def generate_config(self) -> None:
        """
        Generate the TOML configuration structure.
        """
        print_status("Generating configuration structure...")
        
        # Add AMBERHOME if available
        amberhome = os.environ.get('AMBERHOME')
        if amberhome:
            self.config_dict['AMBERHOME'] = amberhome
        
        # Create TOML structure - matches GUI structure
        self.toml_dict = {
            'paths': self.config_dict
        }
        
        print_success("Configuration structure generated")
    
    def save_config(self) -> bool:
        """
        Save configuration to configs.toml file.
        
        Returns:
            True if saved successfully, False otherwise.
        """
        print_status("Saving configuration to configs.toml...")
        
        try:
            # Ensure core directory exists
            self.config_file.parent.mkdir(parents=True, exist_ok=True)
            
            # Save configuration
            with open(self.config_file, 'w') as f:
                toml.dump(self.toml_dict, f)
            
            print_success(f"Configuration saved to: {self.config_file}")
            return True
            
        except Exception as e:
            print_error(f"Failed to save configuration: {e}")
            return False
    
    def display_summary(self, failed_deps: List[str]) -> None:
        """
        Display configuration summary.
        
        Args:
            failed_deps: List of failed dependencies.
        """
        print_header("Configuration Summary")
        
        print(f"{Colors.CYAN}Project Root:{Colors.NC} {self.project_root}")
        print(f"{Colors.CYAN}Config File:{Colors.NC} {self.config_file}")
        print()
        
        print(f"{Colors.CYAN}Configured Paths:{Colors.NC}")
        for key, value in self.config_dict.items():
            status = "✓" if os.path.exists(value) or key == 'AMBERHOME' else "✗"
            print(f"  {status} {key}: {value}")
        print()
        
        if failed_deps:
            print(f"{Colors.YELLOW}Failed Dependencies:{Colors.NC}")
            for dep in failed_deps:
                print(f"  ✗ {dep}")
            print()
            print_warning("Some dependencies failed validation")
            print_warning("These tools may not work correctly until properly installed")
        else:
            print_success("All dependencies validated successfully!")
        
        print()
        print(f"{Colors.CYAN}Next Steps:{Colors.NC}")
        print("1. Activate conda environment: conda activate alchemd")
        print("2. Test core modules:")
        print("   python core/prepare_file.py --help")
        print("   python core/run_md.py --help")
        print("   python core/analyze_result.py --help")
        print("3. Test with example data:")
        print("   python core/prepare_file.py -p example/protein.pdb -l example/ligands.sdf -o test_output")
    
    def run(self, auto_mode: bool = False, interactive: bool = False) -> bool:
        """
        Run the configuration initialization process.
        
        Args:
            auto_mode: If True, run in automatic mode without user interaction.
            interactive: If True, run interactive setup for path customization.
            
        Returns:
            True if configuration was successful, False otherwise.
        """
        print_header("Alchemd Configuration Initialization")
        
        print_status(f"Project root detected: {self.project_root}")
        print_status(f"Configuration file: {self.config_file}")
        
        # Check if config already exists
        if self.config_file.exists() and not auto_mode:
            overwrite = input(f"\nConfiguration file already exists. Overwrite? (y/N): ").strip().lower()
            if overwrite != 'y':
                print_status("Configuration cancelled by user")
                return False
        
        try:
            # Setup default paths
            self.setup_default_paths()
            
            # Interactive setup if requested
            if interactive and not auto_mode:
                self.interactive_setup()
            
            # Validate dependencies
            failed_deps = self.validate_dependencies()
            
            # Check AMBER environment
            amber_ok = self.check_amber_environment()
            if not amber_ok and not auto_mode:
                continue_anyway = input("\nAMBER environment issues detected. Continue anyway? (y/N): ").strip().lower()
                if continue_anyway != 'y':
                    print_error("Configuration cancelled due to AMBER environment issues")
                    return False
            
            # Generate and save configuration
            self.generate_config()
            if not self.save_config():
                return False
            
            # Display summary
            if not auto_mode:
                self.display_summary(failed_deps)
            
            return True
            
        except KeyboardInterrupt:
            print_error("\nConfiguration cancelled by user")
            return False
        except Exception as e:
            print_error(f"Configuration failed: {e}")
            return False


def main():
    """Main entry point for the configuration script."""
    parser = argparse.ArgumentParser(
        description="Initialize Alchemd configuration",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  %(prog)s                    # Interactive setup
  %(prog)s --auto            # Automatic setup with defaults
  %(prog)s --interactive     # Interactive path customization
  %(prog)s --project-root /path/to/alchemd  # Custom project root
        """
    )
    
    parser.add_argument(
        '--auto', 
        action='store_true',
        help='Run in automatic mode without user interaction'
    )
    
    parser.add_argument(
        '--interactive', 
        action='store_true',
        help='Run interactive setup for path customization'
    )
    
    parser.add_argument(
        '--project-root',
        type=str,
        help='Specify custom project root directory'
    )
    
    parser.add_argument(
        '--version',
        action='version',
        version='Alchemd Configuration Tool v1.0'
    )
    
    args = parser.parse_args()
    
    try:
        # Initialize configuration
        initializer = ConfigInitializer(project_root=args.project_root)
        
        # Run configuration process
        success = initializer.run(
            auto_mode=args.auto,
            interactive=args.interactive
        )
        
        if success:
            print_success("Configuration initialization completed successfully!")
            sys.exit(0)
        else:
            print_error("Configuration initialization failed")
            sys.exit(1)
            
    except Exception as e:
        print_error(f"Unexpected error: {e}")
        sys.exit(1)


if __name__ == "__main__":
    main()