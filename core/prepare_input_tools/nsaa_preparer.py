import os
import subprocess
import argparse
import shutil
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem
import parmed as pmd
import copy

# --- Configuration ---
STANDARD_AMINO_ACIDS = [
    "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "HIE",
    "HID", "HIP", "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR",
    "TRP", "TYR", "VAL",
]
EXCLUDED_RESIDUES = ["ACE", "NME", "HOH", "NMA", "CYX", "HID", "HIE", "HIP",
                     "ASH", "GLH", "LYN", ]
# AmberTools executable paths (modify if not in system PATH)
ANTECHAMBER_EXEC = "antechamber"
PARMCHK2_EXEC = "parmchk2"
TLEAP_EXEC = "tleap"
# GAFF version: "gaff" or "gaff2"
GAFF_VERSION = "gaff2"

# --- Helper functions ---
def run_command(command_list, working_dir="."):
    """Execute external command and return its output"""
    print(f"Executing: {' '.join(command_list)} in {working_dir}")
    try:
        process = subprocess.Popen(command_list, stdout=subprocess.PIPE, stderr=subprocess.PIPE, cwd=working_dir, text=True)
        stdout, stderr = process.communicate(timeout=300) # 5 minute timeout
        if process.returncode != 0:
            print(f"Error executing command: {' '.join(command_list)}")
            print(f"STDOUT:\n{stdout}")
            print(f"STDERR:\n{stderr}")
            return False, stdout, stderr
        print(f"STDOUT:\n{stdout}")
        if stderr: # Some tools output normal information to stderr
             print(f"STDERR:\n{stderr}")
        return True, stdout, stderr
    except subprocess.TimeoutExpired:
        print(f"Command timed out: {' '.join(command_list)}")
        process.kill()
        stdout, stderr = process.communicate()
        return False, stdout, stderr
    except Exception as e:
        print(f"Failed to run command {' '.join(command_list)}: {e}")
        return False, "", str(e)

def calculate_residue_charge(pdb_path, residue_number, residue_name):
    """Use RDKit to calculate FormalCharge for specified residue number
    
    Args:
        pdb_path: PDB file path
        residue_number: Residue number (starting from 1)
        residue_name: Residue name
    """
    try:
        # Read PDB file without sanitization
        mol = Chem.MolFromPDBFile(pdb_path, sanitize=False)
        if mol is None:
            print(f"Warning: Could not read PDB file: {pdb_path}")
            return 0
            
        # Get formal charges of all atoms
        total_charge = 0
        for atom in mol.GetAtoms():
            # Get residue information for the atom
            res_info = atom.GetPDBResidueInfo()
            if res_info is not None and res_info.GetResidueNumber() == residue_number and res_info.GetResidueName() == residue_name:
                total_charge += atom.GetFormalCharge()
            
        return total_charge
    except Exception as e:
        print(f"Warning: Could not calculate charge for residue number {residue_number}: {e}")
        return 0

def validate_antechamber_params(pdb_path, charge, gaff_version):
    """Validate antechamber input parameters"""
    if not os.path.exists(pdb_path):
        return False, f"PDB file not found: {pdb_path}"
    
    if not isinstance(charge, (int, float)):
        return False, f"Invalid charge value: {charge}"
    
    if gaff_version not in ["gaff", "gaff2"]:
        return False, f"Invalid GAFF version: {gaff_version}"
    
    return True, ""

def write_pdb_with_ordered_atoms(input_pdb, output_pdb):
    """Reorder atoms in PDB file to ensure hydrogen atoms are continuous with their connected residues
    
    Args:
        input_pdb: Input PDB file path
        output_pdb: Output PDB file path
    """
    try:
        # Read all ATOM lines
        atom_lines = []
        with open(input_pdb, 'r') as f:
            for line in f:
                if line.startswith('ATOM') or line.startswith('HETATM'):
                    atom_lines.append(line)
        
        # Group by residue
        residue_atoms = {}
        for line in atom_lines:
            # Parse residue information
            res_name = line[17:20].strip()
            res_num = int(line[22:26])
            chain_id = line[21]
            ins_code = line[26]
            atom_name = line[12:16].strip()
            is_hydrogen = atom_name.startswith('H')
            
            # Use residue information as key
            res_key = (res_name, res_num, chain_id, ins_code)
            if res_key not in residue_atoms:
                residue_atoms[res_key] = {'heavy': [], 'hydrogen': []}
            
            # Add atom line to corresponding group
            if is_hydrogen:
                residue_atoms[res_key]['hydrogen'].append(line)
            else:
                residue_atoms[res_key]['heavy'].append(line)
        
        # Write atoms in residue order
        with open(output_pdb, 'w') as f:
            atom_idx = 1
            for res_key in sorted(residue_atoms.keys()):
                # Write heavy atoms first
                for line in residue_atoms[res_key]['heavy']:
                    # Update atom number
                    new_line = line[:6] + f"{atom_idx:5d}" + line[11:]
                    f.write(new_line)
                    atom_idx += 1
                # Then write hydrogen atoms
                for line in residue_atoms[res_key]['hydrogen']:
                    # Update atom number
                    new_line = line[:6] + f"{atom_idx:5d}" + line[11:]
                    f.write(new_line)
                    atom_idx += 1
            # Write END line
            f.write("END\n")
            
    except Exception as e:
        print(f"Warning: Error reordering atoms in PDB file: {e}")
        # If error occurs, directly copy original file
        shutil.copy2(input_pdb, output_pdb)

def extract_residue_with_neighbors(struct, residue_idx):
    """Extract residue and its neighboring residues
    
    Args:
        struct: parmed.Structure object
        residue_idx: Target residue index
    
    Returns:
        tuple: (parmed.Structure, dict) New structure containing target residue and its neighboring residues, and atom number mapping
    """
    try:
        # Create a new Structure object
        new_struct = pmd.Structure()
        # Used to record mapping from original atom numbers to new atom numbers
        atom_number_map = {}
        
        # Get target residue
        target_residue = struct.residues[residue_idx]
        
        # Look forward for one neighboring residue
        start_idx = residue_idx
        if residue_idx > 0:
            prev_res = struct.residues[residue_idx - 1]
            # If not TER or HOH, and residue numbers are consecutive, include previous residue
            if prev_res.name != "HOH" and not (hasattr(prev_res, 'ter') and prev_res.ter):
                if (prev_res.number == target_residue.number - 1):
                    start_idx = residue_idx - 1
                elif (prev_res.number == target_residue.number):
                    start_idx = residue_idx - 1
        
        # Look backward for one neighboring residue
        end_idx = residue_idx
        if residue_idx < len(struct.residues) - 1:
            next_res = struct.residues[residue_idx + 1]
            # If not TER or HOH, and residue numbers are consecutive, include next residue
            if next_res.name != "HOH" and not (hasattr(next_res, 'ter') and next_res.ter) and \
               next_res.number == target_residue.number + 1:
                end_idx = residue_idx + 1
        
        # Add residue and its neighboring residues
        new_atom_number = 1
        for i in range(start_idx, end_idx + 1):
            residue = struct.residues[i]
            for atom in residue.atoms:
                # Record atom number mapping
                atom_number_map[new_atom_number] = atom.number
                # Directly copy original atom
                new_atom = copy.copy(atom)
                # Only modify atom number
                new_atom.number = new_atom_number
                # Add to new structure
                new_struct.add_atom(new_atom, residue.name, residue.number)
                new_atom_number += 1
        
        return new_struct, atom_number_map
    except Exception as e:
        print(f"Warning: Error extracting residue with neighbors: {e}")
        return None, None

def read_pdb_charges(pdb_file, target_residues):
    """Read atomic charge information for specified residues from PDB file
    
    Args:
        pdb_file: PDB file path
        target_residues: List of residues that need charge recording, each element is a (resname, resnum) tuple
    
    Returns:
        dict: Mapping from atom number to charge
    """
    charges = {}
    try:
        with open(pdb_file, 'r') as f:
            for line in f:
                if line.startswith('ATOM') or line.startswith('HETATM'):
                    # In PDB format, atom number is in columns 7-11
                    atom_num = int(line[6:11].strip())
                    # Residue name is in columns 18-20
                    resname = line[17:20].strip()
                    # Residue number is in columns 23-26
                    resnum = int(line[22:26].strip())
                    
                    # Only process atoms of target residues
                    if (resname, resnum) in target_residues:
                        charge_str = line[-3:].strip()
                        if len(charge_str) < 2:
                            continue
                        try:
                            positive = charge_str.endswith('+')
                            if positive:
                                charge = int(float(charge_str[:-1]))
                            else:
                                charge = -int(float(charge_str[:-1]))
                            charges[atom_num] = charge
                            print(f"Atom {atom_num} in {resname}{resnum} has {charge_str} charge {charge}")
                        except ValueError:
                            charges[atom_num] = 0
    except Exception as e:
        print(f"Warning: Error reading charges from PDB file: {e}")
    return charges

def fix_pdb_charges(pdb_file, charges, atom_number_map):
    """Fix charge information in PDB file
    
    Args:
        pdb_file: PDB file path that needs fixing
        charges: Mapping from original atom number to charge
        atom_number_map: Mapping from new atom number to original atom number
    """
    try:
        with open(pdb_file, 'r') as f:
            lines = f.readlines()
        
        with open(pdb_file, 'w') as f:
            for line in lines:
                if line.startswith('ATOM') or line.startswith('HETATM'):
                    new_atom_num = int(line[6:11].strip())
                    if new_atom_num in atom_number_map:
                        original_atom_num = atom_number_map[new_atom_num]
                        if original_atom_num in charges:
                            # Keep most of original line content, only modify last 3 columns
                            new_line = line[:-3]
                            charge = charges[original_atom_num]
                            if charge > 0:
                                charge_str = f"{charge}+"
                            else:
                                charge_str = f"{abs(charge)}-"
                            new_line += charge_str.ljust(3)
                            f.write(new_line + '\n')
                        else:
                            f.write(line)
                    else:
                        f.write(line)
                else:
                    f.write(line)
    except Exception as e:
        print(f"Warning: Error fixing charges in PDB file: {e}")

def find_residue_number_in_pdb(pdb_file, target_resname):
    """Find residue number for specified Residue name in PDB file
    
    Args:
        pdb_file: PDB file path
        target_resname: Target Residue name
    
    Returns:
        int: Found residue number, returns None if not found
    """
    try:
        with open(pdb_file, 'r') as f:
            for line in f:
                if line.startswith('ATOM') or line.startswith('HETATM'):
                    resname = line[17:20].strip()
                    if resname == target_resname:
                        return int(line[22:26].strip())
        return None
    except Exception as e:
        print(f"Warning: Error finding residue number in PDB file: {e}")
        return None

# --- Main logic ---
def process_pdb_for_nsaa(pdb_file, output_dir, keep_temp=False):
    """
    Process PDB file, identify non-standard amino acids, and generate force field parameters for them.
    
    Args:
        pdb_file: Input PDB file path
        output_dir: Output directory
        keep_temp: Whether to keep temporary files
    """
    if not os.path.exists(pdb_file):
        print(f"Error: PDB file '{pdb_file}' not found.")
        return

    if os.path.exists(output_dir):
        print(f"Warning: Output directory '{output_dir}' already exists. Files might be overwritten.")
    else:
        os.makedirs(output_dir, exist_ok=True)
        print(f"Created output directory: {output_dir}")

    # 1. Use parmed to load PDB and identify NSAA
    try:
        # Load structure
        struct = pmd.load_file(pdb_file)
        
    except Exception as e:
        print(f"Error loading PDB file '{pdb_file}' with parmed: {e}")
        return

    unique_nsaa = {} # Store {resname: (residx, resnum)}
    print("\nScanning for Non-Standard Amino Acids (NSAAs)...")
    for i, res in enumerate(struct.residues):
        if res.name == "HOH":
            continue
        if res.name not in STANDARD_AMINO_ACIDS and res.name not in EXCLUDED_RESIDUES:
            if res.name not in unique_nsaa:
                # Store index and residue number of first encountered NSAA
                unique_nsaa[res.name] = (i, res.number)
                print(f"  Found NSAA: {res.name} (Index: {i}, Residue number: {res.number})")

    if not unique_nsaa:
        print("No non-standard amino acids found in the PDB file.")
        return

    print(f"\nFound {len(unique_nsaa)} unique NSAA(s): {', '.join(unique_nsaa.keys())}")
    print("--- IMPORTANT NOTES ---")
    print("1. This script extracts NSAAs AS IS. For peptide-linked NSAAs, manual capping is CRUCIAL before parameterization.")
    print("2. Ensure the NSAA structure (hydrogens, geometry) is reasonable before proceeding.")
    print("3. AM1-BCC charges will be used by default. For higher accuracy, consider RESP charges (requires QM).")
    print("4. Generated parameters should be carefully validated.")
    print("-----------------------\n")

    # 2. Generate parameters for each unique NSAA
    for nsaa_name, (nsaa_res_idx, nsaa_res_num) in unique_nsaa.items():
        print(f"\nProcessing NSAA: {nsaa_name}...")
        
        # Create a temporary working directory for current NSAA to avoid filename conflicts
        nsaa_work_dir = os.path.join(output_dir, f"temp_{nsaa_name}")
        os.makedirs(nsaa_work_dir, exist_ok=True)

        # a. Extract NSAA to separate PDB file
        try:
            # First extract residue and its neighboring residues
            print(f"  Extracting {nsaa_name} with neighboring residues...")
            nsaa_struct, atom_number_map = extract_residue_with_neighbors(struct, nsaa_res_idx)
            if nsaa_struct is None:
                print(f"Error: Could not extract residue {nsaa_name} with neighbors. Skipping.")
                if not keep_temp:
                    shutil.rmtree(nsaa_work_dir)
                continue

            # Get list of residues that need charge recording
            target_residues = []
            for residue in nsaa_struct.residues:
                target_residues.append((residue.name, residue.number))
            
            # Read charge information for these residues
            atom_charges = read_pdb_charges(pdb_file, target_residues)
            
            # Save extracted structure
            extracted_pdb_path = os.path.join(nsaa_work_dir, f"{nsaa_name}_extracted.pdb")
            nsaa_struct.save(extracted_pdb_path, overwrite=True)
            
            # Fix charge information in saved PDB file
            fix_pdb_charges(extracted_pdb_path, atom_charges, atom_number_map)
            
            # Find renumbered residue number
            new_res_num = find_residue_number_in_pdb(extracted_pdb_path, nsaa_name)
            if new_res_num is None:
                print(f"  Error: Could not find residue {nsaa_name} in extracted PDB file. Skipping.")
                if not keep_temp:
                    shutil.rmtree(nsaa_work_dir)
                continue
            
            print(f"  Extracted {nsaa_name} with neighbors to {extracted_pdb_path}")
            print(f"  New residue number for {nsaa_name} in extracted PDB: {new_res_num}")
        except Exception as e:
            print(f"  Error extracting {nsaa_name} with parmed: {e}. Skipping.")
            if not keep_temp:
                shutil.rmtree(nsaa_work_dir)
            continue

        # b. Calculate net charge
        charge = calculate_residue_charge(extracted_pdb_path, new_res_num, nsaa_name)
        print(f"  Calculated charge for {nsaa_name}: {charge}")

        # c. Copy and modify template file
        template_path = os.path.join(os.path.dirname(os.path.dirname(__file__)), "templates", "nsaa_prep.txt")
        if not os.path.exists(template_path):
            print(f"Error: Template file not found at {template_path}")
            if not keep_temp:
                shutil.rmtree(nsaa_work_dir)
            continue

        # Copy template file to working directory
        prep_file = os.path.join(nsaa_work_dir, f"{nsaa_name}_prep.txt")
        shutil.copy2(template_path, prep_file)
        if 'GAUSS_SCRDIR' not in os.environ:
            raise EnvironmentError("Environment variable GAUSS_SCRDIR is not set. Please set the GAUSS_SCRDIR environment variable to the directory of Gaussian temporary files.")
        temp_g03_dir = os.environ['GAUSS_SCRDIR']
        os.makedirs(temp_g03_dir, exist_ok=True)

        # Read and modify template file content
        with open(prep_file, 'r') as f:
            content = f.read()

        # Replace placeholders
        content = content.replace("__NSAA_PDB__", os.path.basename(extracted_pdb_path))
        content = content.replace("__NSAA_RESNAME__", nsaa_name)
        content = content.replace("__NSAA_RESIDX__", str(new_res_num))
        content = content.replace("__NSAA_CHG__", str(charge))
        content = content.replace("__TEMP_DIR__", temp_g03_dir)

        # Write back to file
        with open(prep_file, 'w') as f:
            f.write(content)

        print(f"  Created and modified preparation file: {prep_file}")

        # d. Execute pyautomd command
        try:
            # Save current working directory
            current_dir = os.getcwd()
            # Switch to working directory
            os.chdir(nsaa_work_dir)
            
            # Execute pyautomd command
            success, stdout, stderr = run_command(["pyautomd", "-i", f"{nsaa_name}_prep.txt"])
            if not success:
                print(f"  Error executing pyautomd for {nsaa_name}:")
                print(f"  STDOUT: {stdout}")
                print(f"  STDERR: {stderr}")
                os.chdir(current_dir)
                if not keep_temp:
                    shutil.rmtree(nsaa_work_dir)
                continue
            
            # Copy generated files to output_dir
            frcmod_file = os.path.join(nsaa_work_dir, f"{nsaa_name}.frcmod")
            prepi_file = os.path.join(nsaa_work_dir, f"{nsaa_name}.prepi")
            nsaa_out_dir = os.path.join(output_dir, nsaa_name)
            os.makedirs(nsaa_out_dir, exist_ok=True)
            
            if os.path.exists(frcmod_file):
                shutil.copy2(frcmod_file, os.path.join(nsaa_out_dir, f"{nsaa_name}.frcmod"))
                print(f"  Copied {nsaa_name}.frcmod to output directory")
            else:
                print(f"  Warning: {nsaa_name}.frcmod not found")
                
            if os.path.exists(prepi_file):
                shutil.copy2(prepi_file, os.path.join(nsaa_out_dir, f"{nsaa_name}.prepi"))
                print(f"  Copied {nsaa_name}.prepi to output directory")
            else:
                print(f"  Warning: {nsaa_name}.prepi not found")
            
            # Switch back to original directory
            os.chdir(current_dir)
            
        except Exception as e:
            print(f"  Error during pyautomd execution for {nsaa_name}: {e}")
            if not keep_temp:
                shutil.rmtree(nsaa_work_dir)
            continue

        # f. Clean up temporary working directory
        if not keep_temp:
            print(f"  Cleaning up temporary files for {nsaa_name}...")
            shutil.rmtree(nsaa_work_dir)
        else:
            print(f"  Keeping temporary files in {nsaa_work_dir}")
        print(f"  Finished processing {nsaa_name}.")

    print("\nAll processing finished.")
    print(f"Generated files are in: {os.path.abspath(output_dir)}")
    print("Please carefully review all generated files and validate parameters before use in simulations.")
    print("\nTo use the generated parameters in tleap, add the following lines to your tleap script:")
    print("source leaprc.gaff2  # or leaprc.gaff for GAFF")
    for nsaa_name in unique_nsaa.keys():
        print(f"loadamberprep {nsaa_name.upper()}.prep")
        print(f"loadamberparams {nsaa_name.upper()}.frcmod")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Process a PDB file to identify Non-Standard Amino Acids (NSAA) "
                                                 "and generate Amber force field parameters (.frcmod, .lib) for them.")
    parser.add_argument("pdb_file", help="Path to the input PDB file.")
    parser.add_argument("-o", "--output_dir", default="nsaa_parameters",
                        help="Directory to save the generated parameter files (default: nsaa_parameters).")
    parser.add_argument("-nd", "--no-delete", action="store_true",
                        help="Do not delete temporary files after processing.")
    
    args = parser.parse_args()

    # Process output_dir
    if args.output_dir == "nsaa_parameters":
        output_dir = os.path.abspath(os.path.join(os.getcwd(), "nsaa_parameters"))
    else:
        output_dir = os.path.abspath(args.output_dir)

    # Ensure AmberTools executables exist
    missing_tools = []
    for tool in [ANTECHAMBER_EXEC, PARMCHK2_EXEC, TLEAP_EXEC]:
        if shutil.which(tool) is None:
            missing_tools.append(tool)
    if missing_tools:
        print(f"Error: The following AmberTools executables were not found in your PATH: {', '.join(missing_tools)}")
        print("Please ensure AmberTools is correctly installed and its bin directory is in your PATH.")
        exit(1)
        
    process_pdb_for_nsaa(args.pdb_file, output_dir, args.no_delete)