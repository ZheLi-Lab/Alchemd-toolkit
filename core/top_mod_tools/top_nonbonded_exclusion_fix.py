import parmed
import argparse
import os

def read_atom_list(filename):
    with open(filename, 'r') as f:
        lines = f.readlines()
    
    line0 = lines[0].strip().replace(":LGA@", "").replace(":LGB@", "")
    line1 = lines[1].strip().replace(":LGA@", "").replace(":LGB@", "")

    lam_atoms = [i for i in line0.split(",") if i != ""]
    lbm_atoms = [i for i in line1.split(",") if i != ""]

    return lam_atoms, lbm_atoms

def find_atoms_in_structure(structure, res_name, atom_names):
    found_atoms = []
    for atom in structure.atoms:
        if atom.residue.name == res_name and atom.name in atom_names:
            found_atoms.append(atom)
    
    return found_atoms

def check_and_fix_exclusions(structure, lam_atoms_list, lbm_atoms_list):

    lam_atoms = find_atoms_in_structure(structure, 'LAM', lam_atoms_list)
    lbm_atoms = find_atoms_in_structure(structure, 'LBM', lbm_atoms_list)

    # print(f"Found {len(lam_atoms)} LAM atoms and {len(lbm_atoms)} LBM atoms")
    # print(f"LAM atoms: {[i.idx for i in lam_atoms]}")
    # print(f"LBM atoms: {[i.idx for i in lbm_atoms]}")

    fixed_count = 0
    for lam_atom in lam_atoms:
        for lbm_atom in lbm_atoms:
            lam_exlcusions = [i.idx for i in lam_atom.exclusion_partners]
            lbm_exclusions = [i.idx for i in lbm_atom.exclusion_partners]
            if lbm_atom.idx in lam_exlcusions or lam_atom.idx in lbm_exclusions:
                continue
            else:
                lam_atom.exclude(lbm_atom)
                fixed_count += 1
    return fixed_count


def check_and_fix_exclusions_from_file(top_file, mutant_list_file, output_file=None,
                                       verbose=False, if_save=True):
    if output_file is None:
        output_file = top_file.replace(".prmtop", "_fixed.prmtop")

    structure = parmed.amber.LoadParm(top_file)

    lam_atoms_list, lbm_atoms_list = read_atom_list(mutant_list_file)

    fixed_count = check_and_fix_exclusions(structure, lam_atoms_list, lbm_atoms_list)

    if fixed_count > 0:
        if verbose:
            print(f"Moded {fixed_count} LAM-LBM exclusions")
        if if_save:
            structure.save(output_file, overwrite=True)
        return True

    if verbose:
        print("All exclusions set properly")
    return False


def main():
    parser = argparse.ArgumentParser(description='Fix nonbonded_exclusions in Amber topology file')
    parser.add_argument('prmtop', help='Input Amber topology file (.prmtop)')
    parser.add_argument('atom_list', help='Text file containing LAM and LBM residue atom names')
    parser.add_argument('--output', '-o', help='Output fixed topology file name')

    args = parser.parse_args()
    args.prmtop = os.path.abspath(args.prmtop)
    args.atom_list = os.path.abspath(args.atom_list)
    output_file = args.output if args.output else f"{os.path.splitext(args.prmtop)[0]}_fixed.prmtop"

    check_and_fix_exclusions_from_file(args.prmtop, args.atom_list, output_file,
                                       verbose=True, if_save=True)


if __name__ == "__main__":
    main()