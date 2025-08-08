from rdkit import Chem
from rdkit import rdBase
import re
import os
import copy
import sys
import logging
from io import StringIO
from pdbfixer import PDBFixer
from openmm.app import PDBFile

logger = logging.getLogger(__name__)

def _deal_valence_error(_atom):
    symbol = _atom.GetSymbol()
    chg = _atom.GetFormalCharge()
    if symbol == 'C':
        expected_valence = 4
    elif symbol == 'N':
        expected_valence = 3
    elif symbol == 'O':
        expected_valence = 2
    elif symbol == 'H':
        expected_valence = 1
    else:
        return False

    neighbors_info = set()
    for b in _atom.GetBonds():
        neighbor = b.GetOtherAtom(_atom)
        neighbors_info.add((b.GetBondTypeAsDouble(), neighbor.GetSymbol()))
    sum_bond_valence = sum([i[0] for i in neighbors_info])
    new_chg = int(expected_valence - sum_bond_valence)
    logger.debug(f"{_atom.GetIdx()} {symbol} {chg} -> {new_chg}: {expected_valence} {neighbors_info}")
    if new_chg != chg:
        # logger.info(f"Deal valence error: {_atom.GetIdx()} {symbol} {chg} -> {new_chg}")
        _atom.SetFormalCharge(new_chg)
        return True

    for b in _atom.GetBonds():
        if b.GetBondTypeAsDouble() > 1:
            b.SetBondType(Chem.rdchem.BondType.SINGLE)
        neighbor = b.GetOtherAtom(_atom)
        neighbors_info.add((b.GetBondTypeAsDouble(), neighbor.GetSymbol()))
    sum_bond_valence = sum([i[0] for i in neighbors_info])
    new_chg = int(expected_valence - sum_bond_valence)
    if new_chg != chg:
        _atom.SetFormalCharge(new_chg)
        return True

    return False

def protein_pdb_reader(pdb_file, max_attempts=9999):
    """
    NOTICE!!! WARNING!!!
    this function would only care about the 3d coords of atoms.
    This function is designed for compatibility with RDKit input, so it will automatically remove some hydrogens and adjust certain bonds.
    This function should only be used for collision detection and simple docking, not for MD calculations.
    """

    ori_pdb_mol = Chem.MolFromPDBFile(pdb_file, removeHs=False, sanitize=False)

    def _raise_error(msg):
        raise RuntimeError(msg)

    attempts = 0
    while attempts < max_attempts:
        try:
            # Suppress RDKit error messages during sanitization
            rdBase.DisableLog('rdApp.error')
            Chem.SanitizeMol(ori_pdb_mol)
            rdBase.EnableLog('rdApp.error')

            return ori_pdb_mol
        except ValueError as e:
            # Re-enable RDKit error messages after exception
            rdBase.EnableLog('rdApp.error')
            error_msg = str(e)
            match = re.search(r'atom # (\d+)', error_msg, re.IGNORECASE)
            if not match:
                _raise_error(f"Unable to handle: {error_msg}")

            atom_idx = int(match.group(1))
            atom = ori_pdb_mol.GetAtomWithIdx(atom_idx)

            h_idxes = []
            for bond in atom.GetBonds():
                neighbor = bond.GetOtherAtom(atom)
                if neighbor.GetAtomicNum() == 1:
                    h_idxes.append(neighbor.GetIdx())

            if not h_idxes:
                if _deal_valence_error(atom):
                    pass
                else:
                    neighbors_info = set()
                    for bond in atom.GetBonds():
                        neighbor = bond.GetOtherAtom(atom)
                        neighbors_info.add((bond.GetBondTypeAsDouble(), neighbor.GetSymbol()))
                    _raise_error(f"Unable to handle: {atom_idx} {atom.GetSymbol()} {atom.GetFormalCharge()} {neighbors_info}")

            h_idxes.sort(reverse=True)
            emol = Chem.EditableMol(ori_pdb_mol)
            for h_idx in h_idxes:
                emol.RemoveAtom(h_idx)
            ori_pdb_mol = emol.GetMol()
            attempts += 1

    _raise_error("Failed to read pdb file.")

def standardize_pdb(input_pdb):
    # Generate output file name
    base_name = os.path.splitext(input_pdb)[0]
    output_pdb = f"{base_name}_std.pdb"
    
    # Fix PDB
    fixer = PDBFixer(filename=input_pdb)
    
    # Standardize atom names
    fixer.findNonstandardResidues()
    fixer.replaceNonstandardResidues()
    
    # List of atom types to standardize
    standard_atoms = ['C', 'N', 'O', 'H', 'S', 'P']
    
    # Standardize atom names
    for residue in fixer.topology.residues():
        for atom in residue.atoms():
            # Get current atom name
            current_name = atom.name
            # Check if standardization is needed
            for std_atom in standard_atoms:
                if current_name.startswith(std_atom) and len(current_name) > 1:
                    if current_name[1].isupper():
                        atom.name = std_atom
                        break
    
    # Save fixed file
    with open(output_pdb, 'w') as f:
        PDBFile.writeFile(fixer.topology, fixer.positions, f)
    
    print(f"Stanardized: {output_pdb}")


def pdb_info_clean(pdb_file, res_name=None):
    tmp_content = ''
    with open(pdb_file, 'r') as f:
        for line in f:
            if not (line.startswith('HETATM') or line.startswith('ATOM')):
                continue
            if res_name is not None:
                line = line.replace('UNL', res_name)
                line = f"HETATM{line[6:]}"
            tmp_content += line

    with open(pdb_file, 'w') as f:
        f.write(tmp_content)


def pdb_formatter(pdb_file, atom_name_formatter: str = 'title') -> str:
    new_pdbx_content = ""

    def _get_char(_str):
        return ''.join(c for c in _str if c.isalpha())

    def _get_num(_str):
        return ''.join(c for c in _str if c.isnumeric())

    def _format_char(_str: str, _format):
        if _format.lower() == 'title':
            return _str.title()
        elif _format.lower() == 'upper':
            return _str.upper()

    with open(pdb_file, 'r') as f:
        for line in f.readlines():
            if line.startswith('HETATM') or line.startswith('ATOM'):
                # tmp_line = line
                # while tmp_line.find("  ") != -1:
                #     tmp_line = tmp_line.replace("  ", ' ')
                atom_text = line[12:17].strip()
                atom_name = _get_char(atom_text)
                atom_id = _get_num(atom_text)
                formatted_atom_name = _format_char(atom_name, atom_name_formatter)
                boo_line = (f'{line[:12]}{formatted_atom_name:>2}{atom_id:<2}{" "}{line[17:54]}'
                        f'{1.00:>6}{0.00:>6}{" "*6}{" "*4}{formatted_atom_name:>2}')
                if len(formatted_atom_name) == 1:
                    boo_line = f'{boo_line}{line[-3:]}'
                else:
                    boo_line = f'{boo_line}{line[-2:]}'
                # TODO: here we discard the occupancy\tempFactor\segID info in pdb.
                new_pdbx_content += boo_line
            else:
                new_pdbx_content += line

    with open(pdb_file, 'w') as f:
        f.write(new_pdbx_content)

    return os.path.abspath(pdb_file)


def calculate_formal_charge(_mol: Chem.Mol):
    # TODO: Add charge information to pdb file.
    _full_H_mol = Chem.AddHs(copy.deepcopy(_mol))
    _common_valence_dict = {'C': 4, 'N': 3, 'O': 2}
    formal_charge_list = []
    given_formal_charge_list = []
    for atom in _mol.GetAtoms():
        given_formal_charge_list.append(atom.GetFormalCharge())
    total_given_formal_charge = sum(given_formal_charge_list)
    return int(total_given_formal_charge)


def pdb_reindex(pdb_file):
    with open(pdb_file, 'r') as f:
        lines = f.readlines()
    new_content = ""

    chain_res_num_mapping = {}
    chain_first_res_num = {}

    for line in lines:
        if line.startswith('ATOM') or line.startswith('HETATM'):
            try:
                chain_id = line[21:22].strip()
                res_num = int(line[22:26].strip())

                if chain_id not in chain_res_num_mapping:
                    chain_res_num_mapping[chain_id] = {}
                    chain_first_res_num[chain_id] = res_num

                if res_num not in chain_res_num_mapping[chain_id]:
                    chain_res_num_mapping[chain_id][res_num] = len(chain_res_num_mapping[chain_id]) + chain_first_res_num[chain_id]
            except ValueError:
                pass

    for line in lines:
        if line.startswith('ATOM') or line.startswith('HETATM'):
            try:

                chain_id = line[21:22].strip()
                res_num = int(line[22:26].strip())

                if chain_id in chain_res_num_mapping and res_num in chain_res_num_mapping[chain_id]:
                    new_line = line[:22] + f"{chain_res_num_mapping[chain_id][res_num]:4d}" + line[26:]
                    new_content += new_line
                else:
                    new_content += line
            except ValueError:
                new_content += line
        else:
            new_content += line
    
    with open(pdb_file, 'w') as f:
        f.write(new_content)
    
    return os.path.abspath(pdb_file)


if __name__ == "__main__":
    pass