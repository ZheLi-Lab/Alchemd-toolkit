import numpy as np
import warnings
import logging
import argparse
from scipy.spatial import distance
from rdkit import Chem
from typing import List, Union, Tuple
import os, subprocess, shutil

try:
    from pdb_func import pdb_reindex, protein_pdb_reader
except:
    from .pdb_func import pdb_reindex, protein_pdb_reader

logger = logging.getLogger(__name__)

class CollisionChecker:
    def __init__(self, protein_mol: Chem.Mol) -> None:
        self.protein_mol = protein_mol
        self.protein_coords = []
        _conf = self.protein_mol.GetConformer()
        for idx in range(_conf.GetNumAtoms()):
            p3d = _conf.GetAtomPosition(idx)
            self.protein_coords.append([p3d.x, p3d.y, p3d.z])
        self.protein_coords = np.array(self.protein_coords)

    def get_collision_idxes(self, ligand_mol: Chem.Mol, threshold: float = 0.7) -> List[Tuple[int, int, int]]:
        _conf = ligand_mol.GetConformer()
        ligand_coords = []
        for idx in range(_conf.GetNumAtoms()):
            p3d = _conf.GetAtomPosition(idx)
            ligand_coords.append([p3d.x, p3d.y, p3d.z])
        ligand_coords = np.array(ligand_coords)

        dist_matrix = distance.cdist(self.protein_coords, ligand_coords)
        collision_mask = dist_matrix < threshold
        collision_count = np.sum(collision_mask)

        if collision_count == 0:
            return []

        coll_nonlig_idxes, coll_lig_idxes = np.where(collision_mask)
        flat_idxes = np.argsort(dist_matrix.flatten())[:len(coll_lig_idxes)]
        lig_idxes, nonlig_idxes = np.unravel_index(flat_idxes, dist_matrix.shape)
        return_list = []
        for i in range(len(lig_idxes)):
            l_idx = coll_lig_idxes[i]
            p_idx = coll_nonlig_idxes[i]
            dist = dist_matrix[p_idx, l_idx]
            # if (self.protein_mol.GetAtomWithIdx(int(p_idx)).GetSymbol() == "H" and
            #     ligand_mol.GetAtomWithIdx(int(l_idx)).GetSymbol() == "H"):
            #     if dist > 0.5:
            #         continue
            return_list.append((p_idx, l_idx, dist))

        return return_list
    
    @classmethod
    def from_pdb_file(cls, pdb_file: str) -> "CollisionChecker":
        protein_mol = protein_pdb_reader(pdb_file)
        return cls(protein_mol)
    
    @staticmethod
    def format_output(coll_idxes: List[Tuple[int, int]]) -> str:
        total_coll_count = len(coll_idxes)
        if total_coll_count == 0:
            return "No collision detected."
        return f"Found {total_coll_count}. " + \
               f"Closest: protein atom {coll_idxes[0][0]} -> ligand atom {coll_idxes[0][1]} ({coll_idxes[0][2]:.2f} A)"
    
    def get_moles_collision_idxes(self, mol_list: List[Chem.Mol]) -> List[List[Tuple[int, int, int]]]:
        return_list = []
        for mol in mol_list:
            return_list.append(self.get_collision_idxes(mol))
        return return_list
  
    @staticmethod  
    def get_water_remove_pdb_from_file(pdb_file: str, moles: List[Chem.Mol], output_pdb_file: str) -> bool:
        tmp_collision_checker = CollisionChecker.from_pdb_file(pdb_file)
        return tmp_collision_checker.get_water_remove_pdb(moles, output_pdb_file)

    def get_water_remove_pdb(self, moles: List[Chem.Mol], output_pdb_file: str,
                             if_reprepare: bool = False) -> bool:
        moles_collision_idxes = self.get_moles_collision_idxes(moles)
        if_collision = False
        for collision_idxes in moles_collision_idxes:
            if len(collision_idxes) > 0:
                if_collision = True
                break
        if not if_collision:
            return False

        #remove collide water
        idxes_to_remove = set()
        for mol_collision_idxes in moles_collision_idxes:
            if len(mol_collision_idxes) == 0:
                continue
            for collision_idx in mol_collision_idxes:
                protein_collide_idx = int(collision_idx[0])
                atom = self.protein_mol.GetAtomWithIdx(protein_collide_idx)
                residue_info = atom.GetPDBResidueInfo()
                if residue_info:
                    residue_name = residue_info.GetResidueName()
                    if residue_name == "HOH" or residue_name == "WAT":
                        idxes_to_remove.add(protein_collide_idx)
                        if atom.GetSymbol() == "H":
                            idxes_to_remove.add(atom.GetNeighbors()[0].GetIdx())
                            idxes_to_remove.update([i.GetIdx() for i in atom.GetNeighbors()[0].GetNeighbors()])
                        else:
                            # atom is O
                            idxes_to_remove.update([i.GetIdx() for i in atom.GetNeighbors()])

        if len(idxes_to_remove) == 0:
            # collision but not WAT
            warnings.warn("Ligand might collide with protein, but not a water molecule.")
            return False

        logger.info(f'remove {int(len(idxes_to_remove)/3)} water')
        rwmol = Chem.RWMol(self.protein_mol)
        rwmol.BeginBatchEdit()
        for idx in idxes_to_remove:
            rwmol.RemoveAtom(idx)
        rwmol.CommitBatchEdit()
        Chem.MolToPDBFile(rwmol, output_pdb_file, flavor=0b001010)

        pdb_reindex(output_pdb_file)
        self._add_ter(output_pdb_file)

        if if_reprepare:
            logger.debug("Reprepare.")
            re_prepared_pdb = output_pdb_file.replace(".pdb, _amb.pdb")
            cmd = f'pdb2pqr {output_pdb_file} out.pdb\
            --titration-state-method=propka\
            --with-ph=7.45 --ffout=AMBER\
            --pdb-output={re_prepared_pdb}\
            --display-coupled-residues'
            result = subprocess.run(cmd, shell=True, 
                                    stdout=subprocess.PIPE, 
                                    stderr=subprocess.PIPE,
                                    universal_newlines=True,
                                    )
            shutil.copy2(re_prepared_pdb, output_pdb_file)

        return True
        
    @staticmethod
    def _add_ter(pdb_file: str) -> None:
        with open(pdb_file, 'r') as f:
            lines = f.readlines()

        temp_file = pdb_file + '.temp'
        with open(temp_file, 'w') as f:
            in_hoh_block = False
            prev_residue = None
            last_line_was_ter = False

            for i, line in enumerate(lines):
                if line.startswith(('ATOM', 'HETATM')):
                    residue_name = line[17:20].strip()

                    if residue_name == 'HOH':
                        if not in_hoh_block and prev_residue != 'HOH' and not last_line_was_ter:
                            f.write('TER\n')
                            last_line_was_ter = True
                            in_hoh_block = True
                        else:
                            last_line_was_ter = False
                        f.write(line)
                    else:
                        if in_hoh_block and not last_line_was_ter:
                            f.write('TER\n')
                            last_line_was_ter = True
                            in_hoh_block = False
                        else:
                            last_line_was_ter = False
                        f.write(line)
                    
                    prev_residue = residue_name
                else:
                    if line.strip() == 'TER':
                        last_line_was_ter = True
                    else:
                        last_line_was_ter = False
                    f.write(line)

        os.replace(temp_file, pdb_file)
        logger.info(f"Added TER markers around HOH residues in {pdb_file}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='detect collision between protein and ligand')
    parser.add_argument('pdb_file', type=str, help='protein pdb file')
    parser.add_argument('ligand_file', type=str, help='ligand file, (pdb, sdf, mol2)')
    parser.add_argument('-t', '--threshold', type=float, default=0.7, help='collision threshold (default: 0.7)')
    parser.add_argument('-o', '--output_pdb_file', type=str, default=None, help='output pdb file')
    args = parser.parse_args()
    
    tp = CollisionChecker.from_pdb_file(args.pdb_file)

    if args.ligand_file.endswith('.sdf'):
        mol_list = Chem.SDMolSupplier(args.ligand_file, sanitize=False, removeHs=False)
    elif args.ligand_file.endswith('.pdb'):
        mol_list = [Chem.MolFromPDBFile(args.ligand_file, removeHs=False, sanitize=False)]
    elif args.ligand_file.endswith('.mol2'):
        mol_list = [Chem.MolFromMol2File(args.ligand_file, removeHs=False, sanitize=False)]
    else:
        raise ValueError(f"Unsupported ligand file format: {args.ligand_file}")
    check_dict = {}
    for mol_count, mol in enumerate(mol_list):
        collision_idxes = tp.get_collision_idxes(mol, threshold=args.threshold)
        try:
            mol_name = mol.GetProp("_Name")
        except:
            mol_name = f"ligand_{mol_count}"
        check_dict[mol_name] = collision_idxes
    counter = 0
    for key, value in check_dict.items():
        if len(value) > 0:
            counter += 1
            print(f"Collision detected for {key}: {tp.format_output(value)}")
    print(f"Total collided molecules: {counter}")

    if args.output_pdb_file is not None:
        tp.get_water_remove_pdb(mol_list, args.output_pdb_file)
