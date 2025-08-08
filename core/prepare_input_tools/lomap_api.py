import copy
import os
from contextlib import contextmanager
from rdkit import Chem
from rdkit.Chem import rdForceFieldHelpers
from rdkit.Chem import rdDistGeom
from rdkit.Chem import rdPartialCharges
from typing import List
import lomap

try:
    from pair_func import Pair
except:
    from .pair_func import Pair

@contextmanager
def working_dir(mkdir_,):
    pwd = os.getcwd()
    if os.path.exists(mkdir_) and os.path.isdir(mkdir_):
        pass
    else:
        os.makedirs(mkdir_)
    os.chdir(mkdir_)
    yield
    os.chdir(pwd)


class Mol_2_mol2:
    def __init__(self, filename):
        self._filename = filename
        self._3d_mol = self.mol_preparation_stream(self._filename)

    @staticmethod
    def mol_preparation_stream(filename):
        """
        Prepare molecule from file, add 3D coordinates and optimize geometry
        
        Parameters:
        filename: Input file path (.smi format)
        
        Returns:
        RDKit molecule object with 3D coordinates
        """
        if os.path.exists(filename) and os.path.isfile(filename):
            if filename.endswith('.smi'):
                with open(filename, 'r') as f:
                    smiles = f.readline().strip()
                    try:
                        _mol = Chem.MolFromSmiles(smiles)
                        if _mol is None:
                            raise ValueError('Invalid SMILES string: {}'.format(smiles))
                    except Exception:
                        raise ValueError('Invalid SMILES string: {}'.format(smiles))
        _update_mol = copy.deepcopy(_mol)
        _update_mol = Chem.RemoveAllHs(_update_mol)
        _update_mol = Chem.AddHs(_update_mol)

        # Generate 3D conformation
        rdDistGeom.EmbedMolecule(_update_mol)
        
        # Optimize geometry using UFF
        rdForceFieldHelpers.UFFOptimizeMolecule(_update_mol)
        
        # Further optimize using MMFF94 if parameters are available
        if rdForceFieldHelpers.MMFFHasAllMoleculeParams(_update_mol):
            rdForceFieldHelpers.MMFFOptimizeMolecule(_update_mol)
        return _update_mol
    
    @staticmethod
    def determine_tripos_atom_type(atom):
        """
        Determine Tripos atom type based on RDKit atom information
        
        Parameters:
        atom: RDKit atom object
        
        Returns:
        str: Tripos atom type
        """
        atomic_num = atom.GetAtomicNum()
        atom_idx = atom.GetIdx()
        
        # Hydrogen
        if atomic_num == 1:
            return "H"
        
        # Halogens
        if atomic_num == 9:
            return "F"
        if atomic_num == 17:
            return "Cl"
        if atomic_num == 35:
            return "Br"
        if atomic_num == 53:
            return "I"
        
        # Carbon
        if atomic_num == 6:
            if atom.GetIsAromatic():
                return "C.ar"
            
            # Check hybridization state
            hyb = atom.GetHybridization()
            if hyb == Chem.HybridizationType.SP3:
                return "C.3"
            elif hyb == Chem.HybridizationType.SP2:
                # Check if it's a carbonyl carbon
                for bond in atom.GetBonds():
                    if (bond.GetBondType() == Chem.BondType.DOUBLE and
                            bond.GetOtherAtom(atom).GetAtomicNum() == 8):
                        return "C.2"
                return "C.2"
            elif hyb == Chem.HybridizationType.SP:
                return "C.1"
        
        # Nitrogen
        if atomic_num == 7:
            if atom.GetIsAromatic():
                return "N.ar"
            
            # Check if it's an amide nitrogen
            is_amide = False
            for bond in atom.GetBonds():
                neighbor = bond.GetOtherAtom(atom)
                if neighbor.GetAtomicNum() == 6:  # Carbon atom
                    for neighbor_bond in neighbor.GetBonds():
                        if (neighbor_bond.GetBondType() == Chem.BondType.DOUBLE and
                                neighbor_bond.GetOtherAtom(neighbor).GetAtomicNum() == 8):
                            is_amide = True
                            break
            if is_amide:
                return "N.am"
            
            hyb = atom.GetHybridization()
            if hyb == Chem.HybridizationType.SP3:
                return "N.3"
            elif hyb == Chem.HybridizationType.SP2:
                return "N.2"
            elif hyb == Chem.HybridizationType.SP:
                return "N.1"
        
        # Oxygen
        if atomic_num == 8:
            if atom.GetTotalDegree() == 1:  # One connection
                for bond in atom.GetBonds():
                    if bond.GetBondType() == Chem.BondType.DOUBLE:
                        return "O.2"
                        
            # Check if it's a carboxyl oxygen
            is_carboxyl = False
            for bond in atom.GetBonds():
                neighbor = bond.GetOtherAtom(atom)
                if neighbor.GetAtomicNum() == 6:  # Carbon atom
                    carboxyl_count = 0
                    for neighbor_bond in neighbor.GetBonds():
                        other_atom = neighbor_bond.GetOtherAtom(neighbor)
                        if other_atom.GetAtomicNum() == 8:
                            carboxyl_count += 1
                    if carboxyl_count == 2:
                        is_carboxyl = True
                        break
            if is_carboxyl:
                return "O.co2"
                
            return "O.3"
        
        # Sulfur
        if atomic_num == 16:
            hyb = atom.GetHybridization()
            if hyb == Chem.HybridizationType.SP3:
                return "S.3"
            elif hyb == Chem.HybridizationType.SP2:
                return "S.2"
        
        # Phosphorus
        if atomic_num == 15:
            hyb = atom.GetHybridization()
            if hyb == Chem.HybridizationType.SP3:
                return "P.3"
        
        # Default: return atom symbol
        return atom.GetSymbol()
    
    @staticmethod
    def assign_tripos_types(mol):
        """
        Assign Tripos types for all atoms in the molecule
        
        Parameters:
        mol: RDKit molecule object
        
        Returns:
        list: List of Tripos atom types
        """
        tripos_types = []
        for atom in mol.GetAtoms():
            tripos_type = Mol_2_mol2.determine_tripos_atom_type(atom)
            tripos_types.append(tripos_type)
            # Optional: Store type as atom property
            atom.SetProp('TriposType', tripos_type)

    @staticmethod
    def write_mol2_file(mol: Chem.Mol, output_file):
        """
        Write molecule to MOL2 file format with Tripos atom types and charges
        
        Parameters:
        mol: RDKit molecule object
        output_file: Output file path
        """
        Chem.SanitizeMol(mol)
        mol.UpdatePropertyCache()
        # Assign Tripos atom types
        Mol_2_mol2.assign_tripos_types(mol)
        
        # compute Gasteiger charges
        rdPartialCharges.ComputeGasteigerCharges(mol)
        
        # Get conformer for 3D coordinates
        conf = mol.GetConformer()
        
        with open(output_file, 'w') as f:
            # Write header
            f.write("@<TRIPOS>MOLECULE\n")
            f.write(f"{mol.GetProp('_Name') if mol.HasProp('_Name') else 'Molecule'}\n")
            f.write(f"{mol.GetNumAtoms()} {mol.GetNumBonds()} 0 0 0\n")
            f.write("SMALL\nGASTEIGER\n\n")
            
            # Write atom block
            f.write("@<TRIPOS>ATOM\n")
            for i, atom in enumerate(mol.GetAtoms()):
                pos = conf.GetAtomPosition(i)
                atom_line = "{:>7d} {:<10s}{:>8.4f}  {:>8.4f}  {:>8.4f} {:<6s}{:>3d}  MOL        {:>7.4f}\n".format(
                    i+1,
                    atom.GetSymbol(),
                    pos.x,
                    pos.y,
                    pos.z,
                    atom.GetProp('TriposType'),
                    1,
                    float(atom.GetProp('_GasteigerCharge'))
                )
                f.write(atom_line)
            
            # Write bond block
            f.write("\n@<TRIPOS>BOND\n")
            for i, bond in enumerate(mol.GetBonds()):
                bond_type = {
                    Chem.BondType.SINGLE: '1',
                    Chem.BondType.DOUBLE: '2',
                    Chem.BondType.TRIPLE: '3',
                    Chem.BondType.AROMATIC: 'ar'
                }.get(bond.GetBondType(), '1')
                
                f.write("{:>6d}{:>6d}{:>6d}{:>5s}\n".format(
                    i+1,
                    bond.GetBeginAtomIdx()+1,
                    bond.GetEndAtomIdx()+1,
                    bond_type
                ))


class Lomap_api():
    def __init__(self, lomap_path,
                 lomap_out_path='lomap_out',
                 parallel=4,
                 verbose='info',
                 time=20,
                 ecrscore=0.0, # For the alchemical transformation are usually performed between molecules with the same charges, ecrscore=0.0 
                 threed=False,
                 max3d=1000,
                 element_change=True, # Always True
                 output=True, # Always True
                 name='out', # Always 'out'
                 output_no_images=False, # Always False
                 output_no_graph=False, # Always False
                 display=False, # Always False
                 allow_tree=False,
                 _max=6,
                 max_dist_from_actives=2,
                 cutoff=0.4,
                 radial=False,
                 hub=None, # For non-radial graph, no need to set
                 fast=False, # For non-radial graph, no need to set
                 links_file='',
                 known_actives_file='',
                 # chunk_mode=False, # For the fastlomap, may not need to set this, keep the default value
                 # chunk_scale = 10, # For the fastlomap, may not need to set this
                 # chunk_terminate_factor = 2.0, # For the fastlomap, may not need to set this
                 # node_mode=False, # For the fastlomap, may not need to set this
                 ):
        """
        Initialize Lomap API

        Parameters
        ----------
        lomap_path : str
           the mol2/sdf directory file name
        parallel : int
           the number of cores used to generate the similarity score matrices
        verbose : str
           verbose mode, one of 'off'/'info'/'pedantic'
        time : int
           the maximum time in seconds used to perform the MCS search
        ecrscore: float
           the electrostatic score to be used (if != 0) if two molecule have
           different charges
        threed: bool
           If true, symmetry-equivalent MCSes are filtered to prefer the one
           with the best real-space alignment
        max3d: float
           The MCS is filtered to remove atoms which are further apart than
           this threshold. The default of 1000 is effectively "no filter"
        element_change: bool
           Whether to allow changes in elements between two mappings.
           Defaults to True
        output : bool
           a flag used to generate or not the output files
        name : str
           the file name prefix used to produce the output files
        output_no_images : bool
           a flag used to disable the generation of the output image files
        output_no_graph : bool
           a flag used to disable the generation of the output graph (.dot)
           file
        display : bool
           a flag used to display or not a network made by using matplotlib
        allow_tree: bool
           if set, then the final graph does not need a cycle covering and will
           be a tree
        _max : int
           the maximum diameter of the resulting graph
        max_dist_from_actives: int
           The maximum distance of any molecule from an active (requires known_actives_file)
        cutoff : float
           the Minimum Similarity Score (MSS) used to build the graph
        radial: bool
            If using the radial option to build the graph
        hub: str or None
            Using a radial graph approach with a manually specified hub compound
        links_file : str
           the name of a file containing links to seed the graph with; 'Specify a filename listing the pairs of molecule files that should be initialised as linked.'
                              'Each line can be "mol1 mol2", which indicates that Lomap should compute the score and mapping '
                              'but that this link must be used in the final graph, or "mol1 mol2 score", which indicates that '
                              'Lomap should use the provided score, or "mol1 mol2 score force", which indicates that Lomap '
                              'should use the provided score and force this link to be used in the final graph.'
        known_actives_file : str
           the name of a file containing mols whose activity is known
        """
        self._lomap_path = os.path.abspath(lomap_path)
        # print(f"Lomap input files path: {self._lomap_path}")
        self._lomap_out_path = lomap_out_path
        with working_dir(self._lomap_out_path):
            self._db_mol = lomap.DBMolecules(self._lomap_path, parallel=parallel, verbose=False, time=time,
                                             ecrscore=ecrscore, threed=threed, max3d=max3d,
                                             element_change=element_change, output=output, name=name,
                                             output_no_images=output_no_images, output_no_graph=output_no_graph,
                                             display=display, allow_tree=allow_tree, max=_max,
                                             max_dist_from_actives=max_dist_from_actives, cutoff=cutoff, radial=radial,
                                             hub=hub, fast=fast, links_file=links_file,
                                             known_actives_file=known_actives_file,)
                                             # chunk_mode=chunk_mode, chunk_scale=chunk_scale, chunk_terminate_factor=chunk_terminate_factor, node_mode=node_mode)
            # Generate the strict and loose symmetric similarity 
            # score matrices
            strict, loose = self._db_mol.build_matrices()
            # Generate the NetworkX graph and output the results
            nx_graph = self._db_mol.build_graph()

    def output_pairs_lst(self, output_pairs_lst_file=None):
        """
        Generate pairs list for Lomap result
        
        Parameters:
        output_pairs_lst_file: Output pairs list file path
        """
        input_file = os.path.join(self._lomap_out_path, 'out_score_with_connection.txt')
        with open(input_file, 'r') as f:
            lines = f.readlines()
        
        # Skip header line
        pairs = []
        for line in lines[1:]:
            parts = line.strip().split(',')
            if len(parts) < 8:
                continue
                
            connect = parts[7].strip()
            if connect == 'Yes':
                file1 = parts[2].strip().split('.')[0]  # Remove .mol2
                file2 = parts[3].strip().split('.')[0]  # Remove .mol2
                pairs.append(Pair(file1, file2))

        if output_pairs_lst_file is None:
            output_pairs_lst_file = os.path.join(self._lomap_out_path, 'pairs.lst')
        else:
            output_pairs_lst_file = os.path.abspath(output_pairs_lst_file)
        # Write pairs to output file
        os.makedirs(os.path.dirname(output_pairs_lst_file), exist_ok=True)
        with open(output_pairs_lst_file, 'w') as f:
            for pair in pairs:
                f.write(f"{str(pair)}\n")

        return pairs


def GeneratePairMap(moles: List[Chem.Mol, ], output_pairs_file) -> List[Pair]:
    import tempfile
    origin_cwd = os.getcwd()
    use_moles = copy.deepcopy(moles)

    with tempfile.TemporaryDirectory() as tmp_dir:
        os.chdir(tmp_dir)
        with working_dir('mol_files'):
            for _mol in use_moles:

                _update_mol = copy.deepcopy(_mol)
                # _update_mol = Chem.RemoveAllHs(_update_mol)
                # _update_mol = Chem.AddHs(_update_mol)

                # Generate 3D conformation
                rdDistGeom.EmbedMolecule(_update_mol)

                # Optimize geometry using UFF
                rdForceFieldHelpers.UFFOptimizeMolecule(_update_mol)

                # Further optimize using MMFF94 if parameters are available
                if rdForceFieldHelpers.MMFFHasAllMoleculeParams(_update_mol):
                    rdForceFieldHelpers.MMFFOptimizeMolecule(_update_mol)

                # Mol_2_mol2.write_mol2_file(_update_mol, f'{_update_mol.GetProp("_Name")}.mol2')
                Chem.MolToMolFile(_update_mol, f'{_update_mol.GetProp("_Name")}.sdf')
        _lomap_path = 'mol_files'
        _lomap_api = Lomap_api(_lomap_path)
        pairs = _lomap_api.output_pairs_lst(output_pairs_file)
        for i in range(0, len(pairs), 5):
            row = pairs[i:i+5]
            formatted_row = [f"{str(p):^10}" for p in row]
            print(" ".join(formatted_row), flush=True)
        os.chdir(origin_cwd)
        return pairs


# Example usage
if __name__ == "__main__":
    smi_files = os.listdir()
    for file in smi_files:
        if file.endswith('.smi'):
            _input_file = os.path.join(os.getcwd(), file)
            # print(_input_file)
            _output_file = file.replace('.smi', '.mol2')
            with working_dir('mol2_files'):
                ts = Mol_2_mol2(_input_file)
                mol = ts._3d_mol
                ts.write_mol2_file(mol, _output_file)
    # Test Lomap API
    lomap_path = 'mol2_files'
    # known_actives_file = 'known_mol.lst'
    # known_actives_file = os.path.join(os.getcwd(), known_actives_file)
    # max_dist_from_actives = 2
    # lomap_api = Lomap_api(lomap_path, known_actives_file=known_actives_file, max_dist_from_actives=max_dist_from_actives)
    lomap_api = Lomap_api(lomap_path)
    lomap_api.output_pairs_lst()
