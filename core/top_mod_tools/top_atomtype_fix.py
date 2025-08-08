import os
import parmed
# from parmed.constants import PrmtopPointers
import logging
import numpy as np

logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)


def _del_in_ICO(remove_atomtype_idx: int, input_list):
    """
    remove_numtype = 3
    a a | b
    a a | b
    - - + -
    c c | d

    to

    a a b
    a a b
    c c d
    :param remove_atomtype_idx: the idx of numtype to remove
    :return:
    """

    numtype_nums = int(np.sqrt(len(input_list)))
    return_list = []
    for i in range(numtype_nums):
        if i == remove_atomtype_idx - 1:
            continue
        return_list.extend([input_list[i * numtype_nums + j]
                            for j in range(numtype_nums) if j != remove_atomtype_idx - 1])

    return return_list


def COEF_reindex(ico_list, coef_list, ntypes):
    new_coef_list = [0.0 for _ in range(sum(range(ntypes + 1)))]

    def _order_num(i, j):
        return sum(range(0, j + 1)) + i + 1

    for i in range(ntypes):
        for j in range(ntypes):
            if i > j:
                continue

            old_idx = ico_list[i * ntypes + j] - 1
            if old_idx == -1:
                continue
            new_coef_list[_order_num(i, j) - 1] = coef_list[old_idx]

    new_ico_list = []
    for i in range(ntypes):
        for j in range(ntypes):
            if i > j:
                new_ico_list.append(new_ico_list[j * ntypes + i])
                continue

            old_idx = ico_list[i * ntypes + j]
            if old_idx == -1:
                new_ico_list.append(-1)
                continue

            new_ico_list.append(_order_num(i, j))
    return new_ico_list, new_coef_list


def top_atomtype_fix(top_file, output_file=None):
    amb_p = parmed.amber.AmberParm(top_file)
    if output_file is None:
        output_file = os.path.abspath(top_file)

    old_atomtype_dict = {}

    for atom in amb_p.atoms:
        old_atomtype_dict[atom.nb_idx] = atom.type

    sorted_atomtype = sorted(old_atomtype_dict.items(), key=lambda x: x[0], reverse=False)
    logger.debug(f'NonBondIdx: {sorted_atomtype}')
    new_atomtype_map = {}

    if len(sorted_atomtype) < sorted_atomtype[-1][0]:
        logger.info('Non-continuous AtomType idxes')
        unoccupied_atomtype = set(range(1, sorted_atomtype[-1][0] + 1)) - set([i[0] for i in sorted_atomtype])
        unoccupied_atomtype_sorted = sorted(unoccupied_atomtype, reverse=True)
        logger.info(f"Discard {unoccupied_atomtype_sorted}")

        new_ico_list = amb_p.parm_data['NONBONDED_PARM_INDEX']
        for discard_numtype in unoccupied_atomtype_sorted:
            new_ico_list = _del_in_ICO(discard_numtype, new_ico_list)

        nico, nAcoef = COEF_reindex(new_ico_list,
                                     amb_p.parm_data['LENNARD_JONES_ACOEF'],
                                     len(sorted_atomtype))
        nico, nBcoef = COEF_reindex(new_ico_list,
                                     amb_p.parm_data['LENNARD_JONES_BCOEF'],
                                     len(sorted_atomtype))

        amb_p.parm_data['NONBONDED_PARM_INDEX'] = nico
        amb_p.parm_data['LENNARD_JONES_ACOEF'] = nAcoef
        amb_p.parm_data['LENNARD_JONES_BCOEF'] = nBcoef
        # PrmtopPointers.NTYPES = 1
        amb_p.parm_data['POINTERS'][1] = len(sorted_atomtype)
        amb_p.load_pointers()

        for i in range(len(sorted_atomtype)):
            new_atomtype_map[sorted_atomtype[i][0]] = i + 1

        for atom in amb_p.atoms:
            atom.nb_idx = new_atomtype_map[atom.nb_idx]

        amb_p.remake_parm()
        amb_p.load_atom_info()
        amb_p.load_structure()
        amb_p.save(output_file, format='amber', overwrite=True)


if __name__ == '__main__':
    pass
