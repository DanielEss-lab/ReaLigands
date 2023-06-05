from openbabel import openbabel
import BondLengthSwitch


def count_bonds(a):
    bond_count = 0
    for _ in openbabel.OBAtomBondIter(a):
        bond_count += 1
    return bond_count


def find_water_like(oxygen):
    o = oxygen.OBAtom
    oxygen_bond_count = count_bonds(o)
    if oxygen_bond_count > 2 and (oxygen.atomicnum == 8 or oxygen.atomicnum == 16):
        return True
    return False


def find_ammonia(nitrogen):
    n = nitrogen.OBAtom
    nitrogen_bond_count = count_bonds(n)
    if nitrogen_bond_count > 3 and (nitrogen.atomicnum == 7 or nitrogen.atomicnum == 15):
        return True
    return False


def find_oxygen_double_bond(oxygen):
    o = oxygen.OBAtom
    bonds = count_bonds(o)
    if bonds == 1 and oxygen.atomicnum == 8:
        return True
    return False


class LigandChargeFinder:
    def __init__(self, molecule):
        self.molecule = molecule

    def find_co(self, carbon):
        c = carbon.OBAtom
        carbon_bond_count = count_bonds(c)
        for bond in openbabel.OBAtomBondIter(c):
            atomic_num = self.molecule.atoms[bond.GetNbrAtomIdx(c) - 1].atomicnum
            if atomic_num != 8 and atomic_num != 16:
                continue
            if carbon_bond_count == 2:
                oxygen = self.molecule.atoms[bond.GetNbrAtomIdx(c) - 1].OBAtom
                oxygen_bond_count = count_bonds(oxygen)
                if oxygen_bond_count == 1:
                    return True
                return False
            return False
        return False

    def find_cn(self, carbon):
        c = carbon.OBAtom
        carbon_bond_count = count_bonds(c)
        for bond in openbabel.OBAtomBondIter(c):
            atom_x = self.molecule.atoms[bond.GetNbrAtomIdx(c) - 1].atomicnum
            if atom_x != 7:
                continue
            if carbon_bond_count == 2:
                nitrogen = self.molecule.atoms[bond.GetNbrAtomIdx(c) - 1].OBAtom
                nitrogen_bond_count = count_bonds(nitrogen)
                if nitrogen_bond_count > 1:
                    return True
                return False
            return False
        return False

    def find_carbon_double_bond(self, carbon):
        c = carbon.OBAtom
        bonds = count_bonds(c)
        if bonds == 3 and carbon.atomicnum == 6:
            for bond in openbabel.OBAtomBondIter(c):
                if self.molecule.atoms[bond.GetNbrAtomIdx(c) - 1].OBAtom.IsMetal():
                    return BondLengthSwitch.switch(self.molecule.atoms[bond.GetNbrAtomIdx(c) - 1].atomicnum, carbon.atomicnum, bond.GetLength())
        if bonds == 2 and carbon.atomicnum == 6:
            if self.find_co(carbon) or self.find_cn(carbon):
                return False
            return True
        return False

    def help_find_aromatic_ring(self, atom, iteration):
        num_bonds = count_bonds(atom.OBAtom)
        if (atom.atomicnum == 6 and num_bonds != 3) or (atom.atomicnum == 7 and num_bonds > 3):
            return False
        if iteration > 6:
            return False
        a = atom.OBAtom
        is_ring = False
        for bond in openbabel.OBAtomBondIter(a):
            other_end = self.molecule.atoms[bond.GetNbrAtomIdx(a) - 1]
            if other_end.OBAtom.GetId() == 3 and bond.GetId() == 1:
                bond.SetId(2)
                return True
            elif bond.GetId() == 1:
                bond.SetId(2)
                is_ring = self.help_find_aromatic_ring(other_end, iteration + 1)
            if is_ring:
                break

        return is_ring

    def find_aromatic_ring(self, atom):
        if atom.atomicnum != 6 and atom.atomicnum != 7:     # The aromatic ring should start with either a nitrogen or a carbon
            return False
        num_bonds = count_bonds(atom.OBAtom)
        if (atom.atomicnum == 6 and num_bonds != 3) or (atom.atomicnum == 7 and num_bonds > 3):
            return False
        for bond in openbabel.OBMolBondIter(self.molecule.OBMol):
            bond.SetId(1)
        for element in self.molecule:
            element.OBAtom.SetId(1)
        atom.OBAtom.SetId(3)
        a = atom.OBAtom
        is_ring = False
        for bond in openbabel.OBAtomBondIter(a):
            if self.molecule.atoms[bond.GetNbrAtomIdx(a) - 1].OBAtom.IsMetal():
                bond.SetId(2)
            elif bond.GetId() == 1:
                bond.SetId(2)
                is_ring = self.help_find_aromatic_ring(self.molecule.atoms[bond.GetNbrAtomIdx(a) - 1], 1)
            if is_ring:
                return True

        return is_ring

    def find_aromatic_nitrogen(self, atom):
        if atom.atomicnum != 7:
            return False
        is_ring = self.find_aromatic_ring(atom)
        if not is_ring:
            return False

        return True

    def find_carbon_two_nitrogens(self, atom):
        if atom.atomicnum != 6:
            return False
        a = atom.OBAtom
        if count_bonds(a) != 3:
            return False
        nitrogen_count = 0
        for bond in openbabel.OBAtomBondIter(a):
            if self.molecule.atoms[bond.GetNbrAtomIdx(a) - 1].atomicnum == 7:
                nitrogen_count += 1

        if nitrogen_count == 2:
            nitrogen_1_bonds = 0
            nitrogen_2_bonds = 0
            for bond in openbabel.OBAtomBondIter(a):
                if self.molecule.atoms[bond.GetNbrAtomIdx(a) - 1].atomicnum == 7:
                    if nitrogen_1_bonds == 0:
                        nitrogen_1_bonds = count_bonds(self.molecule.atoms[bond.GetNbrAtomIdx(a) - 1].OBAtom)
                    else:
                        nitrogen_2_bonds = count_bonds(self.molecule.atoms[bond.GetNbrAtomIdx(a) - 1].OBAtom)
            if nitrogen_1_bonds == 3 and nitrogen_2_bonds == 3:
                return True

        return False

    def change_charge(self, atom):
        neutral_ligand = self.find_cn(atom) or self.find_co(atom) or find_water_like(atom) or \
                         find_ammonia(atom) or self.find_aromatic_nitrogen(atom) or \
                         self.find_carbon_two_nitrogens(atom)
        double_bond = find_oxygen_double_bond(atom) or self.find_carbon_double_bond(atom)
        triple_bond = False
        if neutral_ligand:
            return 0
        elif double_bond:
            return -2
        elif triple_bond:
            return -3
        else:                   # Negative ligands
            return -1
