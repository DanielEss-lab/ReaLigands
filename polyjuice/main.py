import time

from openbabel import openbabel
from openbabel import pybel
import glob
import LigandExtract

def find_nitrogen_multi_bond(atom):
    num_bonds = 0
    for _ in openbabel.OBAtomBondIter(atom.OBAtom):
        num_bonds += 1

    if num_bonds < 3:
        return True
    return False


class Finder:
    def __init__(self, mol):
        self.multi_bond = False
        self.metal_hit = 0
        self.mol = mol

        # If using tmQm dataset, uncomment these lines. Bonds to the metal 
        # are not always correct in Openbabel, so the correct number of bonds is 
        # grabbed from the comment line of file and missing bonds are manually created.
        # This step is not necessary for the tmQmg dataset. 
        # self.num_metal_bonds = int(re.findall("MND = (\d+)", str(mol))[0])
        # self.find_nearest_atoms()

        # During testing we found that molecules containing Si were often missing bonds.
        # This issue was corrected by manually adding in missing bonds for the closest
        # non-H atoms for each Si atom.
        self.adjust_si_bonds()

        # Resetting bond ID values because they are altered in the find_nearest_atoms() and 
        # adjust_si_bonds()
        self.set_bond_id()
    
    def find_ligands(self, atom):
        a = atom.OBAtom

        for bond in openbabel.OBAtomBondIter(a):
            if bond.GetId() == 0:
                bond.SetId(1)

                # If the atom on the other side of the bond is H or metal we will
                # not traverse across bond. If it is a metal, we add to a counter in order to 
                # classify ligand denticity. Otherwise, find_ligands is recursively called to 
                # continue ligand traversal
                if self.mol.atoms[bond.GetNbrAtomIdx(a) - 1].atomicnum == 1:
                    continue
                elif self.mol.atoms[bond.GetNbrAtomIdx(a) - 1].OBAtom.IsMetal():
                    self.metal_hit += 1
                    continue
                else:
                    self.find_ligands(self.mol.atoms[bond.GetNbrAtomIdx(a) - 1])

    def start(self, metal):
        bond_iter = 0
        total_metal_hits = 0
        m = metal.OBAtom
        ligand_count = 1

        # bond IDs are used to track visited atoms/bonds here and in find_ligands. Bond IDs are 
        # initialized to 0 and set to 1 upon visiting. This conditional helps in ligand traversal itself
        # as well as in preventing the visiting of polydentate ligands multiple times.

        for bond in openbabel.OBAtomBondIter(m):
            if bond.GetId() == 0:
                bond.SetId(1)

                # bond.GetNmbAtomIdx takes in a specified atom and returns its neighbor
                # in this case we pass in our metal base and we get back the 
                # atom on the other side which serves as a start point to the ligand
                ligand_start = self.mol.atoms[bond.GetNbrAtomIdx(m) - 1]

                # Nitrogen bond special case
                if ligand_start.atomicnum == 7:
                    if find_nitrogen_multi_bond(ligand_start):
                        self.metal_hit += 1  # If the start of the ligand is a nitrogen with a double or triple bond
                        # the metal, then add one to the metal count

                if not do_ligand_filter:
                    # identify the rest of the ligand, using our identified start point
                    self.find_ligands(ligand_start)
                elif do_ligand_filter and ligand_start.atomicnum == ligand_filter:
                    self.find_ligands(ligand_start)
                else:
                    bond_iter += 1
                    continue

                total_metal_hits += self.metal_hit
                # A copy of the molecule is passed to a new LigandExtractor instance as extracting 
                # a ligand in place causes issues with atom indices and other features.
                if self.metal_hit == 0:
                    copy_molecule = pybel.Molecule(openbabel.OBMol(self.mol.OBMol))
                    sub = LigandExtract.LigandExtractor(copy_molecule, bond_iter, mol_num, 1, m.GetAtomicNum(), ligand_count, self.csd_id)
                    sub.extract_ligand()
                else:
                    copy_molecule = pybel.Molecule(openbabel.OBMol(self.mol.OBMol))
                    sub = LigandExtract.LigandExtractor(copy_molecule, bond_iter, mol_num, self.metal_hit+1, m.GetAtomicNum(), ligand_count, self.csd_id)
                    sub.extract_ligand()
                    self.metal_hit = 0
                
                ligand_count += 1

            bond_iter += 1
             
    def adjust_si_bonds(self):
        si_count = 0
        si_bonded_lists = []

        # Marking all Si atoms in mol and calulcating distances from current Si atom to all other atoms
        for atom in self.mol:
            if atom.atomicnum == 14:
                atom.OBAtom.SetId(si_count)
                distance_list = []
                for otherAtom in self.mol:
                    if not otherAtom.atomicnum == 14:
                        distance_list.append((otherAtom, atom.OBAtom.GetDistance(otherAtom.OBAtom)))
                si_bonded_lists.append(distance_list)
                si_count += 1

        # For each Si atom we found, we examine the corresponding sorted distance list.
        # If bonds are not found to nearby non-H atoms, bonds are manually inserted.
        for i in range(si_count):
            current_list = sorted(si_bonded_lists[i], key=lambda x: x[1])
            for atom in self.mol:
                if atom.atomicnum == 14 and atom.OBAtom.GetId() == i:
                    for j in range(4):
                        if self.mol.OBMol.GetBond(atom.OBAtom, current_list[j][0].OBAtom) is None:
                            if current_list[j][0].atomicnum == 1:
                                j -= 1
                            else:
                                bond = openbabel.OBBond()
                                bond.SetBegin(atom.OBAtom)
                                bond.SetEnd(current_list[j][0].OBAtom)
                                self.mol.OBMol.AddBond(bond)
                        else:
                            pass
                    break

    def find_nearest_atoms(self):
        metal_bonded_list = []
        for atom in self.mol:
            if atom.OBAtom.IsMetal():
                for otherAtom in self.mol:
                    if not otherAtom.OBAtom.IsMetal():
                        metal_bonded_list.append((otherAtom, atom.OBAtom.GetDistance(otherAtom.OBAtom)))

        new_metal_bonded_list = sorted(metal_bonded_list, key=lambda x: x[1])

        for atom in self.mol:
            if atom.OBAtom.IsMetal():
                for j in range(self.num_metal_bonds):
                    if self.mol.OBMol.GetBond(atom.OBAtom, new_metal_bonded_list[j][0].OBAtom) is None:
                        if new_metal_bonded_list[j][0].atomicnum == 1:
                            i = 0
                            for _ in openbabel.OBAtomBondIter(new_metal_bonded_list[j][0].OBAtom):
                                i += 1
                            if i > 0:
                                break
                            else:
                                bond = openbabel.OBBond()
                                bond.SetBegin(atom.OBAtom)
                                bond.SetEnd(new_metal_bonded_list[j][0].OBAtom)
                                self.mol.OBMol.AddBond(bond)
                        else:
                            bond = openbabel.OBBond()
                            bond.SetBegin(atom.OBAtom)
                            bond.SetEnd(new_metal_bonded_list[j][0].OBAtom)
                            self.mol.OBMol.AddBond(bond)
                    
    def set_bond_id(self):
        for atom in self.mol:
            a = atom.OBAtom
            for bond in openbabel.OBAtomBondIter(a):
                bond.SetId(0)

    def set_csd_id(self, id):
        self.csd_id = id


filter_input = input("Do you want to only identify ligands for a certain metal? ").lower()
do_metal_filter = filter_input == 'y' or filter_input == 'yes'
if do_metal_filter:
    metal_filter = int(input("Please enter the atomic number of the metal you want: "))

filter_input = input("Do you want to only identify ligands with a certain starting atom? ").lower()
do_ligand_filter = filter_input == 'y' or filter_input == 'yes'
if do_ligand_filter:
    ligand_filter = int(input("Please enter the atomic number of the ligand's starting atom: "))

start_time = time.time()
mol_num = 0
# Current file path looks at all xyz files in current directory
for file in glob.glob("./*.xyz"):
    for molecule in pybel.readfile("xyz", file):
        mol_num += 1
        finder = Finder(molecule)

        # tmQmg files are all named based on CSD ID. We want to preserve that information 
        # in order to match ligands to original structures. These lines will need to be modified
        # based on current working directory and input files. 
        name = file.lstrip('/home/')
        finder.set_csd_id(name.rstrip('.xyz'))

        # look for metal atom to start search
        for metal_atom in molecule:
            if metal_atom.OBAtom.IsMetal() and not do_metal_filter:
                finder.start(metal_atom)
            elif metal_atom.OBAtom.IsMetal() and metal_atom.atomicnum == metal_filter:
                finder.start(metal_atom)
            else:
                break

end_time = time.time()
elapsed_time = end_time - start_time
minutes = elapsed_time // 60
seconds = elapsed_time % 60
print(f"Ran through {mol_num} molecules in {minutes:.0f} minutes and {seconds:.2f} seconds.")