import os
from openbabel import openbabel

import HeaderReplacer
import LigandChargeFinder

dentate_map = {1:"Monodentate", 2:"Bidentate", 3:"Tridentate"}
atomic_nums_to_elem = {1:"H", 2:"He", 3:"Li", 4:"Be", 5:"B", 6:"C", 7:"N", 8:"O", 9:"F", 11:"Na",
                       12:"Mg", 13:"Al", 14:"Si", 15:"P", 16:"S", 17:"Cl", 19:"K", 20:"Ca", 21:"Sc",  
                       22:"Ti", 23:"V", 24:"Cr", 25:"Mn", 26:"Fe", 27:"Co", 28:"Ni", 29:"Cu", 
                       30:"Zn", 31:"Ga", 32:"Ge", 33:"As", 34:"Se", 35:"Br", 37:"Rb", 38:"Sr", 
                       39:"Y", 40:"Zr", 41:"Nb", 42:"Mo", 43:"Tc", 44:"Ru", 45:"Rh", 46:"Pd", 
                       47:"Ag", 48:"Cd", 49:"In", 50:"Sn", 51:"Sb", 52:"Te", 53:"I", 55:"Cs", 56:"Ba",
                       57:"La", 72:"Hf", 73:"Ta", 74:"W", 75:"Re", 76:"Os", 77:"Ir", 78:"Pt", 79:"Au",
                       80:"Hg", 104:"Rf", 105:"Db", 106:"Sg", 107:"Bh", 108:"Hs"}

class LigandExtractor:
    def __init__(self, mol, bond_iter, num_atom, denticity, center_metal, ligand_num, csd_id):
        self.mol = mol
        self.bond_iter = bond_iter
        self.num_atom = num_atom
        self.metal_ind = 0
        self.new_carbon_ind = 0
        self.start_index = 0

        for atom in self.mol:
            atom.OBAtom.SetId(0)
            if atom.OBAtom.IsMetal():
                self.type = atom.type

        self.set_bond_id()
        self.new_charge = self.modify_charge()
        self.set_bond_id()

        self.denticity = denticity
        self.atom_connections = []
        self.connecting_indices = []
        self.center_metal = center_metal
        self.ligand_num = ligand_num
        self.csd_id = csd_id

    def modify_charge(self):
        for atom in self.mol:
            if atom.OBAtom.IsMetal():
                counter = 0
                for bond in openbabel.OBAtomBondIter(atom.OBAtom):
                    if self.bond_iter == counter:
                        ligand_atom = self.mol.atoms[bond.GetNbrAtomIdx(atom.OBAtom) - 1]
                        charge_finder = LigandChargeFinder.LigandChargeFinder(self.mol)
                        new_charge = charge_finder.change_charge(ligand_atom)

                        return new_charge

                    counter += 1

    def set_bond_id(self):
        for atom in self.mol:
            a = atom.OBAtom
            for bond in openbabel.OBAtomBondIter(a):
                bond.SetId(0)

    def convert_connections(self):
        self.atom_connections.sort()
        output = ""
        for i in range(len(self.atom_connections)):
            output += atomic_nums_to_elem[self.atom_connections[i]] 
        return output

    def find_ligands(self, atom):
        # Atoms are typically given an ID value of 2 but we set values to 3 
        # and save the atomic number if it is one of the connecting atoms to
        # the metal
        a = atom.OBAtom
        a.SetId(2)
        for bond in openbabel.OBAtomBondIter(a):
            if bond.GetId() == 0:
                bond.SetId(1)

                if self.mol.atoms[bond.GetNbrAtomIdx(a) - 1].OBAtom.IsMetal():
                    if a.GetId() < 3:
                        a.SetId(3)
                        self.atom_connections.append(a.GetAtomicNum())
                    
                    continue
                else:
                    self.find_ligands(self.mol.atoms[bond.GetNbrAtomIdx(a) - 1])
            else:
                pass
            
    def getCSDid(self, file_path):
        with open(file_path, 'r') as file:
            content = file.read()
            index = content.find('CSD_code = ') + 11 
            csd_id = content[index:index+6]
            return csd_id

    def extract_ligand(self):
        for atom in self.mol:
            if atom.OBAtom.IsMetal():
                counter = 0
                for bond in openbabel.OBAtomBondIter(atom.OBAtom):
                    bond.SetId(1)
                    # The bond_iter value passed in from main.py indicates which bond we want to start our extraction process on.
                    if self.bond_iter == counter:
                        self.find_ligands(self.mol.atoms[bond.GetNbrAtomIdx(atom.OBAtom) - 1])
                        # Sets the first atom in the ligand's ID to one so it doesn't get deleted
                        self.mol.atoms[bond.GetNbrAtomIdx(atom.OBAtom) - 1].OBAtom.SetId(3)
                        self.atom_connections.append(self.mol.atoms[bond.GetNbrAtomIdx(atom.OBAtom) - 1].atomicnum)
                    
                    counter += 1

        for atom in self.mol:
            if atom.OBAtom.GetId() < 2:
                self.mol.OBMol.DeleteAtom(atom.OBAtom)

        # Indices are changed after atom deletions in the previous loop, so this loop 
        # iterates over the mol again to grab the updated indices of the connecting atoms
        for atom in self.mol:
            if atom.OBAtom.GetId() == 3:
                self.connecting_indices.append(atom.OBAtom.GetIndex() + 1)

        connecting_atoms = self.convert_connections()

        # creating directories to write ligands to
        base_path = 'ligandstmqmg'
        if not os.path.exists(f'{base_path}'):
            os.makedirs(f'{base_path}')

        if not self.denticity in dentate_map.keys():
            self.dentate_name = "TetrasAndGreater"
        else:
            self.dentate_name = dentate_map[self.denticity]

        if not os.path.exists(f'{base_path}/{self.dentate_name}'):
            os.makedirs(f'{base_path}/{self.dentate_name}')
        
        if not os.path.exists(f'{base_path}/{self.dentate_name}/{connecting_atoms}'):
            os.makedirs(f'{base_path}/{self.dentate_name}/{connecting_atoms}')

        self.mol.OBMol.SetTotalCharge(self.new_charge)

        mol_file_path = f"{base_path}/{self.dentate_name}/{connecting_atoms}/{self.type}{self.num_atom}-{self.bond_iter}.mol2"

        self.mol.write("mol2", mol_file_path, True)

        # File name is rewritten based on CSD ID, metal core, which ligand num 
        # it is from the original mol, num of connecting atoms, and connecting atoms.

        # If using original tmQm data, uncomment this line to get CSD ID from file. 
        # csd_id = self.getCSDid(mol_file_path)

        center_metal_symbol = atomic_nums_to_elem[self.center_metal]
        ligand_num = str(self.ligand_num)
        coordination_num = str(len(self.atom_connections))
        new_file_name = f'./{base_path}/{self.dentate_name}/{connecting_atoms}/'
        new_file_name += self.csd_id + '_' + center_metal_symbol + '_' + \
                        ligand_num + '_' + coordination_num + '_'+ connecting_atoms + '.mol2'
        
        os.rename(mol_file_path, new_file_name)

        # Replacing first line of header to save the connecting indices of the atoms
        header_replacer = HeaderReplacer.HeaderReplacer(new_file_name)
        header_replacer.replace_header(self.connecting_indices)
