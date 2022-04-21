import re
from pymatgen.io.gaussian import GaussianOutput


class ParseGausLog:
    """
    Class to process Gaussian logfiles.
    Copyright 2021, University of Kentucky
    """

    def __init__(self, filepath, ):
        self.log_path = filepath
        self.pmgmol = GaussianOutput(filepath)
        self.code_used = 'Gaussian'

        self.parse_file()

    def parse_file(self):
        self.functional = self.pmgmol.functional
        self.basis_set = self.pmgmol.basis_set
        self.charge = self.pmgmol.charge
        self.spin_multiplicity = self.pmgmol.spin_multiplicity
        self.geometry = self.pmgmol.final_structure.as_dict()['sites']
        self.scf_total_energy = self.pmgmol.final_energy * 27.2114  # convert to eV
        self.homo = self.homo_lumo(self.pmgmol)["homo"]
        self.lumo = self.homo_lumo(self.pmgmol)["lumo"]

    @staticmethod
    def homo_lumo(pymatgen_mol):
        """
        Get homo and lumo energies from a Gaussian molecule
        :return: [homo, lumo] - list containing homo then lumo in eV
        """
        num_electrons = pymatgen_mol.electrons[0]
        eigens = list(pymatgen_mol.eigenvalues.values())[0]
        homo = eigens[num_electrons - 1] * 27.2114  # convert to eV
        lumo = eigens[num_electrons] * 27.2114  # convert to eV

        return {"homo": homo, "lumo": lumo}


class ParsePsi4Log:
    """
    Class to process Gaussian logfiles.
    Copyright 2021, University of Kentucky
    """

    def __init__(self, filepath, ):
        self.log_path = filepath
        self.code_used = 'Psi4'

        self.parse_file()

    def parse_file(self):
        self.functional = self.find_in_file('Composite Functional:', 3)
        self.basis_set = self.find_in_file('Basis Set:', 2)
        self.charge = int(re.sub("[^0-9]", "", self.find_in_file('charge', 5)))
        self.spin_multiplicity = int(re.sub("[^0-9]", "", self.find_in_file('multiplicity', 8)))
        self.geometry = []
        self.scf_total_energy = float(self.find_in_file('Final Energy:', 3)) * 27.2114  # convert to eV

        self.homo = self.homo_lumo["homo"]
        self.lumo = self.homo_lumo["lumo"]

    @property
    def homo_lumo(self):
        """
        Get homo and lumo energies from a Gaussian molecule
        :return: [homo, lumo] - list containing homo then lumo in eV
        """
        return {"homo": 0, "lumo": 0}
        with open(self.log_path) as f:
            num_electrons: float = self.find_in_file('Electrons', 3)[0]
            eigenvalues = {}

        eigens = eigenvalues.values()
        homo = eigens[num_electrons - 1] * 27.2114  # convert to eV
        lumo = eigens[num_electrons] * 27.2114  # convert to eV

        return {"homo": homo, "lumo": lumo}

    def find_in_file(self, search_text, target_position):
        with open(self.log_path) as f:
            value = [line.strip().split()[target_position] for line in f if re.search(r"{}".format(search_text), line)]
        return value[0]
