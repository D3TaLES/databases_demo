import re
from pymatgen.io.gaussian import GaussianOutput


class ParseGausLog:
    """
    Class to process Gaussian logfiles.
    Copyright 2021, University of Kentucky
    """

    def __init__(self, filepath, ):
        self.log_path = filepath
        self.code_used = 'Gaussian'
        self.pmgmol = GaussianOutput(filepath)

        self.functional = self.pmgmol.functional
        self.basis_set = self.pmgmol.basis_set
        self.charge = self.pmgmol.charge
        self.spin_multiplicity = self.pmgmol.spin_multiplicity
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

        self.homo = self.homo_lumo["homo"]
        self.lumo = self.homo_lumo["lumo"]

    @property
    def functional(self):
        with open(self.log_path) as f:
            values = [line.strip().split()[3] for line in f if re.search(r"{}".format('Composite Functional:'), line)]
        return values[0]

    @property
    def basis_set(self):
        with open(self.log_path) as f:
            values = [line.strip().split()[2] for line in f if re.search(r"{}".format('Basis Set:'), line)]
        return values[0]

    @property
    def charge(self):
        with open(self.log_path) as f:
            values = [line.strip().split()[5] for line in f if re.search(r"{}".format('charge'), line)]
        target_value = values[0]
        return int(re.sub("[^0-9]", "", target_value))

    @property
    def spin_multiplicity(self):
        with open(self.log_path) as f:
            values = [line.strip().split()[8] for line in f if re.search(r"{}".format('multiplicity'), line)]
        target_value = values[0]
        return int(re.sub("[^0-9]", "", target_value))

    @property
    def scf_total_energy(self):
        with open(self.log_path) as f:
            values = [line.strip().split()[3] for line in f if re.search(r"{}".format('Final Energy:'), line)]
        target_value = values[0]
        return float(target_value) * 27.2114  # convert to eV

    @property
    def homo_lumo(self):
        """
        Get homo and lumo energies from a Gaussian molecule
        :return: [homo, lumo] - list containing homo then lumo in eV
        """
        return {"homo": 0, "lumo": 0}
        num_electrons = None
        with open(self.log_path) as f:
            for line in f.readlines():
                if re.search(r"{}".format('Electrons'), line) and not num_electrons:
                    num_electrons = int(line.strip().split()[3])
            eigenvalues = {}

        eigens = list(eigenvalues.values())
        homo = eigens[num_electrons - 1] * 27.2114  # convert to eV
        lumo = eigens[num_electrons] * 27.2114  # convert to eV

        return {"homo": homo, "lumo": lumo}

