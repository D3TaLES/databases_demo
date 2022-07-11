import re
from pymatgen.io.gaussian import GaussianOutput


class ParseGausLog:
    """
    Class to process Gaussian logfiles.
    Copyright 2021, University of Kentucky
    Args:
        filepath (str) : filepath to data file
    """

    def __init__(self, filepath):
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
        self.singlet_excitations = self.get_tddft_excitations(filepath)["Singlet"]
        self.singlet_excitations_energy = [e[0] for e in self.singlet_excitations]  # extract eV energy value

    @staticmethod
    def homo_lumo(pymatgen_mol):
        """
        Get homo and lumo energies from a Gaussian molecule
        return
            homo_lumo (dict) : dictionary containing homo then lumo in eV
        """
        num_electrons = pymatgen_mol.electrons[0]
        eigens = list(pymatgen_mol.eigenvalues.values())[0]
        homo = eigens[num_electrons - 1] * 27.2114  # convert to eV
        lumo = eigens[num_electrons] * 27.2114  # convert to eV

        return {"homo": homo, "lumo": lumo}

    @staticmethod
    def get_tddft_excitations(log_path):
        """
        Read excitation energies after a TD-DFT calculation.
        Returns:
            A list: A list of tuple for each transition such as
                    [(energy (eV), lambda (nm), oscillatory strength), ... ]
        """

        float_patt = re.compile(r"\s*([+-]?\d+\.\d+)")
        state_patt = re.compile(r"[a-zA-Z]*let")
        transitions = {"Singlet": [],
                       "Doublet": [],
                       "Triplet": [],
                       }

        # read in file
        with open(log_path, "r") as f:
            line = f.readline()
            td = False
            while line != "":
                if re.search(r"^\sExcitation energies and oscillator strengths:", line):
                    td = True

                if td:
                    if re.search(r"^\sExcited State\s*\d", line):
                        val = [float(v) for v in float_patt.findall(line)]
                        try:
                            state = state_patt.findall(line)[0]
                        except Exception:
                            state_val = val.pop(0)
                            if round(state_val) == 1:
                                state = "Singlet"
                            elif round(state_val) == 2:
                                state = "Doublet"
                            elif round(state_val) == 3:
                                state = "Triplet"
                            else:
                                raise ValueError(
                                    "Calculation has a spin state greater than triplet -- spin {}".format(state_val))
                        transitions[state].append(tuple(val))
                line = f.readline()
        return transitions


class ParsePsi4Log:
    """
    Class to process Gaussian logfiles.
    Copyright 2021, University of Kentucky
    Args:
        filepath (str) : filepath to data file
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
        return
            homo_lumo (dict) : dictionary containing homo then lumo in eV
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

