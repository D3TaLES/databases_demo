import re
import json
import uuid
import pandas as pd
from datetime import date

import pubchempy as pcp
from scipy.signal import find_peaks
from rdkit.Chem.rdchem import Mol
from rdkit.Chem.rdmolops import AddHs
from pymatgen.io.gaussian import GaussianOutput
from rdkit.Chem import MolFromSmiles, MolToSmiles, AllChem
from rdkit.Chem.rdMolDescriptors import CalcMolFormula, CalcExactMolWt


class GenerateMolInfo:
    """
    Generate json object for insertion from smiles string
    Copyright 2021, University of Kentucky
    Args:
        smiles (str) : smiles string
        source (str) : which group the molecule comes from
        names (list) : list of names for molecule
    """

    def __init__(self, smiles, source="", date_made="", names=None, sql=True):
        smiles = smiles
        names = names or []
        # Generate rdkit mol and final (cleaned) smiles
        rdkmol = MolFromSmiles(smiles)
        clean_smile = MolToSmiles(rdkmol)
        rdkmol_hs = AddHs(rdkmol)
        AllChem.EmbedMolecule(rdkmol_hs)
        pcpmol = pcp.get_compounds(clean_smile, namespace="smiles")[0]

        # Populate class
        self.sql = sql
        self._id = str(pcpmol.iupac_name)
        self.source = source
        self.date_made = date_made or date.today()
        self.smiles = clean_smile
        self.molecular_formula = CalcMolFormula(rdkmol)
        self.number_of_atoms = Mol.GetNumAtoms(rdkmol)
        self.molecular_weight = CalcExactMolWt(rdkmol)
        try:
            self.synonyms = names + pcpmol.synonyms
        except TypeError:
            self.synonyms = names

    @property
    def data(self):
        """
        Returns molecule information in a dictionary that matches the No-SQL or SQl schema
        """
        data_dict = {
            "smiles": self.smiles,
            "molecular_formula": self.molecular_formula,
            "number_of_atoms": self.number_of_atoms,
            "molecular_weight": self.molecular_weight,
            "source": self.source,
            "date_made": self.date_made,
        }
        if self.sql:
            data_dict.update({"mol_id": self._id})
        else:
            data_dict.update({"_id": self._id, "synonyms": self.synonyms})

        json_data = json.dumps(data_dict, default=str)
        return json.loads(json_data)

    @property
    def synonym_data(self):
        """
        Returns synonym information in a dictionary that matches the SQL Synonyms Table schema
        """
        return [{"synonym": synonym, "mol_id": self._id} for synonym in self.synonyms]


class ProcessDFT:
    """
    Class to process DFT Gaussian logfiles.
    Copyright 2021, University of Kentucky
    Args:
        filepath (str) : filepath to data file
        mol_id (str) : identifier for the molecule this data belongs to
    """

    def __init__(self, filepath, mol_id=None, sql=True):
        self.log_path = filepath
        self.mol_id = mol_id
        self.uuid = str(uuid.uuid4())
        self.sql = sql
        self.code_used = "Gaussian"

        self.DFTData = GaussianOutput(filepath)
        self.get_homo_lumo_data()
        self.singlet_excitations = self.get_tddft_excitations(filepath)["Singlet"]

    @property
    def data(self):
        """
        Returns DFT information in a dictionary that matches the No-SQL schema
        """
        data_dict = {
            "code_used": self.code_used,
            "functional": self.DFTData.functional,
            "basis_set": self.DFTData.basis_set,
            "charge": self.DFTData.charge,
            "spin_multiplicity": self.DFTData.spin_multiplicity,
            "scf_total_energy": self.DFTData.final_energy * 27.2114,  # convert to eV
            "homo": self.homo,
            "lumo": self.lumo,
            "first_excitation": [e[0] for e in self.singlet_excitations]  # extract eV energy value
        }
        if self.sql:
            data_dict.update({"calculation_id": self.uuid, "mol_id": self.mol_id})

        json_data = json.dumps(data_dict, default=str)
        return json.loads(json_data)

    def get_homo_lumo_data(self):
        """
        Get homo and lumo energies from a Gaussian molecule
        return
            homo_lumo (dict) : dictionary containing homo then lumo in eV
        """
        num_electrons = self.DFTData.electrons[0]
        eigens = list(self.DFTData.eigenvalues.values())[0]
        self.homo = eigens[num_electrons - 1] * 27.2114  # convert to eV
        self.lumo = eigens[num_electrons] * 27.2114  # convert to eV

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


class ProcessUvVis:
    """
    Class to process UV-Vis data files.
    Copyright 2021, University of Kentucky
    Args:
        filepath (str) : filepath to data file
        mol_id (str) : identifier for the molecule this data belongs to
        metadata (dict) : dictionary containing any metadata for this molecule, e.g., {"solvent": "acetonitrile"}
    """

    def __init__(self, filepath, mol_id, metadata=None, sql=True,):
        self.mol_id = mol_id
        self.filepath = filepath
        self.uuid = str(uuid.uuid4())
        self.sql = sql

        metadata = metadata or {}
        self.instrument = metadata.get("instrument", '')
        self.solvent = metadata.get("solvent", '')

        self.parse_file()

    def parse_file(self):
        """
        Use Pandas to parse the raw data file
        """
        df = pd.read_csv(self.filepath, header=None, names=['col1', 'col2'])
        data_df = df.iloc[4:, :].astype(float, errors='ignore')
        self.data_df = data_df.rename(columns={'col1': 'wavelength', 'col2': 'absorbance'})
        self.string_data = df.iloc[:3, :]

    @property
    def integration_time(self):
        query = self.string_data[self.string_data["col1"].str.contains('Integration Time')]['col2'].values
        return query[0] if query else None

    @property
    def date_recorded(self):
        query = self.string_data[self.string_data["col1"].str.contains('Timestamp')]['col2'].values
        return query[0] if query else ''

    @property
    def raw_absorbance_data(self):
        return self.data_df.to_dict('records')

    @property
    def first_peak(self):
        """
        Get first peak from absorption data
        Returns: float wavelength of the earliest peak
        """
        data = self.raw_absorbance_data
        wavelengths = [x.get('wavelength') for x in data]
        absorbances = [x.get('absorbance') for x in data]
        peaks, _ = find_peaks(absorbances, height=0.3)
        peaks_wavelength = [wavelengths[p]for p in peaks]
        return min(peaks_wavelength)

    @property
    def data(self):
        """
        Returns UV-Vis information in a dictionary that matches the No-SQL or SQL schema
        """
        data_dict = {
            "date_recorded": self.date_recorded,
            "solvent": self.solvent,
            "instrument": self.instrument,
            "integration_time": self.integration_time,
            "optical_gap": round(1240 / self.first_peak, 3),  # convert from nm to eV
        }
        if self.sql:
            data_dict.update({"uvvis_id": self.uuid, "mol_id": self.mol_id})
        else:
            data_dict.update({"absorbance_data": self.raw_absorbance_data})

        json_data = json.dumps(data_dict, default=str)
        return json.loads(json_data)

    @property
    def absorbance_data(self):
        """
        Returns UV-Vis information in a dictionary that matches the SQL AbsorbanceData Table schema
        """
        # .update({"absorbance_id": str(uuid.uuid4()), "uvvis_id": self.uuid, "mol_id": self.mol_id})
        data = self.raw_absorbance_data
        for entry in data:
            entry.update({"absorbance_id": str(uuid.uuid4()), "uvvis_id": self.uuid, "mol_id": self.mol_id})
        return data
