import json
import uuid
from datetime import date

import pubchempy as pcp
from rdkit.Chem import MolFromSmiles, MolToSmiles, AllChem
from rdkit.Chem.rdMolDescriptors import CalcMolFormula, CalcExactMolWt
from rdkit.Chem.rdchem import Mol
from rdkit.Chem.rdmolops import AddHs

from databases_demo.processors.dft_parser import *
from databases_demo.processors.uvvis_parser import *


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
            "synonyms": self.synonyms,
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

    def synonym_data(self):
        """
        Returns synonym information in a dictionary that matches the SQL Synonyms Table schema
        """
        return [{"synonym": synonym, "mol_id": self._id} for synonym in self.synonyms]


class ProcessDFT:
    """
    Class to process DFT logfiles.
    Copyright 2021, University of Kentucky
    Args:
        filepath (str) : filepath to data file
        mol_id (str) : identifier for the molecule this data belongs to
        parsing_class (class) : a DFT parsing class with which to parse the data
    """

    def __init__(self, filepath, mol_id=None, sql=True, parsing_class=ParseGausLog):
        self.log_path = filepath
        self.mol_id = mol_id
        self.uuid = str(uuid.uuid4())
        self.sql = sql

        self.DFTData = parsing_class(filepath)

    @property
    def data(self):
        """
        Returns DFT information in a dictionary that matches the No-SQL schema
        """
        data_dict = {
            "code_used": self.DFTData.code_used,
            "functional": self.DFTData.functional,
            "basis_set": self.DFTData.basis_set,
            "charge": self.DFTData.charge,
            "spin_multiplicity": self.DFTData.spin_multiplicity,
            "scf_total_energy": self.DFTData.scf_total_energy,
            "homo": self.DFTData.homo,
            "lumo": self.DFTData.lumo,
            "homo_lumo_gap": self.DFTData.lumo,
        }
        if self.sql:
            data_dict.update({"calculation_id": self.uuid, "mol_id": self.mol_id})

        json_data = json.dumps(data_dict, default=str)
        return json.loads(json_data)


class ProcessUvVis:
    """
    Class to process UV-Vis data files.
    Copyright 2021, University of Kentucky
    Args:
        filepath (str) : filepath to data file
        mol_id (str) : identifier for the molecule this data belongs to
        metadata (dict) : dictionary containing any metadata for this molecule, e.g., {"solvent": "acetonitrile"}
        parsing_class (class) : a UV-Vis parsing class with which to parse the data
    """

    def __init__(self, filepath, mol_id, metadata=None, sql=True, parsing_class=ParseExcel):
        self.mol_id = mol_id
        self.uuid = str(uuid.uuid4())
        self.sql = sql

        metadata = metadata or {}
        self.instrument = metadata.get("instrument", '')
        self.solvent = metadata.get("solvent", '')

        self.UvVisData = parsing_class(filepath)

    @property
    def data(self):
        """
        Returns UV-Vis information in a dictionary that matches the No-SQL or SQL schema
        """
        data_dict = {
            "date_recorded": self.UvVisData.date_recorded,
            "solvent": self.solvent,
            "instrument": self.instrument,
            "integration_time": self.UvVisData.integration_time,
        }
        if self.sql:
            data_dict.update({"uvvis_id": self.uuid, "mol_id": self.mol_id})
        else:
            data_dict.update({"absorbance_data": self.UvVisData.absorbance_data})

        json_data = json.dumps(data_dict, default=str)
        return json.loads(json_data)

    @property
    def absorbance_data(self):
        """
        Returns UV-Vis information in a dictionary that matches the SQL AbsorbanceData Table schema
        """
        return [data.update({"uvvis_id": self.uuid, "mol_id": self.mol_id}) for data in self.UvVisData.absorbance_data]
