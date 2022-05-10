import uuid
import json
from datetime import date
from databases_demo.processors.uvvis_parser import *
from databases_demo.processors.dft_parser import *
import pubchempy as pcp
from rdkit.Chem.rdchem import Mol
from rdkit.Chem.rdmolops import AddHs
from rdkit.Chem import MolFromSmiles, MolToSmiles, AllChem
from rdkit.Chem.rdMolDescriptors import CalcMolFormula, CalcExactMolWt


class GenerateMolInfo:
    """
        Generate json object for insertion from smiles string
        :param smiles: smiles string
        :param source: which group the molecule comes from
        :param names: list of names for molecule
        Copyright 2021, University of Kentucky
        """
    def __init__(self, smiles, source="", date_made="", names=[]):
        smiles = smiles
        names = names
        # Generate rdkit mol and final (cleaned) smiles
        rdkmol = MolFromSmiles(smiles)
        clean_smile = MolToSmiles(rdkmol)
        rdkmol_hs = AddHs(rdkmol)
        AllChem.EmbedMolecule(rdkmol_hs)
        pcpmol = pcp.get_compounds(clean_smile, namespace="smiles")[0]

        # Populate class
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
    def no_sql_data(self):
        data_dict = {
            "_id": self._id,
            "smiles": self.smiles,
            "synonyms": self.synonyms,
            "molecular_formula": self.molecular_formula,
            "number_of_atoms": self.number_of_atoms,
            "molecular_weight": self.molecular_weight,
            "source": self.source,
            "date_made": self.date_made,
        }
        json_data = json.dumps(data_dict, default=str)
        return json.loads(json_data)

    @property
    def sql_synonym_data(self):
        return [{"synonym": synonym, "mol_id": self._id} for synonym in self.synonyms]

    @property
    def sql_data(self):
        data_dict = {
            "mol_id": self._id,
            "smiles": self.smiles,
            "molecular_formula": self.molecular_formula,
            "number_of_atoms": self.number_of_atoms,
            "molecular_weight": self.molecular_weight,
            "source": self.source,
            "date_made": self.date_made,
        }
        json_data = json.dumps(data_dict, default=str)
        return json.loads(json_data)


class ProcessDFT:
    """
    Class to process DFT logfiles.
    Copyright 2021, University of Kentucky
    """

    def __init__(self, filepath, mol_id=None, parsing_class=ParseGausLog):
        self.log_path = filepath
        self.mol_id = mol_id
        self.uuid = str(uuid.uuid4())

        self.DFTData = parsing_class(filepath)

    @property
    def no_sql_data(self):
        data_dict = {
            "code_used": self.DFTData.code_used,
            "functional": self.DFTData.functional,
            "basis_set": self.DFTData.basis_set,
            "charge": self.DFTData.charge,
            "spin_multiplicity": self.DFTData.spin_multiplicity,
            "scf_total_energy": self.DFTData.scf_total_energy,
            "homo": self.DFTData.homo,
            "lumo": self.DFTData.lumo,
        }
        json_data = json.dumps(data_dict, default=str)
        return json.loads(json_data)

    @property
    def sql_data(self):
        data_dict = {
            "calculation_id": self.uuid,
            "mol_id": self.mol_id,
            "code_used": self.DFTData.code_used,
            "functional": self.DFTData.functional,
            "basis_set": self.DFTData.basis_set,
            "charge": self.DFTData.charge,
            "spin_multiplicity": self.DFTData.spin_multiplicity,
            "scf_total_energy": self.DFTData.scf_total_energy,
            "homo": self.DFTData.homo,
            "lumo": self.DFTData.lumo,
        }
        json_data = json.dumps(data_dict, default=str)
        return json.loads(json_data)


class ProcessUvVis:
    """
    Class to process Gaussian logfiles.
    Copyright 2021, University of Kentucky
    """

    def __init__(self, filepath, mol_id, metadata=None, parsing_class=ParseExcel):
        self.mol_id = mol_id
        self.uuid = str(uuid.uuid4())

        metadata = metadata or {}
        self.instrument = metadata.get("instrument", '')
        self.solvent = metadata.get("solvent", '')

        self.UvVisData = parsing_class(filepath)

    @property
    def no_sql_data(self):
        all_data_dict = {
            "date_recorded": self.UvVisData.date_recorded,
            "solvent": self.solvent,
            "instrument": self.instrument,
            "integration_time": self.UvVisData.integration_time,
            "absorbance_data": self.UvVisData.absorbance_data,
        }
        json_data = json.dumps(all_data_dict, default=str)
        return json.loads(json_data)

    @property
    def sql_data(self):
        all_data_dict = {
            "cv_id": self.uuid,
            "mol_id": self.mol_id,
            "date_recorded": self.CVData.date_recorded,
            "working_electrode": self.working_electrode,
            "counter_electrode": self.counter_electrode,
            "reference_electrode": self.reference_electrode,
            "solvent": self.solvent,
            "electrolyte": self.electrolyte,
            "ionic_liquid": self.ionic_liquid,
            "instrument": self.instrument,
            "working_electrode_surface_area": self.working_electrode_surface_area,
            "redox_mol_concentration": self.redox_mol_concentration,
            "scan_rate": self.CVData.scan_rate,
            "num_scans": self.CVData.num_scans,
            "quiet_time": self.CVData.quiet_time,
            "sensitivity": self.CVData.sensitivity,
            "comp_r": self.CVData.comp_R,
            "scan_data": self.CVData.scan_data,
        }
        json_data = json.dumps(all_data_dict, default=str)
        return json.loads(json_data)
