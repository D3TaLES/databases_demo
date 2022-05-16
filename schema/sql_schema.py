import datetime
from sqlalchemy import (
    Column,
    Float,
    ForeignKey,
    Integer,
    String,
)
from qcfractal.storage_sockets.models.sql_base import Base


class Molecules(Base):
    __tablename__ = "molecules"

    mol_id = Column(String, primary_key=True)
    smiles = Column(String, nullable=False)
    molecular_formula = Column(String, nullable=False)
    number_of_atoms = Column(Integer, nullable=False)
    molecular_weight = Column(Float, nullable=False)
    source = Column(String)
    date_made = Column(String, default=datetime.datetime.utcnow, index=True)


class DftData(Base):
    __tablename__ = "dft_data"

    calculation_id = Column(String, primary_key=True)
    mol_id = Column(String, ForeignKey("molecules.mol_id"))
    date_added = Column(String, default=datetime.datetime.utcnow, index=True)
    code_used = Column(String)
    functional = Column(String)
    basis_set = Column(String)
    charge = Column(Integer)
    spin_multiplicity = Column(Integer)
    scf_total_energy = Column(Float)
    homo = Column(Float)
    lumo = Column(Float)
    homo_lumo_gap = Column(Float)


class UvVisData(Base):
    __tablename__ = "uvvis_data"

    uvvis_id = Column(String, primary_key=True)
    mol_id = Column(String, ForeignKey("molecules.mol_id"))
    date_recorded = Column(String)
    solvent = Column(String)
    instrument = Column(String)
    integration_time = Column(Float)
    optical_gap = Column(Float)


class AbsorbanceData(Base):
    __tablename__ = "absorbance_data"

    absorbance_id = Column(String, primary_key=True)
    uvvis_id = Column(String, ForeignKey("uvvis_data.uvvis_id"))
    mol_id = Column(String, ForeignKey("molecules.mol_id"))
    wavelength = Column(Float)
    absorbance = Column(Float)


class Synonyms(Base):
    __tablename__ = "synonyms"

    synonym = Column(String, primary_key=True)
    mol_id = Column(String, ForeignKey("molecules.mol_id"))
