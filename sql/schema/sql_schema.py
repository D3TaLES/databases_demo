import datetime
from sqlalchemy import (
    Column,
    DateTime,
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
    date_made = Column(DateTime, default=datetime.datetime.utcnow, index=True)


class DftData(Base):
    __tablename__ = "dft_data"

    calculation_id = Column(String, primary_key=True)
    mol_id = Column(String, ForeignKey("molecules.mol_id"))
    date_recorded = Column(DateTime, default=datetime.datetime.utcnow, index=True)
    code_used = Column(String)
    functional = Column(String)
    basis_set = Column(String)
    charge = Column(Integer)
    spin_multiplicity = Column(Integer)
    scf_dipole_moment = Column(Float)
    scf_total_energy = Column(Float)
    homo = Column(Float)
    lumo = Column(Float)


class Geometry(Base):
    __tablename__ = "geometry"

    geometry_id = Column(String, primary_key=True)
    mol_id = Column(String, ForeignKey("molecules.mol_id"))
    # .....


class CvData(Base):
    __tablename__ = "cv_data"

    cv_id = Column(String, primary_key=True)
    mol_id = Column(String, ForeignKey("molecules.mol_id"))
    date_recorded = Column(DateTime, default=datetime.datetime.utcnow, index=True)
    working_electrode = Column(String)
    counter_electrode = Column(String)
    reference_electrode = Column(String)
    solvent = Column(String)
    electrolyte = Column(String)
    ionic_liquid = Column(String)
    instrument = Column(String)
    scan_rate = Column(Float)
    num_scans = Column(Integer)
    initial_potential = Column(Float)
    working_electrode_surface_area = Column(Float)
    redox_mol_concentration = Column(Float)
    quiet_time = Column(Float)
    sensitivity = Column(Float)
    comp_r = Column(Float)


class Synonyms(Base):
    __tablename__ = "synonyms"

    synonym = Column(String, primary_key=True)
    mol_id = Column(String, ForeignKey("molecules.mol_id"))
