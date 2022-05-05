CREATE TABLE [Molecules] (
	'mol_id' VARCHAR(255)  NOT NULL PRIMARY KEY,
	'smiles' VARCHAR(255),
	'molecular_formula' VARCHAR(255),
	'number_of_atoms' INTEGER,
	'molecular_weight' FLOAT,
	'source' VARCHAR(255),
	'date_made' DATE,
);

CREATE TABLE [DftData] (
	'calculation_id' VARCHAR(255)  NOT NULL PRIMARY KEY,
    'mol_id' VARCHAR(255),
	'date_recorded' DATE,
	'code_used' VARCHAR(255),
	'functional' VARCHAR(255),
	'basis_set' VARCHAR(255),
	'charge' INTEGER,
	'spin_multiplicity' INTEGER,
	'scf_dipole_moment' FLOAT,
	'scf_total_energy' FLOAT,
	'homo' FLOAT,
	'lumo' FLOAT,
);

CREATE TABLE [Geometry] (
    'geometry_id' VARCHAR(255)  NOT NULL PRIMARY KEY,
	'mol_id' VARCHAR(255),
    ...
    FOREIGN KEY(mol_id) REFERENCES album(Molecules)
);

CREATE TABLE [CvData] (
	'cv_id' VARCHAR(255)  NOT NULL PRIMARY KEY,
	'mol_id' VARCHAR(255),
    'date_recorded' DATE,
    'working_electrode' VARCHAR(255),
    'counter_electrode' VARCHAR(255),
    'reference_electrode' VARCHAR(255),
    'solvent' VARCHAR(255),
    'electrolyte' VARCHAR(255),
    'ionic_liquid' VARCHAR(255),
    'instrument' VARCHAR(255),
    'scan_rate' FLOAT,
    'num_scans' INTEGER,
    'initial_potential' FLOAT,
    'working_electrode_surface_area' FLOAT,
    'redox_mol_concentration' FLOAT,
    'quiet_time' FLOAT,
    'sensitivity' FLOAT,
    'comp_r' FLOAT,
);

CREATE TABLE [ScanData] (
	'point_id' VARCHAR(255) NOT NULL PRIMARY KEY,
    'cv_id' VARCHAR(255),
    'point_sequence' INTEGER,
    'x_coord' FLOAT,
    'y_coord' FLOAT,
    FOREIGN KEY(cv_id) REFERENCES album(CvData)
);

CREATE TABLE [Synonyms] (
	'synonym' VARCHAR(255) PRIMARY KEY NOT NULL,
	'mol_id' VARCHAR(255),
    FOREIGN KEY(mol_id) REFERENCES album(Molecules)
);