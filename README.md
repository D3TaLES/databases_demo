# Database Demo

This repository gives simple, chemistry-based demonstrations of both an SQL
and a No-SQL database structure. For each database structure, there exists 
an example notebook that shows how to: 
1. Initiate the database 
2. Load a schema and use it to validate example data
3. Insert validated data into the database
4. Query the database.

The repository contains four key directories: 

* `sql`: Contains the schema and all examples for a SQL database. 
* `no-sql`: Contains the schema and all examples for a No-SQL database. 
* `raw_data`: Contains the raw data files for data that is inserted into the databases 
in the examples. There are two computational data files (one from the software 
Gaussian and one from the software Psi4) and one experimental potetiostat file. 
* `processors`:  Contains python code for extracting key values from the raw
computational and experimental data files. These parsing codes are used in the 
examples. 

## SQL Example
[Notebook](https://github.com/D3TaLES/databases_demo/blob/main/sql/sql_notebook.ipynb)

## SQL Example
[Notebook](https://github.com/D3TaLES/databases_demo/blob/main/no-sql/no-sql_notebook.ipynb)

