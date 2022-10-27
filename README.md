# Database Demo

This repository gives simple, chemistry-based demonstrations of both an SQL
and a No-SQL database structure. For each database structure, there exists 
an example notebook that shows how to: 
1. Initialize the database 
2. Load a schema and use it to validate example data
3. Insert validated data into the database
4. Query the database.

The repository contains five key parts: 

* [`notebooks`](notebooks): Directory containing the example notebooks that demonstrate file processing and using SQL and No-SQL databases. 
* [`schema`](schema): Directory containing the schema for both SQL and No-SQL database examples. 
* [`raw_data`](raw_data): Directory containing the raw data files for data that is inserted into the databases 
in the examples. There are two computational data files (one from the software 
Gaussian and one from the software Psi4) and one experimental potetiostat file. 
* [`file_parser.py`](file_parser.py):  File containing python code for extracting key values from the raw
computational and experimental data files. These parsing codes are used in the 
examples. The [processing_notebook.ipynb](notebooks/processing_demo.ipynb) shows the coding principles behind the code in file_parser.py
* [`external_resources.md`](external_resources.md): A list of external resources that give more specific details for setting up a database.


## Parsing Example
[![Open in Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/D3TaLES/databases_demo/blob/main/notebooks/processing_demo.ipynb)

## SQL Example
[![Open in Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/D3TaLES/databases_demo/blob/main/notebooks/sql_demo.ipynb)

## No-SQL Example
[![Open in Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/D3TaLES/databases_demo/blob/main/notebooks/no_sql_demo.ipynb)

