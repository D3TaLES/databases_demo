# Database Demo

This repository gives simple, chemistry-based demonstrations of both an SQL
and a No-SQL database structure. For each database structure, there exists 
an example notebook that shows how to: 
1. Initialize the database 
2. Load a schema and use it to validate example data
3. Insert validated data into the database
4. Query the database.

The repository contains four key directories: 

* `schema`: Contains the all schema for both SQL and No-SQL database examples. 
* `raw_data`: Contains the raw data files for data that is inserted into the databases 
in the examples. There are two computational data files (one from the software 
Gaussian and one from the software Psi4) and one experimental potetiostat file. 
* `processors`:  Contains python code for extracting key values from the raw
computational and experimental data files. These parsing codes are used in the 
examples. 

## SQL Example
[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/drive/1gQV3LxoQ65NyTFQulHzLTp8IRV3lCICT?usp=sharing)

## No-SQL Example
[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/drive/1LFJUazlB9JYoeqk6U9OVJ3l_m6tu_knD?usp=sharing)

