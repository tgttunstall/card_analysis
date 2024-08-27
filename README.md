# README

## Overview 

This repository provides a comprehensive pipeline for mapping protein accessions and generating knowledge graphs using several antimicrobial resistance (AMR) databases. The project is designed to facilitate the exploration of relationships between protein sequences and AMR mechanisms by leveraging data from CARD, AMRFinderPlus, ResFinder.

![KG example](https://github.com/ebi-uniprot/card-analisys/blob/main/example/Q9K351.png)

## Features

- **Automated Database Download**: Scripts for downloading and setting up required AMR databases.
- **Protein Accession Mapping**: Tools to map various protein accessions to UniProtKB accessions.
- **Knowledge Graph Generation**: Scripts to generate visual knowledge graphs in multiple formats (e.g., HTML, PNG, Cytoscape JSON) from UniProtKB accessions.

## Directory Structure

```plaintext
bioinformatics-amr
├── code/                     # Directory for custom Python scripts
├── databases/                # Directory for downloaded AMR databases
│   ├── card/
│   │   ├── data/
│   │   ├── ontology/
│   ├── amrfinderplus/
│   ├── resfinder/
│   ├── pointfinder/
├── map_tsv/                  # Directory for mapped TSV files
├── config.py                 # Configuration file for directory paths
├── download_databases.py     # Script to download and extract AMR databases
├── map_prot_acc.py           # Script for mapping protein accessions to UniProtKB
├── KGs_acc.py                # Script for generating knowledge graphs
├── common_functions.py       # Script with helper functions
├── environment.yml           # Conda environment configuration
└── README.md                 # This file
```
## Usage

### 1. Set up Conda Environment

Ensure you have [Miniconda](https://docs.conda.io/en/latest/miniconda.html) installed, which is a minimal installer for Conda. Create and activate the Conda environment using the provided `environment.yml` file.

```bash
conda env create -n bioinformatics-amr -f environment.yml
conda activate bioinformatics-amr
```

### 2. Download Required Databases

Run the `download_databases.py` script to download and set up the required databases.

```bash
python download_databases.py
```

### 3. Map Protein Accessions

Use `map_prot_acc.py` to map protein accessions to UniProtKB accessions.

```bash
python map_prot_acc.py
```

### 4. Generate Knowledge Graphs

Generate knowledge graphs using the `KGs_acc.py` script. Specify the UniProtKB accession and output directory.

```bash
python KGs_acc.py --accession <UniProtKB_accession> --outdir <output_directory> --visualization <pyvis|cytoscape|png>
```

## Detailed Descriptions

### `config.py`

Defines the directory paths used across the project. It ensures that all scripts reference consistent directories for data storage and retrieval.

### `download_databases.py`

Automates the download and extraction of AMR databases, including CARD, AMRFinderPlus, ResFinder, and PointFinder. The script handles both HTTP and FTP downloads, as well as Git repository cloning.

### `map_prot_acc.py`

Maps protein accessions from the downloaded AMR databases to UniProtKB accessions using various mapping strategies, including API-based completions for unmatched entries.

### `KGs_acc.py`

Generates knowledge graphs based on UniProtKB accessions using data from the CARD, AMRFinderPlus, and ResFinder databases. Supports output in various formats, including HTML for PyVis visualizations and JSON for Cytoscape.

### `common_functions.py`

A collection of utility functions for downloading files, extracting data, mapping IDs, and generating graphs. It includes functions for interacting with FTP servers, cloning Git repositories, parsing OBO files, and creating visualizations.

### `environment.yml`

Specifies the Conda environment with all required dependencies, including packages for data processing (e.g., `pandas`, `numpy`), network visualization (e.g., `networkx`, `pyvis`), and web requests (e.g., `requests`).

## Dependencies

Key dependencies for this project include:

- `pandas`
- `networkx`
- `pyvis`
- `obonet`
- `numpy`
- `requests`
- `gitpython`
- `ftplib`
- `matplotlib`

For the full list of dependencies, refer to the `environment.yml` file.


