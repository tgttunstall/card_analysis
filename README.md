# README

## Overview

This repository contains scripts for mapping protein accessions and generating knowledge graphs using various antimicrobial resistance databases. The primary functionalities include downloading and processing data from CARD, AMRFinderPlus, and ResFinder databases, mapping protein accessions to UniProtKB accessions, and creating visual knowledge graphs.

## Directory Structure

- `config.py`: Configuration file to define directory paths.
- `map_prot_acc.py`: Script to map protein accessions to UniProtKB accessions.
- `KGs_acc.py`: Script to create knowledge graphs from UniProtKB accessions.
- `common_functions.py`: Contains various helper functions used in the mains scripts.

## Usage

1. **Set up the Conda environment**:
   Create and activate a new Conda environment with the required dependencies using the `environment.yml` file provided.

   ```bash
   conda env create -f environment.yml
   conda activate amr-mapping
   ```

2. **Mapping Protein Accessions**:
   Run `map_prot_acc.py` to map protein accessions to UniProtKB accessions. The results will be saved in TSV files.

   ```bash
   python map_prot_acc.py
   ```

3. **Creating Knowledge Graphs**:
   Run `KGs_acc.py` with the appropriate arguments to create knowledge graphs.

   ```bash
   python KGs_acc.py --accession <UniProtKB_accession> --outdir <output_directory>
   ```

## Dependencies

The project relies on the following Python packages:

- pandas
- requests
- networkx
- pyvis
- obonet
- numpy
- time

---

## File Descriptions

### `config.py`

Defines the directory paths used throughout the scripts.

### `map_prot_acc.py`

Main script to map protein accessions to UniProtKB accessions. It reads the data from the specified directories, performs the mapping, and saves the results in TSV files.

### `KGs_acc.py`

Script to create knowledge graphs from UniProtKB accessions. It uses the mapped data to generate visual knowledge graphs and saves them as HTML files.

### `common_functions.py`

Contains various helper functions used in the main scripts. These include functions for downloading and processing data, parsing JSON and OBO files, generating graphs, and performing ID mappings.