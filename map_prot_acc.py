#!/usr/bin/env python
# coding: utf-8

import pandas as pd
from config import BASE_DIR, DATA_DIR, CODE_DIR
from common_functions import mapping_card_protein_accession
from common_functions import mapping_armfinderplus_protein_accession
from common_functions import mapping_resfinder_protein_accession
 

def main():
    df = mapping_card_protein_accession(card_data_dir=DATA_DIR + '/card/card-data/' )
    df.to_csv(BASE_DIR + '/map_tsv/card_map.tsv', sep='\t', index=False)

    df = mapping_armfinderplus_protein_accession(amrfinderplus_dir=DATA_DIR + '/AMRFinderPlus/')
    df.to_csv(BASE_DIR + '/map_tsv/amrfinderplus_map.tsv', sep='\t', index=False)

    df = mapping_resfinder_protein_accession(resfinder_dir=DATA_DIR + '/resfinder_db/')
    df.to_csv(BASE_DIR + '/map_tsv/resfinder_map.tsv', sep='\t', index=False)

if __name__ == "__main__":
    main()











