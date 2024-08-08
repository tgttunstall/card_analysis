#!/usr/bin/env python
# coding: utf-8

import pandas as pd
from config import MAP_DIR, CARD_DIR, AMRFINDERPLUS_DIR, RESFINDER_DIR, POINTFINDER_DIR
from common_functions import mapping_card_protein_accession
from common_functions import mapping_armfinderplus_protein_accession
from common_functions import mapping_resfinder_protein_accession
 

def main():
    df = mapping_card_protein_accession(card_data_dir=CARD_DIR + 'data/' )
    df.to_csv(MAP_DIR / 'card_map.tsv', sep='\t', index=False)

    df = mapping_armfinderplus_protein_accession(amrfinderplus_dir=AMRFINDERPLUS_DIR)
    df.to_csv(MAP_DIR / 'amrfinderplus_map.tsv', sep='\t', index=False)

    df = mapping_resfinder_protein_accession(resfinder_dir=RESFINDER_DIR)
    df.to_csv(MAP_DIR / 'resfinder_map.tsv', sep='\t', index=False)

if __name__ == "__main__":
    main()











